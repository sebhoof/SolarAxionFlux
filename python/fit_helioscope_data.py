from mpi4py import MPI
import os
import sys
import time
import numpy as np

from datetime import date
from scipy.optimize import minimize
from iminuit import Minuit

os.environ["OMP_NUM_THREADS"] = "1"

gagg0 = 0.6 # = 0.6e-10 GeV^{-1}
if len(sys.argv) > 1:
    gagg0 = float(sys.argv[1])
code_path = "."
output_path = "."
if len(sys.argv) > 3:
    code_path = str(sys.argv[2])
    output_path = str(sys.argv[3])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
# Hacky way for safe read/write of solar model files etc.
time.sleep(rank)

sys.path.append(code_path)
from python.physics import *
from python.tomography.tomography import *
from python.tomography.fitting import *

sys.path.append(code_path+"lib")
from pyaxionflux import SolarModel

model_file = code_path+"/data/solar_models/SolarModel_B16-AGSS09.dat"
sol = SolarModel(model_file)

# Set up the MPI environment and variables.
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

rbins0 = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0])
n_rbins0 = len(rbins0)-1
ndim = 2*n_rbins0+1
n_ebins0 = 4
buffer_size = ndim+(n_ebins0+1)+2
grid_bins0 = 128
dims0 = 2*[grid_bins0]
spot0 = [65, 68, 56]
matrix0 = create_matrix_poly(rbins0)
tru_t = sol.temperature(rbins0[:-1])
tru_k = np.sqrt(sol.kappa_squared(rbins0[:-1]))
tru_x = np.array([gagg0] + list(tru_k) + list(tru_t))

xb0 = [[0.75*gagg0, 1.25*gagg0]] + [[0.01, max(k+5, 3*k)] for k in tru_k] + 2*[[0.5, 2.0]] + [[0.01, max(t+1.0, 1.5*t)] for t in tru_t[2:]]
# Not used anymore:
# de_options = { 'popsize': 100, 'tol': 1e-4, 'polish': False, 'maxiter': 1000, 'init': 'sobol', 'x0': tru_x, 'bounds': xb0 }

flux_P_data = np.genfromtxt(code_path+"/results/flux_P_B16-AGSS09.dat")
flux_P_interp = CubicSpline(flux_P_data[:,0], flux_P_data[:,1])
n_counts_code0 = gagg0*gagg0*flux_to_events(flux_P_interp.integrate(0.3, 15.0), gagg=gagg0)
n_mc_sims = int(n_counts_code0)
mcg = MCGenerator(code_path+"/results/matrix_P_B16-AGSS09.dat")


def task_wrapper(id):
   t0 = time.time()
   events0 = mcg.draw_events(n_mc_sims)
   grids0, ebins0, *_ = grids_and_binning_from_events(events0, dims0, spot0, n_ebins=n_ebins0, n_rbins=n_rbins0, kind='even_bins')
   counts0 = [get_counts_in_rings_simple(g, spot0, rbins0) for g in grids0]
   pm0 = primakoff_moments_integral(ebins0)

   cost = lambda x: fitting_metric_poly(x, rbins0, matrix0, pm0, counts0)

   confusion = np.random.normal(1, 0.02, size=ndim)
   m1 = minimize(cost, confusion*tru_x, method='Nelder-Mead', bounds=xb0, tol=1e-11, options={ 'adaptive': True, 'maxfev': 2e6 })
   m2 = Minuit(cost, tuple(m1.x))
   m2.tol = 1e-8
   m2.limits = [tuple(x) for x in xb0]
   m2.migrad(ncall=int(1e6))
   
   print("Report for task", id, ":\nMinuit? {} Call limit? {} ({:.2f}), NM? {} ({:.2f}).".format(m2.fmin.is_valid, m2.fmin.has_reached_call_limit, m2.fval, m1.success, m1.fun))

   t1 = time.time()
   dt = (t1-t0)/60.0
   return np.array(list(m2.values)+list(ebins0)+[np.sum(counts0), m2.fval]), dt

def run_mpi_job(nsamples, out_file_name_root, save_temp_results=True):
    # Set up the MPI environment and variables.
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()
    start_time = time.time()

    if (rank > 0):
        # Send first result to main process to signal that this proc is ready for more tasks
        res, dt = task_wrapper(rank-1)
        print('Initial job finished on rank {}, took {:.2} mins.'.format(rank, dt), flush=True)
        comm.Send(res, dest=0)
        for _ in range(nsamples):
            task_id = comm.recv(source=0)
            if (task_id > nsamples):
                break
            res, dt = task_wrapper(task_id-1)
            comm.Send(res, dest=0)
        dt = (time.time()-start_time)/3600.0
        print('MPI rank {} finished! MC simulations took {:.3f} hrs.'.format(rank, dt))
        
    # Main process distribute tasks and receive results
    if (rank == 0):
        all_results = []
        out_file_name = out_file_name_root+".dat"
        out_file_temp_name = out_file_name_root+"_temp.dat"
        n_temp_save = int(nsamples//25 + 1) # Save after ~4% progress for many tasks (otherwise every task)
        print('Main process waiting for {} results from {} other processes...'.format(nsamples, ncores-1), flush=True)

        # Receive results and send more work (send out ncores work task too many to finalise jobs).
        for task_id in range(ncores, nsamples+ncores):
            info = MPI.Status()
            res = np.zeros(buffer_size) # Container for the result
            comm.Recv(res, source=MPI.ANY_SOURCE, status=info)
            worker_id = info.Get_source()
            all_results.append(list(res))
            if save_temp_results and (task_id%n_temp_save == 0):
                t0 = time.time()
                a = np.array(all_results)
                np.savetxt(out_file_temp_name, a, fmt="%.6e")
                t1 = time.time()
                dt = t1 - t0
                info = "Job {}/{} finished on rank {}.".format(task_id-ncores+1, nsamples, worker_id)
                if dt > 5:
                    info += " Writing temporary results took {:.1f} s.".format(dt)
                print(info, flush=True)
            comm.send(task_id, dest=worker_id)
        print('All MPI tasks finished after {:.3f} hrs!'.format( (time.time()-start_time)/3600.0 ), flush=True)
        print('Formatting results and saving them to '+out_file_name+'...')
        a = np.array(all_results)
        header = "Fitting results for g_agamma, kappa_s, and T calculated on {}.\n".format(date.today().strftime("%Y-%m-%d"))
        header += "Columns: bfit"
        np.savetxt(out_file_name, a, fmt="%.6e", header=header)
        print('All tasks complete! Finishing MPI routine...', flush=True)


if __name__ == "__main__":
    run_mpi_job(1000, output_path+"/fit_grids_{:.2f}".format(gagg0), save_temp_results=True)
