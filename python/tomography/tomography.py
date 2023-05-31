import sys, os
import numpy as np

from typing import Union
from scipy.interpolate import interp1d, PPoly, CubicSpline, PchipInterpolator
from scipy.stats import poisson, uniform
from scipy.integrate import quad
from scipy.linalg import inv

from ..constants import *
from ..physics import *
from ..grid_integrator.flux_grid_integrator import *


###################################
#  Monte Carlo event generation   #
###################################

class MCGenerator:
    """
    This is a Monte Carlo generator class to generate detector events from solar axions.
    """
    def __init__(self, filename: str) -> None:
        """
        Initialize the MCGenerator with a file for the differential flux on the solar disc.

        Args:
        filename (str): Path to the file containing the differential flux.

        Returns:
        None
        """
        # Read in a file for the differential flux on the solar disc and interpolate it
        data = np.genfromtxt(filename)
        rads = np.sort(np.unique(data[:,0]))
        ergs = np.sort(np.unique(data[:,1]))
        e_min, e_max = ergs[0], ergs[-1]
        self.diff = rads[1] - rads[0]
        rflux, efluxes_interp = [], []
        # Interpolate and integrate the differential flux for each radius
        for r in rads:
            eflux_interp = CubicSpline(ergs, data[data[:,0]==r])
            efluxes_interp.append(eflux_interp)
            rflux.append(eflux_interp.integrate(e_min, e_max)[-1])
        # Generate an interpolating function for the radial flux
        self.rflux_interp = CubicSpline(rads, rflux)
        # Calculate the (conditional) CDFs and invert them by swapping rows and columns
        # First, calculate the CDFs for the radial flux
        rflux_cdf = [self.rflux_interp.integrate(rads[0], r) for r in rads]
        # Find the index of the maximum value of the CDFs, and truncate the arrays at that index
        for i_max in range(len(rflux_cdf)-1):
            if rflux_cdf[i_max+1] <= rflux_cdf[i_max]:
                break
        self.n_expected_ref = rflux_cdf[-1]
        # Normalise the CDFs and generate an interpolating function for them
        rflux_cdf = [cdf/rflux_cdf[-1] for cdf in rflux_cdf[:(i_max+1)]]
        self.rflux_cdf_interp = CubicSpline(rflux_cdf, rads[:(i_max+1)])
        self.efluxes_cdf_interp = []
        # Next, calculate the CDFs for the differential flux for each radius
        for i in range(len(rads)-1):
            eflux_cdf = [efluxes_interp[i+1].integrate(e_min, e)[-1] for e in ergs]
            # Again, find the index of the maximum value of the CDFs, and truncate the arrays at that index
            for i_max in range(len(eflux_cdf)-1):
                if eflux_cdf[i_max+1] <= eflux_cdf[i_max]:
                    break
            # Normalise the CDFs and generate an interpolating function for them
            eflux_cdf = [cdf/eflux_cdf[-1] for cdf in eflux_cdf[:(i_max+1)]]
            self.efluxes_cdf_interp.append(CubicSpline(eflux_cdf, ergs[:(i_max+1)]))

    def draw_r(self, u: np.ndarray) -> np.ndarray:
        """
        Draw a random radius rho.

        Args:
            u (np.ndarray): An array of uniform random numbers between 0 and 1.

        Returns:
            np.ndarray: An array of radial values corresponding to the input array of uniform random numbers.
        """
        return self.rflux_cdf_interp(u)

    def draw_e_given_r(self, u: Union[float, list[float]], r: Union[float, list[float]]) -> Union[float, list[float]]:
        """
        Draw a random axion energy E_a given a radius rho.

        Args:
        - u (float or List[float]): A uniform random variable, or a list of uniform random variables, in the range [0,1].
        - r (float or List[float]): A radius, or a list of radii, in the range [0,1].

        Returns:
        - float or List[float]: A random axion energy E_a, or a list of random axion energies, in units of MeV.
        """
        try:
            i = int(r/self.diff)
            res = 0.5*(self.efluxes_cdf_interp[i](u) + self.efluxes_cdf_interp[i+1](u))
        except TypeError:
            i = [int(rr/self.diff) for rr in r]
            res = [0.5*(self.efluxes_cdf_interp[ii](uu) + self.efluxes_cdf_interp[ii+1](uu)) for uu,ii in zip(u,i)]
        return res

    def draw_events(self, n_expected: int, poisson_draw: bool = True) -> np.ndarray:
        """
        Generate sample of detector events based on n_expected number of events.

        Args:
        - n_expected (int): The number of detector events to generate.
        - poisson_draw (bool): Generate a random number?

        Returns:
        - np.ndarray: An array of shape (n_samples, 3) where each row is a tuple (rho, phi, E_a), representing the radius, azimuthal angle, and axion energy of a detector event.
        """
        if poisson_draw:
            n_samples = poisson.rvs(n_expected)
        else:
            n_samples = n_expected
        rads = uniform.rvs(size=n_samples)
        rads = self.draw_r(rads)
        phis = uniform.rvs(scale=2.0*np.pi, size=n_samples)
        ergs = uniform.rvs(size=n_samples)
        ergs = self.draw_e_given_r(ergs, rads)
        return np.array([rads, phis, ergs]).T

def grids_from_events(events: np.ndarray, dims: tuple[int, int], spot: tuple[float, float, float], ebins: np.ndarray) -> list[np.ndarray]:
    """
    Bins the generated events into bins of fixed size.

    Args:
        events: Numpy array of shape (n_events, 3) containing event coordinates and energies.
        dims: Tuple of (height, width) specifying the dimensions of the grids.
        spot: Tuple of (x, y, r) specifying the center and radius of the spot.
        ebins: Numpy array of shape (n_ebins+1,) containing the energy bin edges.

    Returns:
        A list of numpy arrays of shape (height, width) containing the event counts for each energy bin.

    """
    xc, yc, rc = spot
    n_ebins = len(ebins) - 1
    # Initialize numpy arrays for each energy bin
    grids = [np.zeros(dims, dtype=np.int32) for _ in range(n_ebins)]
    cond = (events[:, 2] > ebins[0])&(events[:, 2] < ebins[-1])
    events_cut = events[cond]
    n_cut = len(events_cut)

    # Do the binning
    indices = np.searchsorted(ebins, events_cut[:, 2], side='left')
    jx = np.floor(xc + rc * events_cut[:, 0] * np.cos(events_cut[:, 1])).astype(np.int32)
    jy = np.floor(yc + rc * events_cut[:, 0] * np.sin(events_cut[:, 1])).astype(np.int32)
    for i in range(n_cut):
        # Add an event to the corresponding energy bin
        grids[indices[i]-1][jx[i], jy[i]] += 1
    return grids

def grids_and_binning_from_events(events: np.ndarray, dims: tuple[int, int], spot: tuple[float, float, float], 
                                   n_ebins: int = 3, n_rbins: int = 7, kind: str = 'even_counts') -> tuple[list[np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    """
    Bin the generated events into energy and radius bins, and calculate the center of the spot.

    Args:
        events: Numpy array of shape (n_events, 3) containing event coordinates and energies.
        dims: Tuple of (height, width) specifying the dimensions of the grids.
        spot: Tuple of (x, y, r) specifying the center and radius of the spot.
        n_ebins: The number of energy bins to generate.
        n_rbins: The number of radius bins to generate.
        kind: If 'even_counts', even out the counts in the radius bins. If 'even_bins', generate evenly spaced radial bins.
              If 'mixed', first even out counts for n_rbins - 3 bins, then add three more bins in the last bin

    Returns:
        A tuple of the following:
        - A list of numpy arrays of shape (height, width) containing the event counts for each energy bin.
        - A numpy array of shape (n_ebins+1,) containing the energy bin edges.
        - A numpy array of shape (n_rbins+1,) containing the radius bin edges.
        - A numpy array of shape (2,) containing the x and y coordinates of the estimated centre of the spot.
    """
    # Get erg bins
    ebins = np.quantile(events[:,2], [i/n_ebins for i in range(1,n_ebins)], method='midpoint')
    ebins = np.insert(ebins, [0, n_ebins-1], [0.3, 15.0])
    # Calculate grids and get rbins
    grids = grids_from_events(events, dims, spot, ebins)
    summed_grids = sum(grids)
    if kind == 'even_counts':
        rbins, *_ = get_counts_in_rings(summed_grids, spot, n_rbins)
    elif kind == 'even_bins':
        rbins = np.linspace(0, 1, n_rbins+1)
    elif kind == 'mixed':
        if n_rbins < 3:
            raise ValueError("Argument 'n_rbins' must be >= 3 for kind='mixed'.")
        rbins, *_ = get_counts_in_rings(summed_grids, spot, n_rbins-2)
        add_rbins = np.linspace(rbins[-2], 1, 4)
        rbins = np.insert(rbins, -1, add_rbins[1:-1])
    else:
        raise ValueError("Argument 'kind' must be one of 'even_counts', 'even_bins', 'mixed'.")
    # Also estimate centre of the circular signal region
    weighted_positions = [
        w * np.array([i + 0.5, j + 0.5])
        for i, ai in enumerate(summed_grids)
        for j, w in enumerate(ai)
    ]
    centre = sum(weighted_positions) / summed_grids.sum()
    return grids, ebins, rbins, centre

#######################################
#  Reconstruction of the Gamma bars   #
#######################################

def dthreehalf(r1: float, r2: float) -> float:
    """
    Calculates the power 3/2 of the difference between two numbers.

    Args:
        r1 (float): The first number.
        r2 (float): The second number.

    Returns:
        float: The power 3/2 of the difference between the two numbers.
    """
    arg = r1*r1 - r2*r2
    return pow(arg, 1.5)

def reconstruction_algorithm(events: np.ndarray, rhos: np.ndarray, exp_settings: dict = iaxo) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Reconstructs the (piecewise-constant) Gamma bars from explicit matrix inversion.

    Args:
        events (np.ndarray): The event counts.
        rhos (np.ndarray): The array of rho values.
        exp_settings (dict): The experimental settings used for the simulation. Defaults to iaxo.

    Returns:
        tuple: A tuple containing the rhos array, Gamma bars and the corresponding errors.
    """
    mat_rescale_factor = flux_to_events(1, gagg=1, exp_settings=exp_settings)*distance_factor
    n_rhobins = len(events)

    # Set up the matrix
    m = np.zeros((n_rhobins, n_rhobins))
    for i in range(n_rhobins):
        ris = rhos[i]*rhos[i]
        for j in range(n_rhobins):
            if i < j:
                rip1s = rhos[i+1]*rhos[i+1]
                rjs = rhos[j]*rhos[j]
                rjp1s = rhos[j+1]*rhos[j+1]
                m[i,j] = - ( (rjp1s - rip1s)**1.5 - (rjp1s - ris)**1.5 - (rjs - rip1s)**1.5 + (rjs - ris)**1.5 ) / 3.0
            elif i==j:
                rip1s = rhos[i+1]*rhos[i+1]
                m[i,j] = pow(rip1s - ris, 1.5) / 3.0
    m *= mat_rescale_factor

    # Solve set of linear equations to recover Gamma
    gbs, gbs_err = np.zeros(n_rhobins), np.zeros(n_rhobins)
    
    for n in range(n_rhobins):
        i = n_rhobins-n-1
        if (n==0):
            gbs[i] = events[i] / m[i,i]
            gbs_err[i] = np.sqrt(events[i]) / m[i,i]
        else:
            gbs[i] = (events[i] - np.dot(m[i,i+1:], gbs[i+1:])) / m[i,i]
            gbs_err[i] = np.sqrt(events[i] + np.dot(m[i,i+1:]*m[i,i+1:], gbs_err[i+1:]*gbs_err[i+1:])) / m[i,i]

    return rhos, gbs, gbs_err

def reconstruct_Gamma_bars(grids: list[np.ndarray], spot: tuple[float, float, float], rbins: np.ndarray, exp_settings: dict = iaxo) -> tuple[np.ndarray, np.ndarray]:
    """
    Reconstructs the (piecewise-constant) Gamma bars from explicit matrix inversion for multiple grids.

    Args:
        grids (list[np.ndarray]): The list of grids.
        spot (tuple[float, float]): The (x,y) position of the spot.
        rbins (np.ndarray): The array of r values.
        exp_settings (dict): The experimental settings used for the simulation. Defaults to iaxo.

    Returns:
        tuple: A tuple containing the Gamma bars array and the corresponding errors.
    """
    gbs, gbs_err = np.array([rbins[:-1]]), np.array([rbins[:-1]])
    for g in grids:
        counts = get_counts_in_rings_simple(g, spot, rbins)
        _, gb, gb_err = reconstruction_algorithm(counts, rbins, exp_settings)
        gbs = np.concatenate((gbs, [gb]))
        gbs_err = np.concatenate((gbs_err, [gb_err]))
    return gbs.T, gbs_err.T

# Reconstruction algorithm with constant functions; numerical matrix inversion
def reconstruction_algorithm_cnst(events: np.ndarray, rbins: np.ndarray, exp_settings: dict = iaxo) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Reconstructs the Gamma function using the constant function approach with numerical matrix inversion.

    Args:
        events: An array of event counts in each radial bin.
        rbins: An array of bin edges for the radial bins.
        exp_settings: A dictionary containing experiment-specific settings.

    Returns:
        A tuple containing the following:
        - rbins: An array of bin edges for the radial bins.
        - ys: An array of the reconstructed Gamma values.
        - cov_ys: An array of the covariance matrix for the reconstructed Gamma values.
        - cov_ys_inv: An array of the inverse covariance matrix for the reconstructed Gamma values.
    """
    rescale_factor = distance_factor*flux_to_events(1.0, gagg=1, exp_settings=exp_settings)
    n = len(rbins)
    m = np.zeros((n,n))

    # Set up the matrix
    for i in range(n-1):
        r1 = rbins[i]
        r2 = rbins[i+1]
        m[i,i] = dthreehalf(r2, r1)/3.0
        for j in range(i+1,n-1):
            r3 = rbins[j]
            r4 = rbins[j+1]
            m[i,j] = (dthreehalf(r4, r1) - dthreehalf(r4, r2) + dthreehalf(r3, r2) - dthreehalf(r3, r1))/3.0
    # Set y_n = 0
    m[n-1,n-1] = 1

    # Solve set of linear equations to recover Gamma from ms and ys
    b = np.zeros(n)
    b[:(n-1)] = events
    db2 = np.zeros((n,n))
    db2[:(n-1),:(n-1)] = np.diag(events)
    idb2 = np.zeros((n,n))
    idb2[:(n-1),:(n-1)] = np.diag(1/events)
    
    rsf2 = rescale_factor*rescale_factor
    b /= rescale_factor
    db2 /= rsf2
    idb2 *= rsf2

    ys = np.linalg.solve(m, b)
    assert np.allclose(m@ys, b)
    m_inv = inv(m)
    cov_ys = m_inv@db2@(m_inv.T)
    cov_ys_inv = (m.T)@idb2@m
    
    return rbins, ys, cov_ys, cov_ys_inv

# Wrapper routine for reconstructing the Gamma bars using the algo above
def reconstruct_Gamma_bars_cnst(grids, spot, rbins):
    gbs, gbs_err, gbs_icvo = [], [], []
    for g in grids:
        counts = get_counts_in_rings_simple(g, spot, rbins)
        _, y, cov_y, cov_y_inv = reconstruction_algorithm_cnst(counts, rbins)
        gbs.append(y)
        gbs_err.append(np.sqrt(np.diag(cov_y)))
        gbs_icvo.append(cov_y_inv)
    return gbs, gbs_err, gbs_icvo


def logterm(r1, r2, r3):
    if (r3 > 0):
        r1s = r1*r1
        r2s = r2*r2
        r3s = r3*r3
        return np.log((r1 - np.sqrt(r1s - r3s))/(r2 - np.sqrt(r2s - r3s)))
    else:
        return 0

# Reconstruction algorithm with linear functions; numerical matrix inversion
def reconstruction_algorithm_linear(events, rbins, exp_settings=iaxo):
    rescale_factor = distance_factor*flux_to_events(1.0, gagg=1, exp_settings=exp_settings)
    n = len(rbins)
    m = np.zeros((n,n))

    # Set up the matrix
    for i in range(n-1):
        r1 = rbins[i]
        r1s = r1*r1
        r2 = rbins[i+1]
        r2s = r2*r2
        rt21 = np.sqrt(r2s-r1s)
        h21 = r2 - r1
        th21 = dthreehalf(r2, r1)
        if (i==0):
            r2c = r2*r2s
            m[i,i] = r2c/12.0
            m[i,i+1] = r2c/4.0
        else:
            mii = (2*r2*th21 + 3*r1s*(r1s*np.log((r2 + rt21)/r1) - r2*rt21))/(24*h21)
            m[i,i] = mii
            m[i,i+1] = th21/3.0 - mii
        for j in range(i+1,n-1):
            r3 = rbins[j]
            r3s = r3*r3
            r4 = rbins[j+1]
            r4s = r4*r4
            h43 = r4 - r3
            th41 = dthreehalf(r4,r1)
            th42 = dthreehalf(r4,r2)
            th31 = dthreehalf(r3,r1)
            th32 = dthreehalf(r3,r2)
            rt41 = np.sqrt(r4s - r1s)
            rt42 = np.sqrt(r4s - r2s)
            rt31 = np.sqrt(r3s - r1s)
            rt32 = np.sqrt(r3s - r2s)
            lt1 = r1s*r1s*logterm(r4, r3, r1)
            lt2 = logterm(r3, r4, r2)
            ltt = r1s*r1s*lt1 + r2s*r2s*lt2
            m[i,j] += (2*r4*(th41 - th42) - (8*r4-3*r3)*(th31 - th32) + 3*(r4*(r2s*rt42-r1s*rt41) + r3*r3s*(rt31-rt32) - ltt))/(24*h43)
            m[i,j+1] += (2*r3*(th31 - th32) - (8*r3-3*r4)*(th41 - th42) + 3*(r3*(r2s*rt32-r1s*rt31) + r4*r4s*(rt41-rt42) + ltt))/(24*h43)
    # Set y_n = 0
    m[n-1,n-1] = 1
    m *= rescale_factor

    b = np.zeros(n)
    b[:(n-1)] = events

    # Solve set of linear equations to recover Gamma from ms and ys
    ys = np.linalg.solve(m, b)
    #print(np.allclose(np.dot(m, ys), b))
    #assert np.allclose(np.dot(m, ys), b)
    
    return rbins, ys

# Wrapper routine for reconstructing the Gamma bars using the algo above
def reconstruct_Gamma_bars_linear(grids, spot, rbins):
    gbs = []
    ys = []
    for g in grids:
        counts = get_counts_in_rings_simple(g, spot, rbins)
        x, y = reconstruction_algorithm_linear(counts, rbins)
        p = interp1d(x, y, kind='linear')
        gbs.append(p)
        ys.append(y)
    return gbs, ys

def create_matrix_spline(rbins, bc='clamped'):
    n = len(rbins)
    m = np.zeros((2*n,2*n))

    # Set up the matrix
    # Set up the upper half (integral conditions; match observed counts)
    for i in range(n-1):
        r1 = rbins[i]
        r1s = r1*r1
        r1q = r1s*r1s
        r2 = rbins[i+1]
        r2s = r2*r2
        r2q = r2s*r2s
        rt21 = np.sqrt(r2s - r1s)
        if (i==0):
            r2c = r2*r2s
            r2quin = r2c*r2s
            m[i,i] = r2c/12.0
            m[i,i+1] = r2c/4.0
            m[i,n+i] = -r2quin/90.0
            m[i,n+i+1] = -r2quin/72.0
        else:
            h21 = r2-r1
            l21 = r1s*np.log((r2 + rt21)/r1)
            th21 = dthreehalf(r2, r1)
            mii = (2*r2*th21 + 3*r1s*(l21 - r2*rt21))/(24*h21)
            m[i,i] = mii
            m[i,i+1] = th21/3.0 - mii
            mini = (2*r2*th21*(23*r1s + 20*r1*r2 - 8*r2s) - 15*r1s*(r1s - 4*r1*r2 - 4*r2s)*(l21 - r2*rt21))/(1440*h21)
            m[i,n+i] = mini
            temp = (r1+r2)*(2*th21*(8*r1 - 3*r2) + 15*r1s*(l21 - r2*rt21))/240.0
            m[i,n+i+1] = temp - mini
        for j in range(i+1,n-1):
            r3 = rbins[j]
            r3s = r3*r3
            r4 = rbins[j+1]
            r4s = r4*r4
            h43 = r4 - r3
            rt41 = np.sqrt(r4s - r1s)
            rt42 = np.sqrt(r4s - r2s)
            rt31 = np.sqrt(r3s - r1s)
            rt32 = np.sqrt(r3s - r2s)
            th41 = dthreehalf(r4,r1)
            th42 = dthreehalf(r4,r2)
            th31 = dthreehalf(r3,r1)
            th32 = dthreehalf(r3,r2)
            lt1 = logterm(r4, r3, r1)
            lt2 = logterm(r3, r4, r2)
            ltt = r1q*lt1 + r2q*lt2
            temp = (th32 - th31 + th41 - th42)/3.0
            mij = (2*((3*r3-4*r4)*(th31 - th32) + r4*(th41 - th42)) + 3*(r1s*(r3*rt31 - r4*rt41) + r2s*(r4*rt42 - r3*rt32) - ltt))/(24*h43)
            m[i,j] += mij
            m[i,j+1] += temp - mij
            aux1 = 4*r4s*(5*r3s - 10*r3*r4 + 4*r4s)
            aux2 = 2*(25*r3s - 50*r3*r4 + 11*r4s)
            term1 = r4*(rt41*(-81*r1q - aux1 + r1s*aux2) + rt42*(81*r2q + aux1 - r2s*aux2))
            aux3 = 4*r3*r3s*(5*r3s -14*r3*r4 + 10*r4s)
            aux4 = 4*r3*(5*r3s - 23*r3*r4 + 25*r4s)
            term2 = rt31*(r1q*(-15*r3 + 96*r4) - aux3 + r1s*aux4) + rt32*(3*r2q*(5*r3 - 32*r4) + aux3 - r2s*aux4)
            aux5 = -2*r3s + 4*r3*r4 + 4*r4s
            term3 = 15*(r1q*(r1s + aux5)*lt1 + r2q*(r2s + aux5)*lt2)
            m[i,n+j] += (term1 + term2 - term3)/(1440*h43)
            aux1 = 4*r4*r4s*(10*r3s - 14*r3*r4 + 5*r4s)
            aux2 = 4*r4*(25*r3s - 23*r3*r4 + 5*r4s)
            term1 = rt41*(3*r1q*(32*r3 - 5*r4) + r1s*aux2 - aux1) + rt42*(r2q*(-96*r3 + 15*r4) - r2s*aux2 + aux1)
            aux3 = 4*r3*r3s*(4*r3s - 10*r3*r4 + 5*r4s)
            aux4 = 2*(11*r3s - 50*r3*r4 + 25*r4s)
            term2 = r3*(rt31*(-81*r1q - aux3 + r1s*aux4) + rt32*(81*r2q + aux3 - r2s*aux4))
            aux5 = 4*r3s + 4*r3*r4 - 2*r4s
            term3 = 15*(r1q*(r1s + aux5)*lt1 + r2q*(r2s + aux5)*lt2)
            m[i,n+j+1] += (term1 + term2 + term3)/(1440*h43)
    # Set y_n = 0
    m[n-1,n-1] = 1
    # Set up the lower half (spline conditions)
    for i in range(n-1):
        r1 = rbins[i]
        r2 = rbins[i+1]
        h1 = r2-r1
        if (i==0):
            # Set left boundary conditions
            if (bc == 'natural'):
                # Natural BC
                m[n,n] = 1.0
            elif ((bc == 'clamped') or (bc == 'mixed')):
                # Hermite (f' = 0)/clamped BC
                m[n,0] = 1.0/h1
                m[n,1] = -1.0/h1
                m[n,n] = h1/3.0
                m[n,n+1] = h1/6.0
            else:
                raise ValueError("Boundary condition '"+bc+"' not supported. Choose one of 'natural', 'clamped', or 'mixed'.")
        else:
            r0 = rbins[i-1]
            h0 = r1 - r0
            m[n+i,i-1] = -1.0/h0
            m[n+i,i] = 1.0/h0 + 1.0/h1
            m[n+i,i+1] = -1.0/h1
            m[n+i,n+i-1] = h0/6.0
            m[n+i,n+i] = (h0+h1)/3.0
            m[n+i,n+i+1] = h1/6.0
    # Set right boundary conditions
    if ((bc == 'natural') or (bc == 'mixed')):
        # Natural BC
        m[2*n-1,2*n-1] = 1.0
    elif (bc == 'clamped'):
        # Hermite (f' = 0)/clamped BC
        h1 = rbins[n-1] - rbins[n-2]
        m[2*n-1,n-2] = -1.0/h1
        m[2*n-1,n-1] = 1.0/h1
        m[2*n-1,2*n-2] = h1/6.0
        m[2*n-1,2*n-1] = h1/3.0
    else:
        raise ValueError("Boundary condition '"+bc+"' not supported. Choose one of 'natural', 'clamped', or 'mixed'.")

    return m

# Reconstruction algorithm with splines; numerical matrix inversion
def reconstruction_algorithm_spline(events, rbins, exp_settings=iaxo, bc='clamped'):
    rescale_factor = distance_factor*flux_to_events(1.0, gagg=1, exp_settings=exp_settings)
    n = len(rbins)
    m = create_matrix_spline(rbins, bc)

    b = np.zeros(2*n)
    b[:(n-1)] = events
    db2 = np.zeros((2*n,2*n))
    db2[:(n-1),:(n-1)] = np.diag(events)
    idb2 = np.zeros((2*n,2*n))
    idb2[:(n-1),:(n-1)] = np.diag(1/events)

    rsc2 = rescale_factor*rescale_factor
    b /= rescale_factor
    db2 /= rsc2
    idb2 *= rsc2

    x = np.linalg.solve(m, b)
    assert np.allclose(m@x, b)

    ys = x[:n]
    ms = x[n:]
    cov_ys_inv = (m.T)@idb2@m
    m_inv = inv(m)
    cov_ys = m_inv@db2@(m_inv.T)
    
    return rbins, ys, ms, cov_ys, cov_ys_inv

# Wrapper routine for reconstructing the Gamma bars using the algo above
def reconstruct_Gamma_bars_spline(grids, spot, rbins, bc='clamped'):
    gbs_spl, gbs, gbs_err, gbs_cov_inv = [], [], [], []
    n = len(rbins)
    for g in grids:
        counts = get_counts_in_rings_simple(g, spot, rbins)
        _, y, m, cov_ys, cov_ys_inv = reconstruction_algorithm_spline(counts, rbins, bc=bc)
        dy = np.sqrt(np.diag(cov_ys))[:n]
        h = [rbins[i+1]-rbins[i] for i in range(n-1)]
        c0 = [y[i] for i in range(n-1)]
        c1 = [(y[i+1]-y[i])/h[i] - h[i]*(m[i+1]+2*m[i])/6.0 for i in range(n-1)]
        c2 = [m[i]/2.0 for i in range(n-1)]
        c3 = [(m[i+1]-m[i])/(6.0*h[i]) for i in range(n-1)]
        c = np.array([c3, c2, c1, c0])
        p = PPoly(c, rbins)
        gbs_spl.append(p)
        gbs.append(y)
        gbs_err.append(dy)
        gbs_cov_inv.append(cov_ys_inv)
    return gbs_spl, np.array(gbs), np.array(gbs_err), gbs_cov_inv

def create_matrix_poly(rbins):
    n_rbins = len(rbins) - 1
    # Need a (n-1) x 4(n-1) matrix for y and 3 more polynomial coefficients
    # I_j = sum_i y_ij + c1,ij (r - r_i) + c2,ij (r - r_i)^2 + c3,ij (r - r_i)^3 ~ n_j
    # -> I_j =  M x_j ~ n_j with x_j = (y_j, c1,j, c2,j, c3,j)
    m = np.zeros((n_rbins,4*n_rbins))

    for i in range(n_rbins):
        # Compute recurring terms for efficiency
        r1 = rbins[i]
        r1s = r1*r1
        r1c = r1s*r1
        r1q = r1s*r1s
        r2 = rbins[i+1]
        r2s = r2*r2
        r2c = r2s*r2
        r2q = r2s*r2s
        if (i==0):
            m[i,i] = r2c/3.0
            m[i,i+n_rbins] = r2q/4.0
            m[i,i+2*n_rbins] = r2*r2q/5.0
            m[i,i+3*n_rbins] = r2s*r2q/6.0
        else:
            rt21 = np.sqrt(r2s - r1s)
            th21 = dthreehalf(r2,r1)
            ltt = np.log(r1/(r2 + rt21))
            m[i,i] = th21/3.0
            m[i,i+n_rbins] = (rt21*(8*r1c - 3*r1s*r2 - 8*r1*r2s + 6*r2c) + 3*r1q*ltt)/24.0
            m[i,i+2*n_rbins] = (rt21*(-28*r1q + 15*r1c*r2 + 16*r1s*r2s - 30*r1*r2c + 12*r2q) - 15*r1*r1q*ltt)/60.0
            m[i,i+3*n_rbins] = (rt21*(176*r1*r1q - 105*r1q*r2 - 32*r1c*r2s + 170*r1s*r2c - 144*r1*r2q + 40*r2*r2q) + 105*r1s*r1q*ltt)/240.0

        for j in range(i+1, n_rbins):
            r3 = rbins[j]
            r3s = r3*r3
            r3c = r3s*r3
            r3q = r3s*r3s
            r4 = rbins[j+1]
            r4s = r4*r4
            r4c = r4s*r4
            r4q = r4s*r4s
            rt41 = np.sqrt(r4s - r1s)
            rt42 = np.sqrt(r4s - r2s)
            rt31 = np.sqrt(r3s - r1s)
            rt32 = np.sqrt(r3s - r2s)
            lt1 = logterm(r4, r3, r1)
            lt2 = logterm(r3, r4, r2)

            m[i,j] = (r1s*(rt31-rt41) + r2s*(rt42-rt32) + r3s*(rt32-rt31) + r4s*(rt41-rt42))/3.0

            tmp1 =  rt31*r3*(-5*r1s + 2*r3s)
            tmp2 = -rt32*r3*(-5*r2s + 2*r3s)
            tmp = -8*r3*r4s + 6*r4c
            tmp3 =  rt41*(8*r1s*r3 - 3*r1s*r4 + tmp)
            tmp4 = -rt42*(8*r2s*r3 - 3*r2s*r4 + tmp)
            ltt = r1q*lt1 + r2q*lt2
            m[i,j+n_rbins] = (tmp1 + tmp2 + tmp3 + tmp4)/24.0 + ltt/8.0

            tmp1 =  rt31*(8*r1q + 9*r1s*r3s - 2*r3q)
            tmp2 = -rt32*(8*r2q + 9*r2s*r3s - 2*r3q)
            tmp = 20*r3s*r4s - 30*r3*r4c + 12*r4q
            tmp3 =  rt41*(-8*r1q - 20*r1s*r3s + 15*r1s*r3*r4 - 4*r1s*r4s + tmp)
            tmp4 = -rt42*(-8*r2q - 20*r2s*r3s + 15*r2s*r3*r4 - 4*r2s*r4s + tmp)
            ltt = -r3*ltt
            m[i,j+2*n_rbins] = (tmp1 + tmp2 + tmp3 + tmp4)/60.0 + ltt/4.0

            tmp1 =  rt31*(-81*r1q*r3 - 28*r1s*r3c + 4*r3*r3q)
            tmp2 = -rt32*(-81*r2q*r3 - 28*r2s*r3c + 4*r3*r3q)
            tmp = -80*r3c*r4s + 180*r3s*r4c - 144*r3*r4q + 40*r4*r4q
            tmp3 =  rt41*(96*r1q*r3 + 80*r1s*r3c - 15*r1q*r4 - 90*r1s*r3s*r4 + 48*r1s*r3*r4s + tmp - 10*r1s*r4c)
            tmp4 = -rt42*(96*r2q*r3 + 80*r2s*r3c - 15*r2q*r4 - 90*r2s*r3s*r4 + 48*r2s*r3*r4s + tmp - 10*r2s*r4c)
            ltt = r1q*lt1*(r1s + 6*r3s) + r2q*lt2*(r2s + 6*r3s)
            m[i,j+3*n_rbins] = (tmp1 + tmp2 + tmp3 + tmp4)/240.0 + ltt/16.0
    return m


#############################
#  Theoretical Gamma bars   #
#############################

def numerical_Gamma_bars(rbins, ebins, gagg, solar_model):
    res = []
    n_ebins = len(ebins)-1
    g2 = gagg*gagg
    integrand = lambda om, r: om*om * solar_model.primakoff_rate(om, r) / (2.0*np.pi*np.pi)
    for r in rbins:
        tmp = [r]
        tmp += [g2*quad(integrand, ebins[i], ebins[i+1], args=(r), epsabs=0)[0] for i in range(n_ebins)]
        res.append(tmp)
    return np.array(res)

def parametric_Gamma_bar(rbins, ebins, gagg, ks, temp):
    grid = []
    n_ebins, n_rbins = len(ebins)-1, len(rbins)-1
    g2 = gagg*gagg
    for j in range(n_rbins):
        temp_grid = [rbins[j]]
        integrand = lambda om: om*om * parametric_primakoff_rate(om, 1, ks[j], temp[j]) / (2.0*np.pi*np.pi)
        temp_grid += [g2*quad(integrand, a=ebins[i], b=ebins[i+1], epsabs=0, limit=100)[0] for i in range(n_ebins)]
        grid.append(temp_grid)
    return np.array(grid)


######################################
#  Calculate the theoretical grids   #
######################################

def calculate_theoretical_grids(x, ebins, rbins, spot, obs):
    n_ebins, n_rbins = len(ebins)-1, len(rbins)-1
    kwargs = { 'points': rbins[1:-1] }
    alt_integrator = lambda func, a, b: default_integrator(func, a, b, **kwargs)
    ny, nx = obs[0].shape

    # Integrate out energies
    def fpp(om, j):
        return distance_factor*om*om*parametric_primakoff_rate(om, x[0], x[j+1], x[n_rbins+j+1])/(2.0*np.pi*np.pi)
    
    primakoff_terms = [[quad(fpp, ebins[i], ebins[i+1], args=(j), epsabs=0, epsrel=1e-7, limit=200)[0] for j in range(n_rbins)] for i in range(n_ebins)]

    def fp(rho, i):
        red_rho = np.array([rho]+list(np.array(rbins)[np.array(rbins)>rho]))
        i0 = len(rbins) - len(red_rho)
        rho2 = rho*rho
        sqrt_terms = np.sqrt(red_rho*red_rho - rho2)
        fps = [primakoff_terms[i][j] * (sqrt_terms[j+1-i0] - sqrt_terms[j-i0]) for j in range(i0,n_rbins)]
        return np.sum(fps)  # no rho factor needed here

    theo = []
    for i in range(n_ebins):
        radial_flux = lambda rho, rs: fp(rho/rs, i) / (2.0*np.pi*rs*rs) # still need to divide by 2pi here
        grid_integral = calculate_counts_on_grid(nx, ny, spot, radial_flux, integrator=alt_integrator)
        theo.append(flux_to_events(grid_integral, gagg=x[0]))
    return theo

def calculate_true_theoretical_grids(solar_model, ebins, spot, grid_dims, gagg=0.5):
    r_max = 0.999*0.97
    n_ebins = len(ebins)-1
    nx, ny = grid_dims

    # Integrate out inner radii and energies to approximate differential fluxes as a function of rho
    def fp1(om, rho):
        inner_integrand = lambda r: r * solar_model.primakoff_rate(om, r) / np.sqrt(r*r - rho*rho)
        inner_integral = quad(inner_integrand, rho, r_max)[0]
        return om*om * inner_integral / (2.0*np.pi*np.pi) # no rho factor needed here

    def fp2(rho, i):
        return distance_factor * quad(fp1, ebins[i], ebins[i+1], args=(rho), limit=200)[0]

    rhos = np.linspace(0, r_max, 101)
    integrands = [[fp2(rho, i) for rho in rhos[:-1]]+[0] for i in range(n_ebins)]
    interp_funs = [CubicSpline(rhos, integrands[i]) for i in range(n_ebins)]

    theo = []
    for i in range(n_ebins):
        radial_flux = lambda rho, rs: interp_funs[i](rho/rs) / (2.0*np.pi*rs*rs) # still need to divide by 2pi here
        grid_integral = calculate_counts_on_grid(nx, ny, spot, radial_flux, r_max=r_max*spot[2])
        theo.append(gagg*gagg*flux_to_events(grid_integral, gagg=gagg))
    return theo