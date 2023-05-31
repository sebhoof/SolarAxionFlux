import pytest
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import dblquad

from ..tomography.tomography import create_matrix_poly

def test_matrix_integration():
   r0 = np.array([0.0, 0.3, 0.64826434, 1.0])
   y0 = 1e-35*np.array([0.71085, 0.3675, 0.15403692, 0.49735])
   
   y0_max = np.max(y0)
   interp = PchipInterpolator(r0, y0/y0_max)
   foo = lambda r, rho: r * rho * interp(r) / np.sqrt(r*r - rho*rho)
   integral1 = []
   for i in range(len(r0)-1):
      res = dblquad(foo, r0[i], r0[i+1], lambda rho: rho, lambda rho: 1, epsabs=1e-9, epsrel=1e-9)[0]
      integral1.append(res)
   integral1 = y0_max*np.array(integral1)

   matrix = create_matrix_poly(r0)
   coeff = np.flip(interp.c, axis=0).reshape(-1)
   integral2 = y0_max*(matrix@coeff)

   r = integral2 - integral1
   r = np.sqrt(np.sum(r*r))/np.sum(integral1)
   print("Rel. difference of integration via matrix muliplication vs numerical integration: {:.2e}".format(r))
   assert r < 1e-8