import pytest
import numpy as np

from ..tomography.fitting import log_prior
from ..constants import *

def test_prior():
   n_rbins = 6
   np1 = n_rbins+1
   rbins = np1*[1]
   n = 2*n_rbins+1
   xb = n*[[-np.inf, np.inf]]
   rng = np.random.default_rng(42)
   x = rng.random(n)
   x[1:np1] = xi_m*x[np1:]
   lp = log_prior(x, rbins, xb)
   assert -2.0*lp < 1e-30 