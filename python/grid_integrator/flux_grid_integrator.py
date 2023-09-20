import numpy as np
from typing import Callable, Tuple

from .grid_integrator import *

def calculate_counts_on_grid(nx: int, ny: int, spot: Tuple[float, float, float], radial_flux: Callable, r_max: float = None, integrator: Callable = default_integrator) -> np.ndarray:
    """
    Compute the counts in each grid bin given a radially symmetric function f(r).

    Args:
        nx (int): Number of bins along x-axis.
        ny (int): Number of bins along y-axis.
        spot (tuple): Tuple of integers representing the center of the grid.
        radial_flux (callable): A function that takes in a float r and returns the value of f(r).
        r_max (float): Cutoff value for the radial itegration. Default is spot[2]
        integrator (callable): A function used to compute the integral. Default is default_integrator (based on scipy.integrate.quad).

    Returns:
        np.ndarray: Array of shape (nx, ny) containing the counts in each bin of the grid.
    """
    r_scale = spot[2]
    r_max = r_scale if r_max==None else r_max
    arr = np.zeros((nx, ny))
    # Define a new radial_flux_rescaled function using lambda to apply r_scale
    radial_flux_rescaled = lambda r: radial_flux(r, r_scale)
    # Loop through the array and compute the counts in each bin
    for i in range(nx):
        for j in range(ny):
            # Calculate the boundaries of the square in coordinates where c = (0,0)
            x1 = i - spot[0]
            x2 = x1 + 1
            y1 = j - spot[1]
            y2 = y1 + 1
            # Calculate the overlap and multiply by the number of counts in the square
            arr[i, j] = rectangular_integral(
                radial_flux_rescaled, x1, x2, y1, y2, integrator, r_cutoff=r_max
            )
    return arr


def calculate_grid_overlaps(nx: int, ny: int, spot: Tuple[float, float, float], r_max: float = None) -> np.ndarray:
    """
    Routine to compute the overlap for each grid pixel. This uses the grid integration routines
    for f(r) = 1.

    Args:
        nx (int): Number of bins along x-axis.
        ny (int): Number of bins along y-axis.
        spot (tuple): Tuple of integers representing the center of the grid.
        r_max (float): Cutoff value for the radial itegration. Default is spot[2].

    Returns:
        np.ndarray: Array of shape (nx, ny) containing the counts in each bin of the grid.
    """
    r_scale = spot[2]
    r_max = r_scale if r_max==None else r_max
    arr = np.zeros((nx, ny))
    for i in range(nx):
        for j in range(ny):
            # Calculate the boundaries of the square in coordinates where c = (0,0)
            x1 = i - spot[0]
            x2 = x1 + 1
            y1 = j - spot[1]
            y2 = y1 + 1
            # Calculate the overlap and multiply by the number of counts in the square
            arr[i, j] = rectangular_integral(
                1, x1, x2, y1, y2, r_cutoff=r_scale*r_max
            )
    return arr


def get_counts_in_rings(grid: np.ndarray, spot: Tuple[float, float, float], n_rbins: int = 15) -> tuple:
    """
    Determine the radial bins from all counts such that they contain approximately the same number of counts.
    Return the radii and the counts in each bin.

    Args:
        grid (numpy.ndarray): A 2D array of counts.
        spot (tuple): A tuple of the form (x, y, r) representing the coordinates of the circle centre and its radius.
        n_rbins (int): The number of radial bins.

    Returns:
        tuple: A tuple containing the following arrays:
            * rebinned_radii (numpy.ndarray): An array of the bin edges for the rebinned data.
            * rebinned_counts (numpy.ndarray): An array of the counts in each bin of the rebinned data.
            * radii (numpy.ndarray): An array of the bin edges for the original data.
            * counts (numpy.ndarray): An array of the counts in each bin of the original data.
    """

    # Calculate the radii for the original data
    r_scale = spot[2]
    radii = np.linspace(0, r_scale, 100 * n_rbins)
    integrals = [0]

    # Calculate the overlap of the circles with the grid squares
    for r in radii[1:]:
        s = 0
        for i, ai in enumerate(grid):
            for j, counts in enumerate(ai):
                if counts > 0:
                    # Calculate the boundaries of the square in coordinates where c = (0,0)
                    x1 = i - spot[0]
                    x2 = x1 + 1
                    y1 = j - spot[1]
                    y2 = y1 + 1
                    s += counts * rectangular_integral(1, x1, x2, y1, y2, r_cutoff=r)
        integrals.append(s)
    events = np.diff(integrals)

    # Re-bin data so that bins contain about an equal number of events
    events_total = events.sum()
    target_events_per_bin = events_total / n_rbins
    rebinned_events, rebinned_radii = np.zeros(n_rbins), np.zeros(n_rbins + 1)
     # Add entries from counts_in_rings to rebinned_counts until the target count is reached
    j, j_max = 0, len(events)
    for i in range(n_rbins):
        while (rebinned_events[i] < target_events_per_bin) and (j < j_max):
            rebinned_events[i] += events[j]
            rebinned_radii[i + 1] = radii[j+1]
            j += 1
        if j >= j_max:
            break
    rebinned_radii, rebinned_events = (
        rebinned_radii[:(i+2)],
        rebinned_events[:(i+1)],
    )
    assert np.abs(np.sum(rebinned_events)/events_total - 1.0) < 1e-10
    assert radii[-1] == rebinned_radii[-1]

    return rebinned_radii/r_scale, rebinned_events, radii/r_scale, events


# Compute the number of counts from a grid given pre-defined radial bins
def get_counts_in_rings_simple(grid: np.ndarray, spot: Tuple[float, float, float], rbins: np.ndarray) -> np.ndarray:
    """
    Compute the number of counts from a grid given pre-defined radial bins.

    Args:
        grid (np.ndarray): An array of shape containing the counts.
        spot (Tuple[float, float, float]): A tuple containing the center (cx, cy) and radius (r) of the circular spot.
        rbins (np.ndarray): An array of shape (K,) containing the radial bin edges.

    Returns:
        np.ndarray: An array of shape (K-1,) containing the counts in each radial bin.
    """
    
    cx, cy, r0 = spot
    # Calculate the overlap of the circles with the grid squares
    integrals = [0]
    for r in rbins[1:]:
        s = 0
        for i, gi in enumerate(grid):
            for j, counts in enumerate(gi):
                if counts > 0:
                    # Calculate the boundaries of the square in coordinates where c = (0,0)
                    x1 = i - cx
                    x2 = x1 + 1
                    y1 = j - cy
                    y2 = y1 + 1
                    s += counts * rectangular_integral(1, x1, x2, y1, y2, r_cutoff=r*r0)
        integrals.append(s)
    counts = np.diff(integrals)

    return counts