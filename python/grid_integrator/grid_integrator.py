import numpy as np
from typing import Callable
from scipy.integrate import quad


def default_integrator(func: Callable, a: float, b: float, **kwargs) -> float:
    """
    Default integrator uses the standard scipy quad routine.
    The user could modify this to e.g. include difficult integration points etc.

    Args:
        func (callable): A function to be integrated.
        a (float): Lower limit of integration.
        b (float): Upper limit of integration.
        **kwargs: Optional arguments for scipy.integrate.quad.

    Returns:
        float: The numerical value of the definite integral.
    """
    epsabs = kwargs.pop("epsabs", 1e-6)
    epsrel = kwargs.pop("epsrel", 1e-6)
    limit = kwargs.pop("limit", 200)
    res, _ = quad(func, a, b, epsabs=epsabs, epsrel=epsrel, limit=limit, **kwargs)
    return res


# Smallest computation unit that only considers the 1st (top-left) quadrant of the integration grid
def rectangular_integral_1stquadrant(
    func: Callable, x1: float, x2:float, y1:float, y2: float, integrator: Callable = default_integrator, r_cutoff:float = np.inf
) -> tuple:
    """
    Computes the integral of a radially symmetric function over the rectangle x1,x2,y1,y2 in the first quadrant.

    The integration region has sharp corners at r = x1*y2 and r = x2*y1, so we need to break up the integration 
    into 3 parts (or less).

    Args:
        func (callable): The function to be integrated.
        x1 (float): Lower limit of integration along the x-axis.
        x2 (float): Upper limit of integration along the x-axis.
        y1 (float): Lower limit of integration along the y-axis.
        y2 (float): Upper limit of integration along the y-axis.
        integrator (callable, optional): Integration routine to use. Defaults to default_integrator.
        r_cutoff (float, optional): Maximum radius to integrate over. Defaults to np.inf.

    Returns:
        tuple: A tuple containing the numerical value of the definite integral and the integration limits.
    """
    # Precompute some quantities to speed up computations
    x1s, x2s, y1s, y2s = x1 * x1, x2 * x2, y1 * y1, y2 * y2
    x1y1, x2y2 = x1 * y1, x2 * y2
    rx1y1s = x1s + y1s
    rx2y2s = x2s + y2s
    r_min = np.sqrt(rx1y1s)
    r_max = min(np.sqrt(rx2y2s), r_cutoff)
    rx1y2 = np.sqrt(x1s + y2s)
    rx2y1 = np.sqrt(x2s + y1s)

    def f1(r):  # r_min to min(rx1y2, rx2y1)
        r2 = r * r
        return (
            r
            * func(r)
            * np.arcsin((r2 - rx1y1s) / (np.sqrt((r2 - x1s) * (r2 - y1s)) + x1y1))
        )

    def f3(r):  # max(rx1y2, rx2y1) to r_max
        r2 = r * r
        return (
            r
            * func(r)
            * np.arcsin((rx2y2s - r2) / (np.sqrt((r2 - x2s) * (r2 - y2s)) + x2y2))
        )

    if r_cutoff <= r_min:
        f2 = None
        return 0, f1, f2, f3, r_min, r_max, rx1y2, rx2y1
    elif rx1y2 < rx2y1:
        diffy2y1 = y2s - y1s

        def f2(r):
            r2 = r * r
            return (
                r
                * func(r)
                * np.arcsin(
                    diffy2y1 / (y2 * np.sqrt(r2 - y1s) + y1 * np.sqrt(r2 - y2s))
                )
            )

        rx1y2 = min(rx1y2, r_cutoff)
        rx2y1 = min(rx2y1, r_cutoff)
        result = (
            (integrator(func=f1, a=r_min, b=rx1y2) if r_min < rx1y2 else 0)
            + (integrator(func=f2, a=rx1y2, b=rx2y1) if rx1y2 < rx2y1 else 0)
            + (integrator(func=f3, a=rx2y1, b=r_max) if rx2y1 < r_max else 0)
        )

    elif rx2y1 < rx1y2:
        diffx2x1 = x2s - x1s

        def f2(r):
            r2 = r * r
            return (
                r
                * func(r)
                * np.arcsin(
                    diffx2x1 / (x2 * np.sqrt(r2 - x1s) + x1 * np.sqrt(r2 - x2s))
                )
            )

        rx1y2 = min(rx1y2, r_cutoff)
        rx2y1 = min(rx2y1, r_cutoff)
        result = (
            (integrator(func=f1, a=r_min, b=rx2y1) if r_min < rx2y1 else 0)
            + (integrator(func=f2, a=rx2y1, b=rx1y2) if rx2y1 < rx1y2 else 0)
            + (integrator(func=f3, a=rx1y2, b=r_max) if rx1y2 < r_max else 0)
        )
    # rx1y2==rx2y1 and only need to do 2 integrals
    else:
        f2 = None
        rx1y2 = min(rx1y2, r_cutoff)
        result = (integrator(func=f1, a=r_min, b=rx1y2) if r_min < rx1y2 else 0) + (
            integrator(func=f3, a=rx1y2, b=r_max) if rx1y2 < r_max else 0
        )

    return result, f1, f2, f3, r_min, r_max, rx1y2, rx2y1


def calculate_circular_segment(ch: float, r: float) -> float:
    """
    Calculate the area of a circular segment of radius r and chord length ch.

    Args:
        ch (float): chord length of the circular segment
        r (float): radius of the circle

    Returns:
        float: area of the circular segment
    """
    theta = 2.0 * np.arcsin(0.5 * ch / r)
    return 0.5 * r * r * (theta - np.sin(theta))


def rectangular_integral_1stquadrant_constant_f(x1: float, x2: float, y1: float, y2: float, r_cutoff: float = np.inf) -> float:
    """
    Calculates the integral of a function f(x,y) over the rectangle [x1,x2] x [y1,y2] in the first quadrant,
    assuming that f(r) is constant, where r is the distance from the origin to (x,y).
    
    Args:
    x1 (float): Lower bound of the x-coordinate.
    x2 (float): Upper bound of the x-coordinate.
    y1 (float): Lower bound of the y-coordinate.
    y2 (float): Upper bound of the y-coordinate.
    r_cutoff (float, optional): Maximum radius of the circular segment. Default is infinity.
    
    Returns:
    float: Integral of f(x,y) over the rectangle [x1,x2] x [y1,y2] in the first quadrant.
    """
    r = r_cutoff
    r2 = r * r
    x1s, x2s, y1s, y2s = x1 * x1, x2 * x2, y1 * y1, y2 * y2
    d1, d2 = np.sqrt(x1s + y1s), np.sqrt(x2s + y2s)
    lx, ly = x2 - x1, y2 - y1
    lxs, lys = lx * lx, ly * ly
    a0 = lx * ly

    # If the circle is too small (large), return 0 (1)
    if r <= d1:
        return 0
    elif r >= d2:
        return a0

    dx1 = np.sqrt(r2 - y1s) - x1
    dy1 = np.sqrt(r2 - x1s) - y1

    xcheck, ycheck = dx1 < lx, dy1 < ly

    if xcheck and ycheck:
        a1 = 0.5 * dx1 * dy1  # = triangle in the left corner
        dx1s, dy1s = dx1 * dx1, dy1 * dy1
        ch = np.sqrt(dx1s + dy1s)
        a2 = calculate_circular_segment(ch, r)
        return a1 + a2
    elif xcheck and (not ycheck):
        dx2 = np.sqrt(r2 - y2s) - x1
        a1 = 0.5 * ly * (dx1 + dx2)  # = rectangle + triangle = ly*dx2 + ly*(ly - dx2)/2
        dx12 = dx1 - dx2
        ch = np.sqrt(dx12 * dx12 + lys)
        a2 = calculate_circular_segment(ch, r)
        return a1 + a2
    elif (not xcheck) and ycheck:
        dy2 = np.sqrt(r2 - x2s) - y1
        a1 = 0.5 * lx * (dy1 + dy2)  # = rectangle + triangle
        dy12 = dy1 - dy2
        ch = np.sqrt(dy12 * dy12 + lxs)
        a2 = calculate_circular_segment(ch, r)
        return a1 + a2
    else:
        dx2 = np.sqrt(r2 - y2s) - x1
        dy2 = np.sqrt(r2 - x2s) - y1
        dx3, dy3 = lx - dx2, ly - dy2
        a1 = a0 - 0.5 * dx3 * dy3  # = full rectangle - triangle in the right corner
        ch = np.sqrt(dx3 * dx3 + dy3 * dy3)
        a2 = calculate_circular_segment(ch, r)
        return a1 + a2


# Compute the integral
def rectangular_integral(
    func, x1, x2, y1, y2, integrator=default_integrator, r_cutoff=np.inf
):
    """
    Integrate func over 2d rectangular region:
         x1 < x < x2
         y1 < y < y2
    func is a radially symmetric function of distance only.
    Does it as 1d integral in r:
         integral of f( sqrt(x^2+y^2) )dxdy from x1 to x2 and y1 to y2
           = int of f(r)*r*deltheta(r)*dr from r=sqrt(x1^2 + y1^2) to sqrt(x2^2 + y2^2),
           where deltheta(r) is the length of the arc of radius r cutting through the rectangle

    integrator = e.g. the trapz or rombint function of a TrapIntegrator instance
         integrator needs to take the following inputs and return a number:
             func, a, b, inittrap, and everything in **integratorkwargs
    integratorkwargs = extra kwargs to pass to integratorfunc (e.g. epsrel, K)
    """

    # Flip rectangle into upper right quadrant;
    # if rectangle crosses an axis then break it in half
    if x1 < 0 and x2 <= 0:
        return rectangular_integral(func, abs(x2), -x1, y1, y2, integrator, r_cutoff)
    if y1 < 0 and y2 <= 0:
        return rectangular_integral(func, x1, x2, abs(y2), -y1, integrator, r_cutoff)
    if x1 < 0 and x2 > 0:
        return rectangular_integral(
            func, x1, 0, y1, y2, integrator, r_cutoff
        ) + rectangular_integral(func, 0, x2, y1, y2, integrator, r_cutoff)
    if y1 < 0 and y2 > 0:
        return rectangular_integral(
            func, x1, x2, y1, 0, integrator, r_cutoff
        ) + rectangular_integral(func, x1, x2, 0, y2, integrator, r_cutoff)

    # if we made it here then rectangle is in upper right quadrant
    if func == 1:
        res = rectangular_integral_1stquadrant_constant_f(x1, x2, y1, y2, r_cutoff)
    else:
        res = rectangular_integral_1stquadrant(
            func, x1, x2, y1, y2, integrator, r_cutoff
        )[0]
    return res