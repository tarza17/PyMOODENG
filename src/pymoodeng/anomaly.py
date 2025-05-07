# Using code from http://www.jgiesen.de/kepler/kepler.html
# Implemnetation based on the following paper: [A Practical Method for Solving the Kepler Equation][1] by Marc A. Murison from the U.S. Naval Observatory
# [1] https://www.researchgate.net/publication/271214938_A_Practical_Method_for_Solving_the_Kepler_Equation

#Compatibility reasons
from __future__ import absolute_import, division, print_function

from collections import namedtuple
from contextlib import contextmanager
from copy import deepcopy
from typing import Literal
from dataclasses import dataclass
import numpy as np

# from numpy import atan2, floor, fmod, ising, isnan
# from numpy import arrcos as acos
# from numpy import cos, dot, sin, sqrt
# from numpy.linalg import norm
# #from math import atan2, floor, fmod, ising, isnan #chatqpt advice
#from scipy.constants import pi

MAX_ITER = 100

# Custom exception Exception

class ConvergenceError(Exception):
    pass


#Anomaly class
'''
class Anomaly(object):
    '''
   # Class representing an anomaly (mean, eccentric, or true)
'''
    def __init__(self, **kwargs):
        super()().__init__()
        
        valid_args = set(['M' , 'E' , 'theta'])

        extra_args = set(kwargs.keys()) - valid_args
        if extra_args:
            raise TypeError('Invalid arguments: {}'.format(extra_args))
        if not kwargs:
            raise TypeError('No arguments provided.')
        
    def M(self, e):
        if self.key == 'M':
            return self.anomaly
        if self.key == 'E':
            return mean_anomaly_from_eccentric(e, self.anomaly)
        if self.key == 'theta':
            return mean_anomaly_from_true(e, self.anomaly)
        
    def E(self, e):
        if self.key == 'M':
            return eccentric_anomaly_from_mean(e, self.anomaly)
        elif self.key == 'E':
            return self.anomaly
        elif self.key == 'theta':
            return eccentric_anomaly_from_true(e, self.anomaly)
        
    def theta(self, e):
        if self.key == 'M':
            return true_anomaly_from_mean(e, self.anomaly)
        elif self.key == 'E':
            return true_anomaly_from_eccentric(e, self.anomaly)
        elif self.key == 'theta':
            return self.anomaly
'''

@dataclass
class Anomaly:
    """
    Represents an orbital anomaly and provides conversions between 
    mean (M), eccentric (E), and true (theta) anomalies.

    Attributes:
        value (float): The value of the anomaly.
        type (Literal['M', 'E', 'theta']): The type of the anomaly:
            - 'M': Mean anomaly
            - 'E': Eccentric anomaly
            - 'theta': True anomaly

    Methods:
        M(e: float) -> float:
            Returns the mean anomaly, computed from the stored value if necessary.

        E(e: float) -> float:
            Returns the eccentric anomaly, computed from the stored value if necessary.

        theta(e: float) -> float:
            Returns the true anomaly, computed from the stored value if necessary.
    """

    value: float
    type: Literal['M', 'E', 'theta']

    def __post_init__(self):
        valid_types = {'M', 'E', 'theta'}
        if self.type not in valid_types:
            raise ValueError(f"Invalid anomaly type: '{self.type}'. Must be one of {valid_types}.")

    def M(self, e: float) -> float:
        if self.type == 'M':
            return self.value
        elif self.type == 'E':
            return mean_anomaly_from_eccentric(e, self.value)
        elif self.type == 'theta':
            return mean_anomaly_from_true(e, self.value)

    def E(self, e: float) -> float:
        if self.type == 'M':
            return eccentric_anomaly_from_mean(e, self.value)
        elif self.type == 'E':
            return self.value
        elif self.type == 'theta':
            return eccentric_anomaly_from_true(e, self.value)

    def theta(self, e: float) -> float:
        if self.type == 'M':
            return true_anomaly_from_mean(e, self.value)
        elif self.type == 'E':
            return true_anomaly_from_eccentric(e, self.value)
        elif self.type == 'theta':
            return self.value


#Anomaly calcs

'''
def eccentric_anomaly_from_mean(e, M, tolerance=1e-10):

    Calc eccenetric anomaly from mean anomaly

    Mnorm = np.fmod(M, 2 * np.pi)
    E0 = M + (-1 / 2 * e**3 + e + (e**2 + 3 / 2 * np.cos(M) * e ** 3) * np.cos(M)) * np.sin(M)
    dE = tolerance + 1
    count = 0
    while dE > tolerance:
        t1 = np.cos(E0)
        t2 = -1 + e * t1
        t3 = np.sin(E0)
        t4 = e * t3
        t5 = -E0 + t4 + Mnorm
        t6 = t5 / (1 / 2 * t5 * t4 / t2 + t2)
        E = E0 - t5 / ((1 / 2 * t3 - 1 / 6 * t1 * t6) * e * t6 + t2)
        dE = abs(E - E0)
        E0 = E
        count += 1
        if count == MAX_ITER:
            raise ConvergenceError('Did not converge after {n} iterations. (e={e!r}, M={M!r})'.format(n=MAX_ITER, e=e, M=M))
        return E
'''

def eccentric_anomaly_from_mean(e, M, tolerance=1e-8):  # Using relaxed tolerance too
    """
    Compute the eccentric anomaly from the mean anomaly using the Newton-Raphson method.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        M (float): Mean anomaly in radians.
        tolerance (float, optional): Convergence tolerance. Default is 1e-8.

    Returns:
        float: Eccentric anomaly in radians.

    Raises:
        ConvergenceError: If the Newton-Raphson iteration fails to converge within MAX_ITER iterations.
    """
    Mnorm = np.fmod(M, 2 * np.pi)

    # Initial guess (Mnorm is usually good for Newton-Raphson)
    E0 = Mnorm
    # Optional refinement for high eccentricity, but Mnorm often suffices
    # if e > 0.8:
    # E0 = np.pi

    dE = tolerance + 1
    count = 0
    E = E0  # Start with E = E0


    while dE > tolerance:
        # Newton-Raphson step: E_new = E - f(E) / f'(E)
        # f(E) = E - e * sin(E) - Mnorm
        # f'(E) = 1 - e * cos(E)
        E_new = E - (E - e * np.sin(E) - Mnorm) / (1.0 - e * np.cos(E))

        dE = abs(E_new - E)
        E = E_new  # Update E for the next iteration
        count += 1

        if count == MAX_ITER:
            raise ConvergenceError(
                f'Newton-Raphson did not converge after {MAX_ITER} iterations. '
                f'(e={e!r}, M={Mnorm!r}, E={E!r}, dE={dE!r})'  # More info
            )
    return E

def eccentric_anomaly_from_true(e, theta):
    """
    Compute the eccentric anomaly from the true anomaly.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        theta (float): True anomaly in radians.

    Returns:
        float: Eccentric anomaly in radians.
    """
    E = np.atan2(np.sqrt(1 - e**2) * np.sin(theta), np.cos(theta) + e)
    E = np.mod(E, 2 * np.pi)
    return E

def mean_anomaly_from_eccentric(e, E):
    """
    Compute the mean anomaly from the eccentric anomaly.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        E (float): Eccentric anomaly in radians.

    Returns:
        float: Mean anomaly in radians.
    """
    M = E - e * np.sin(E)
    M = np.mod(M, 2 * np.pi)
    return M


def mean_anomaly_from_true(e, theta):
    """
    Compute the mean anomaly from the true anomaly.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        theta (float): True anomaly in radians.

    Returns:
        float: Mean anomaly in radians.
    """
    E = eccentric_anomaly_from_true(e, theta)
    M = E - e * np.sin(E)
    M = np.mod(M, 2 * np.pi)
    return M

def true_anomaly_from_eccentric(e, E):
    """
    Compute the true anomaly from the eccentric anomaly.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        E (float): Eccentric anomaly in radians.

    Returns:
        float: True anomaly in radians.
    """
    #theta = np.atan2(np.sqrt(1 - e**2) * np.sin(E), np.cos(E) + e)
    theta = np.arctan2(np.sqrt(1 - e ** 2) * np.sin(E), np.cos(E) + e) #Z
    theta = np.mod(theta, 2 * np.pi)
    return theta

def true_anomaly_from_mean(e, M, tolerance=1e-10):
    """
    Compute the true anomaly from the mean anomaly.

    This function internally computes the eccentric anomaly via Newton-Raphson, 
    then converts it to the true anomaly.

    Args:
        e (float): Orbital eccentricity (0 <= e < 1).
        M (float): Mean anomaly in radians.
        tolerance (float, optional): Convergence tolerance for Newton-Raphson. Default is 1e-10.

    Returns:
        float: True anomaly in radians.

    Raises:
        ConvergenceError: If Newton-Raphson fails to converge within MAX_ITER iterations.
    """
    E = eccentric_anomaly_from_mean(e, M, tolerance)
    theta = true_anomaly_from_eccentric(e, E)
    return theta














