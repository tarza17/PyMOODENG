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
from numpy import atan2, floor, fmod, ising, isnan
from numpy import arrcos as acos
from numpy import cos, dot, sin, sqrt
from numpy.linalg import norm
#from math import atan2, floor, fmod, ising, isnan #chatqpt advice
from scipy.constants import pi

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

def eccentric_anomaly_from_mean(e, M, tolerance=1e-10):
    '''
    Calc eccenetric anomaly from mean anomaly
    '''
    Mnorm = fmod(M, 2 * pi)
    E0 = M + (-1 / 2 * e^3 + e + (e^2 + 3 / 2 * cos(M) * e ** 3) * cos(M)) * sin(M)
    dE = tolerance + 1
    count = 0
    while dE > tolerance:
        t1 = cos(E0)
        t2 = -1 + e * t1
        t3 = sin(E0)
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


def eccentric_anomaly_from_true(e, theta):
    '''
    Calc eccenetric anomaly from true anomaly
    '''
    E = atan2(sqrt(1 - e^2) * sin(theta), cos(theta) + e)
    E = mod(E, 2 * pi)
    return E

def mean_anomaly_from_eccentric(e, E):
    '''
    Calc mean anomaly from eccenetric anomaly
    '''
    M = E - e * sin(E)
    M = mod(M, 2 * pi)
    return M


def mean_anomaly_from_true(e, theta):
    '''
    Calc mean anomaly from true anomaly
    '''
    E = eccentric_anomaly_from_true(e, theta)
    M = E - e * sin(E)
    M = mod(M, 2 * pi)
    return M

def true_anomaly_from_eccentric(e, E):
    '''
    Calc true anomaly from eccenetric anomaly
    '''
    theta = atan2(sqrt(1 - e^2) * sin(E), cos(E) + e)
    theta = mod(theta, 2 * pi)
    return theta

def true_anomaly_from_mean(e, M, tolerance=1e-10):
    '''
    Calc true anomaly from mean anomaly
    '''
    E = eccentric_anomaly_from_mean(e, M, tolerance)
    theta = true_anomaly_from_eccentric(e, E)
    return theta














