'''
This module defines classes and functions for simulating celestial bodies and their orbits.'''
#IMPORTS

from __future__ import absolute_import, division, print_function
from dataclasses import dataclass
from typing import Tuple
from contextlib import contextmanager
from copy import deepcopy

import numpy as np
from numpy import atan2, floor, fmod, isnan
from numpy import cos, dot, sin, sqrt
from numpy.linalg import norm
#from math import atan2, floor, fmod, ising, isnan #chatqpt advice
from scipy.constants import pi


# Generator and context manager

@contextmanager
def saved_state(orbit):
    """Context manager to restore the orbit state upon leaving the block.
    
    This context manager saves the state of the orbit object before executing
    the code inside the block. After the block completes, the orbit object
    is restored to its original state.
    
    Args:
        orbit: The orbit object whose state is to be saved and restored.
    
    Yields:
        None: This function is used within a `with` block and doesn't return a value.
    """
    state = deepcopy(orbit.__dict__)
    yield
    orbit.__dict__.update(state)


def lookahead(collection, fillvalue=None):
    """Generates a series with lookahead to the next item in the collection.
    
    This function iterates over a collection and yields each item alongside the
    next item in the collection. The last item is paired with the `fillvalue` (default `None`).
    
    Args:
        collection: The collection (e.g., list or iterable) to iterate over.
        fillvalue: The value to pair with the last item in the collection if no next item exists.
    
    Yields:
        Tuple: A pair of the current and next items in the collection.
    """
    first = True
    for next_item in collection:
        if first:
            first = False
        else:
            yield current_item, next_item
        current_item = next_item
    yield current_item, fillvalue


# 3d vector class for easy shenanigans #gpt advice

@dataclass
class Vector3D:
    x: float
    y: float
    z: float

    @property
    def array(self):
        """Returns the vector as a NumPy array.
        
        This property converts the 3D vector into a NumPy array for easier
        manipulation with vectorized operations.
        
        Returns:
            np.array: The vector as a NumPy array [x, y, z].
        """
        return np.array([self.x, self.y, self.z])

    def cross(self, other):
        """Calculates the cross product of this vector with another vector.
        
        Args:
            other (Vector3D): The vector to compute the cross product with.
        
        Returns:
            Vector3D: The resulting vector from the cross product.
        """
        result = np.cross(self.array, other.array)
        return Vector3D(*result)

    def __add__(self, other):
        """Adds this vector with another vector.
        
        Args:
            other (Vector3D): The vector to add to this vector.
        
        Returns:
            Vector3D: The resulting vector after addition.
        """
        result = self.array + other.array
        return Vector3D(*result)

    def __sub__(self, other):
        """Subtracts another vector from this vector.
        
        Args:
            other (Vector3D): The vector to subtract from this vector.
        
        Returns:
            Vector3D: The resulting vector after subtraction.
        """
        result = self.array - other.array
        return Vector3D(*result)

    def __repr__(self):
        """String representation of the Vector3D object."""
        return f"{self.__class__.__name__}(x={self.x}, y={self.y}, z={self.z})"

# Position and Velocity types for clarity
@dataclass
class Position(Vector3D):
    """Represents the position vector in 3D space."""
    pass

@dataclass
class Velocity(Vector3D):
    """Represents the velocity vector in 3D space."""
    pass

# Simple state container
@dataclass
class StateVector:
    """Represents the state vector of an object, including its position and velocity.
    
    Attributes:
        position (Position): The position vector of the object.
        velocity (Velocity): The velocity vector of the object.
    """
    position: Position
    velocity: Velocity

# Orbital elements container
@dataclass
class OrbitalElements:
    """Represents the orbital elements of an orbiting object.
    
    Attributes:
        a (float): Semi-major axis.
        e (float): Eccentricity.
        i (float): Inclination.
        raan (float): Right ascension of the ascending node.
        arg_pe (float): Argument of periapsis.
        theta (float): True anomaly.
    """
    a: float        # semi-major axis
    e: float        # eccentricity
    i: float        # inclination
    raan: float     # right ascension of ascending node
    arg_pe: float   # argument of periapsis
    theta: float    # true anomaly

# Angular momentum computation
def angular_momentum(position: Position, velocity: Velocity) -> Vector3D:
    """Calculates the angular momentum vector from position and velocity vectors.
    
    Args:
        position (Position): The position vector.
        velocity (Velocity): The velocity vector.
    
    Returns:
        Vector3D: The resulting angular momentum vector.
    """
    return position.cross(velocity)




#Orbital element 

def orbit_radius(a, e, theta):
    '''
    Calc orbit radius
    '''
    r = a * (1 - e^2) / (1 + e * cos(theta))
    return r


def elements_4_apsides(apo_radius, peri_radius):
    """
    Calc orbital elements from apoapsis and periapsis radii
    """
    a = (apo_radius + peri_radius) / 2
    e = (apo_radius - peri_radius) / (apo_radius + peri_radius)
    return a, e


def create_orbit(i, raan, arg_pe, theta) -> Tuple[Vector3D, Vector3D, Vector3D]:
    """
    Create an orbit from the given elements
    """
    u = arg_pe + theta

    sin_u = sin(u)
    cos_u = cos(u)
    sin_raan = sin(raan)
    cos_raan = cos(raan)
    sin_i = sin(i)
    cos_i = cos(i)

    U = np.array(
        [cos_u * cos_raan - sin_u * sin_raan * cos_i,
         cos_u * sin_raan + sin_u * cos_raan * cos_i,
         sin_u * sin_i]
    )

    V = np.array(
        [-sin_u * cos_raan - cos_u * sin_raan * cos_i,
         -sin_u * sin_raan + cos_u * cos_raan * cos_i,
         cos_u * sin_i]
    )

    W = np.array(
        [sin_raan * sin_i,
         -cos_raan * sin_i,
         cos_i]
    )

    return U, V, W

def node_vector(angular_momentum: Vector3D) -> Vector3D:
    """Calculates the node vector from the angular momentum vector.
    
    The node vector is perpendicular to both the orbital plane and the angular momentum vector.
    
    Args:
        angular_momentum (Vector3D): The angular momentum vector.
    
    Returns:
        Vector3D: The node vector.
    """
    tmp = Vector3D(0, 0, 1)
    return angular_momentum.cross(tmp)

def eccentricity_vector(position: Position, velocity: Velocity, mu: float) -> Vector3D:
    """Calculates the eccentricity vector from position and velocity vectors.
    
    The eccentricity vector points from the center of the orbit toward periapsis.
    
    Args:
        position (Position): The position vector.
        velocity (Velocity): The velocity vector.
        mu (float): The standard gravitational parameter.
    
    Returns:
        Vector3D: The resulting eccentricity vector.
    """
    r = position.array
    v = velocity.array
    h = angular_momentum(position, velocity).array
    e = (np.cross(v, h) / mu) - (r / norm(r))
    return Vector3D(*e)

def specific_energy(position: Position, velocity: Velocity, mu: float) -> float:
    """Calculates the specific orbital energy of the object.
    
    Specific energy is the sum of kinetic and potential energy per unit mass.
    
    Args:
        position (Position): The position vector.
        velocity (Velocity): The velocity vector.
        mu (float): The standard gravitational parameter.
    
    Returns:
        float: The specific orbital energy.
    """
    r = norm(position.array)
    v = norm(velocity.array)
    return (v^2 / 2) - (mu / r)

def elements_from_state(mu: float, position: Position, velocity: Velocity) -> OrbitalElements:
    #erőssen krdőjeles, mert tabbal készült
    """Calculates the orbital elements from state vectors (position and velocity).
    
    Args:
        mu (float): The standard gravitational parameter.
        position (Position): The position vector.
        velocity (Velocity): The velocity vector.
    
    Returns:
        OrbitalElements: The calculated orbital elements (semi-major axis, eccentricity, etc.).
    """
    r = position.array
    v = velocity.array
    h = angular_momentum(position, velocity).array
    n = node_vector(h).array

    # Eccentricity vector
    e_vec = eccentricity_vector(position, velocity, mu).array
    e = norm(e_vec)

    # Semi-major axis
    energy = specific_energy(position, velocity, mu)
    a = -mu / (2 * energy)

    # Inclination
    i = acos(h[2] / norm(h))

    if abs(i - 0) < 1e-9:
        # For non-inclined orbits, raan is undefined;
        # set to zero by convention
        raan = 0
    if abs(e - 0) < 1e-9:
        # For circular orbits, place periapsis
        # at ascending node by convention
        arg_pe = 0
    else:
        # Argument of periapsis is the angle between
        # eccentricity vector and its x component.
        arg_pe = acos(e_vec.x / norm(e_vec))
    '''
        # Right ascension of ascending node
        raan = atan2(n[1], n[0])

        # Argument of periapsis
        arg_pe = atan2(e_vec[2], np.dot(n, e_vec))

        # True anomaly
        theta = atan2(np.dot(r, v), np.dot(r, e_vec))


    '''
    if abs(e - 0) < 1e-9:
        if abs(i - 0) < 1e-9:
            # True anomaly is angle between position
            # vector and its x component.
            f = acos(r.x / norm(r))
            if v.x > 0:
                f = 2 * pi - f
        else:
            # True anomaly is angle between node
            # vector and position vector.
            f = acos(dot(n, r) / (norm(n) * norm(r)))
            if dot(n, v) > 0:
                f = 2 * pi - f
    else:
        if e_vec.z < 0:
            arg_pe = 2 * pi - arg_pe

        # True anomaly is angle between eccentricity
        # vector and position vector.
        theta = acos(dot(e_vec, r) / (norm(e_vec) * norm(r)))

        if dot(r, v) < 0:
            theta = 2 * pi - theta
    return OrbitalElements(a=a, e=e, i=i, raan=raan, arg_pe=arg_pe, theta=theta)




##Conversions

def radius_from_alt(altitude, radius_body):
    """Converts altitude to radius from the center of the body.
    
    Args:
        altitude (float): The altitude above the surface.
        radius_body (float): The radius of the celestial body.
    
    Returns:
        float: The radius from the center of the body.
    """
    return radius_body + altitude

def alt_from_radius(radius, radius_body):
    """Converts radius from the center of the body to altitude.
    
    Args:
        radius (float): The radius from the center of the celestial body.
        radius_body (float): The radius of the celestial body.
    
    Returns:
        float: The altitude above the surface.
    """
    return radius - radius_body


def impulse_from_finite(acceleration, duration):
    """Returns the impulsive velocity change for a finite thrust burn.
    
    Args:
        acceleration (float): The acceleration produced by the thrust.
        duration (float): The duration of the thrust application.
    
    Returns:
        float: The impulsive velocity change.
    """
    return acceleration * duration



