import pytest
import numpy as np
from math import isclose

from pymoodeng.orbital_conversions import (
    Vector3D, Position, Velocity, StateVector,
    orbit_radius, elements_4_apsides, create_orbit,
    angular_momentum, node_vector, eccentricity_vector,
    specific_energy, elements_from_state,
    radius_from_alt, alt_from_radius, impulse_from_finite,
    OrbitalElements
)

@pytest.fixture
def test_vectors():
    return Position(7000, 0, 0), Velocity(0, 7.5, 0)

def test_vector3d_operations():
    v1 = Vector3D(1, 0, 0)
    v2 = Vector3D(0, 1, 0)
    cross = v1.cross(v2)
    assert cross == Vector3D(0, 0, 1)

def test_orbit_radius():
    r = orbit_radius(10000, 0.1, np.pi / 3)
    expected = 10000 * (1 - 0.1**2) / (1 + 0.1 * np.cos(np.pi / 3))
    assert isclose(r, expected, rel_tol=1e-6)

def test_elements_4_apsides():
    a, e = elements_4_apsides(8000, 6000)
    assert isclose(a, 7000)
    assert isclose(e, 0.142857, rel_tol=1e-6)

def test_create_orbit_dimensions():
    U, V, W = create_orbit(np.radians(45), np.radians(30), np.radians(60), np.radians(10))
    assert U.shape == (3,)
    assert V.shape == (3,)
    assert W.shape == (3,)

def test_angular_momentum(test_vectors):
    p, v = test_vectors
    h = angular_momentum(p, v)
    assert isinstance(h, Vector3D)
    assert isclose(h.z, 7000 * 7.5)

def test_node_vector():
    h = Vector3D(0, 0, 1)
    n = node_vector(h)
    assert isinstance(n, Vector3D)
    assert n == Vector3D(1, 0, 0)

def test_eccentricity_vector(test_vectors):
    p, v = test_vectors
    mu = 398600.4418  # Earth's mu
    e_vec = eccentricity_vector(p, v, mu)
    assert isinstance(e_vec, Vector3D)
    assert not np.isnan(e_vec.array).any()

def test_specific_energy(test_vectors):
    p, v = test_vectors
    mu = 398600.4418
    energy = specific_energy(p, v, mu)
    r = np.linalg.norm(p.array)
    v_mag = np.linalg.norm(v.array)
    expected = 0.5 * v_mag**2 - mu / r
    assert isclose(energy, expected, rel_tol=1e-6)

def test_elements_from_state(test_vectors):
    p, v = test_vectors
    mu = 398600.4418
    elements = elements_from_state(mu, p, v)
    assert isinstance(elements, OrbitalElements)
    assert 0 <= elements.e < 1

def test_radius_altitude_conversion():
    R = 6371
    alt = 500
    r = radius_from_alt(alt, R)
    assert r == R + alt
    assert alt_from_radius(r, R) == alt

def test_impulse_from_finite():
    a = 0.01  # m/s^2
    t = 1000  # seconds
    delta_v = impulse_from_finite(a, t)
    assert isclose(delta_v, 10.0)