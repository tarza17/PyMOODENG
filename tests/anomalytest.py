import pytest
import numpy as np
from pymoodeng.anomaly import (
    Anomaly, eccentric_anomaly_from_mean, eccentric_anomaly_from_true,
    mean_anomaly_from_eccentric, mean_anomaly_from_true,
    true_anomaly_from_eccentric, true_anomaly_from_mean,
    ConvergenceError
)

# Tolerance for floating point comparison
TOL = 1e-9

@pytest.mark.parametrize("e, M", [
    (0.0, 0.0),                # Circular orbit
    (0.5, np.pi / 4),          # Typical elliptical orbit
    (0.9, np.pi),              # High eccentricity
])
def test_eccentric_mean_roundtrip(e, M):
    E = eccentric_anomaly_from_mean(e, M)
    M_back = mean_anomaly_from_eccentric(e, E)
    assert np.isclose(M_back, M % (2*np.pi), atol=TOL)

@pytest.mark.parametrize("e, theta", [
    (0.1, 0.0),
    (0.3, np.pi/2),
    (0.7, np.pi),
])
def test_true_eccentric_roundtrip(e, theta):
    E = eccentric_anomaly_from_true(e, theta)
    theta_back = true_anomaly_from_eccentric(e, E)
    assert np.isclose(theta_back, theta % (2*np.pi), atol=TOL)

@pytest.mark.parametrize("e, theta", [
    (0.1, 0.0),
    (0.5, np.pi/3),
    (0.8, 2*np.pi/3),
])
def test_true_mean_roundtrip(e, theta):
    M = mean_anomaly_from_true(e, theta)
    theta_back = true_anomaly_from_mean(e, M)
    assert np.isclose(theta_back, theta % (2*np.pi), atol=TOL)

@pytest.mark.parametrize("value, kind, e", [
    (1.0, 'M', 0.5),
    (1.0, 'E', 0.5),
    (1.0, 'theta', 0.5),
])
def test_anomaly_conversions(value, kind, e):
    a = Anomaly(value=value, type=kind)
    M = a.M(e)
    E = a.E(e)
    theta = a.theta(e)

    # Roundtrips
    assert np.isclose(mean_anomaly_from_eccentric(e, E), M % (2*np.pi), atol=TOL)
    assert np.isclose(true_anomaly_from_mean(e, M), theta % (2*np.pi), atol=TOL)

def test_invalid_anomaly_type():
    with pytest.raises(ValueError):
        Anomaly(value=1.0, type='invalid')

def test_nonconvergence():
    # Force failure with bad initial guess by using extremely high eccentricity
    with pytest.raises(ConvergenceError):
        eccentric_anomaly_from_mean(0.9999999999, 1.0, tolerance=1e-16)
