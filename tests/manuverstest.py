import pytest
from scipy.constants import pi
from pymoodeng.manuvers import (  
    SetApocenterRadiusTo,
    SetApocenterAltitudeTo,
    ChangeApocenterBy,
    SetPericenterRadiusTo,
    SetPericenterAltitudeTo,
    ChangePericenterBy,
    SetInclinationTo,
    ChangeInclanationBy,
    PropagateAnomalyTo,
    PropagateAnomalyBy,
    Circularise,
    SetPerincenterHere,
    Maneuver,
)

# Mock orbit class for testing purposes
class MockOrbit:
    def __init__(self, apocenter_radius=10000, pericenter_radius=7000, inclination=0, arg_pe=0, v=7.5, M=0, f=0, t=0, body=None):
        self.apocenter_radius = apocenter_radius
        self.pericenter_radius = pericenter_radius
        self.inclination = inclination
        self.arg_pe = arg_pe
        self.v = v
        self.M = M
        self.f = f
        self.t = t
        self.body = body if body else MockBody()

    def propagate_anomaly_to(self, **kwargs):
        # Mock implementation
        pass

class MockBody:
    def __init__(self):
        self.mean_radius = 6371  # Example value for Earth's mean radius

# Test SetApocenterRadiusTo
def test_set_apocenter_radius_to():
    orbit = MockOrbit()
    operation = SetApocenterRadiusTo(apocenter_radius=12000)
    operation.__apply__(orbit)
    assert orbit.apocenter_radius == 12000

# Test SetApocenterAltitudeTo
def test_set_apocenter_altitude_to():
    orbit = MockOrbit()
    operation = SetApocenterAltitudeTo(apocenter_altitude=5629)  # Altitude above Earth's mean radius
    operation.__apply__(orbit)
    assert orbit.apocenter_radius == 5629 + orbit.body.mean_radius

# Test ChangeApocenterBy
def test_change_apocenter_by():
    orbit = MockOrbit(apocenter_radius=10000)
    operation = ChangeApocenterBy(delta=2000)
    operation.__apply__(orbit)
    assert orbit.apocenter_radius == 12000

# Test SetPericenterRadiusTo
def test_set_pericenter_radius_to():
    orbit = MockOrbit()
    operation = SetPericenterRadiusTo(pericenter_radius=8000)
    operation.__apply__(orbit)
    assert orbit.pericenter_radius == 8000

# Test SetPericenterAltitudeTo
def test_set_pericenter_altitude_to():
    orbit = MockOrbit()
    operation = SetPericenterAltitudeTo(pericenter_altitude=1629)  # Altitude above Earth's mean radius
    operation.__apply__(orbit)
    assert orbit.pericenter_radius == 1629 + orbit.body.mean_radius

# Test ChangePericenterBy
def test_change_pericenter_by():
    orbit = MockOrbit(pericenter_radius=7000)
    operation = ChangePericenterBy(delta=1000)
    operation.__apply__(orbit)
    assert orbit.pericenter_radius == 8000

# Test SetInclinationTo
def test_set_inclination_to():
    orbit = MockOrbit()
    operation = SetInclinationTo(inclination=30)
    operation.__apply__(orbit)
    assert orbit.inclination == 30

# Test ChangeInclanationBy
def test_change_inclination_by():
    orbit = MockOrbit(inclination=10)
    operation = ChangeInclanationBy(delta=20)
    operation.__apply__(orbit)
    assert orbit.inclination == 30

# Test PropagateAnomalyTo
def test_propagate_anomaly_to():
    orbit = MockOrbit()
    operation = PropagateAnomalyTo(M=pi)
    operation.__apply__(orbit)
    assert orbit.M == pi

# Test PropagateAnomalyBy
def test_propagate_anomaly_by():
    orbit = MockOrbit()
    operation = PropagateAnomalyBy(M=pi)
    operation.__apply__(orbit)
    assert orbit.M == pi

# Test Circularise
def test_circularise():
    orbit = MockOrbit(apocenter_radius=10000, pericenter_radius=7000)
    operation = Circularise(raise_pericenter=True)
    operation.__apply__(orbit)
    assert orbit.apocenter_radius == orbit.pericenter_radius

# Test SetPerincenterHere
def test_set_perincenter_here():
    orbit = MockOrbit()
    operation = SetPerincenterHere()
    operation.__apply__(orbit)
    assert orbit.arg_pe == orbit.f

# Test Maneuver
def test_maneuver():
    orbit = MockOrbit()
    maneuver = Maneuver(
        PropagateAnomalyTo(M=pi),
        ChangeApocenterBy(delta=2000),
        Circularise(raise_pericenter=True)
    )
    maneuver.__apply__(orbit)
    assert orbit.apocenter_radius == orbit.pericenter_radius