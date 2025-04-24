import constants
import anomaly
import numpy as np

def current_radius(r_p, e, theta):
  """
  Calculate current radius from percenter radius, eccentricity, and true anomaly.
  """
  r = r_p*(1 + e) / (1 + e * np.cos(theta))
  return r

class Orbit:
    def __init__(self, peri_r = 0, tp = 0, e = 0):
        """
        Contains orbital parameters.
        Args:
            e (float): Eccentricity of the orbital (0-1).
            peri_r (float): Perihelion radius.
            tp (float): Period time.
        """
        self.e = e # Eccentricity (0-1)
        self.peri_r = peri_r # Radius of perihelion (initial position)
        self.tp = tp # Period time

class Body:
  def __init__(self, name, color, mass, mean_diameter, orbit = None):
        """
        Initializes a Body object for 2D simulation.
        Args:
            name (str): Name of the body.
            mass (float): Mass of the body (kg).
            mean_diameter (float): Mean diameter of the body (m).
            color (str): Color for plotting.
            orbit (Orbit): Orbital parameters.
        """
        self.name = name
        self.mass = mass
        self.mean_diameter = mean_diameter
        self.color = color
        self.orbit = orbit
        if orbit is not None:
          self.x = orbit.peri_r # Initial x-coordinate
        else:
          self.x = 0.0  # Initial x-coordinate
        self.y = 0.0  # Initial y-coordinate

  def update(self, t, c_x = 0.0 , c_y = 0.0):
        mean_anomaly = t * 2 * np.pi / self.orbit.tp # Mean anomaly is always 0 at t=0 (-> no addition of it here)
        true_anomaly = anomaly.Anomaly(mean_anomaly, type = 'M').theta(self.orbit.e)
        r = current_radius(self.orbit.peri_r,self.orbit.e, true_anomaly)
        self.x = r * np.cos(true_anomaly) # Local x coordinate
        self.x = self.x + c_x # Global x coordinate
        self.y = r * np.sin(true_anomaly)    # Local y coordinate
        self.y = self.y + c_y # Global y coordinate

class System:
    def __init__(self, center, orbiting=None):
        """
        Initializes an orbital System.

        Args:
            center (Body): The central body of the system.
            bodies (list of Body, optional): List of bodies orbiting the center. Defaults to None.
        """
        self.center = center
        self.orbiting_objects = orbiting if orbiting is not None else []

    def add_body_or_system(self, body):
        self.orbiting_objects.append(body)

    # def add_system(self, body):
    #     self.bodies.append(body)

    def give_center_body(self, center_body):
        self.center = center_body

    def set_center_pos(self, coords):
        self.center.x = coords[0]
        self.center.y = coords[1]

    def update(self, t, c_x, c_y):
        if self.center.orbit:
          self.center.update(t, c_x, c_y)
        else:
          self.center.x = c_x
          self.center.y = c_y
        for obj in self.orbiting_objects:
          obj.update(t, self.center.x, self.center.y)

class Universe:
    def __init__(self, systems=None):
        """
        Initializes simulation space.

        Args:
            systems : List of systems in the simulation.
        """
        if systems is None:
            systems = []
        self.all_elements = systems

    def update(self, time):
        for elem in self.all_elements:
          elem.update(time, 0.0, 0.0)

Sun = Body(
    name = "Sun",
    color = "yellow",
    mass = constants.Sun_m,
    mean_diameter = constants.Sun_rm)

Earth = Body(
    name = "Earth",
    color = "blue",
    mass = constants.Earth_m,
    mean_diameter = constants.Earth_rm,
    orbit = Orbit(constants.Earth_perihelion, constants.Earth_T, constants.Earth_e))

Moon = Body(
    name = "Moon",
    color = "grey",
    mass = constants.Moon_m,
    mean_diameter = constants.Moon_rm,
    orbit = Orbit(constants.Moon_perihelion, constants.Moon_T, constants.Moon_e))

Mercury = Body(
    name="Mercury",
    color="darkgrey",
    mass=constants.Mercury_m,
    mean_diameter=constants.Mercury_rm,
    orbit=Orbit(constants.Mercury_perihelion, constants.Mercury_T, constants.Mercury_e)
)

Venus = Body(
    name="Venus",
    color="goldenrod",
    mass=constants.Venus_m,
    mean_diameter=constants.Venus_rm,
    orbit=Orbit(constants.Venus_perihelion, constants.Venus_T, constants.Venus_e)
)

Mars = Body(
    name="Mars",
    color="red",
    mass=constants.Mars_m,
    mean_diameter=constants.Mars_rm,
    orbit=Orbit(constants.Mars_perihelion, constants.Mars_T, constants.Mars_e)
)

Jupiter = Body(
    name="Jupiter",
    color="orange",
    mass=constants.Jupiter_m,
    mean_diameter=constants.Jupiter_rm,
    orbit=Orbit(constants.Jupiter_perihelion, constants.Jupiter_T, constants.Jupiter_e)
)

Saturn = Body(
    name="Saturn",
    color="khaki",
    mass=constants.Saturn_m,
    mean_diameter=constants.Saturn_rm,
    orbit=Orbit(constants.Saturn_perihelion, constants.Saturn_T, constants.Saturn_e)
)

Uranus = Body(
    name="Uranus",
    color="lightblue",
    mass=constants.Uranus_m,
    mean_diameter=constants.Uranus_rm,
    orbit=Orbit(constants.Uranus_perihelion, constants.Uranus_T, constants.Uranus_e)
)

Neptune = Body(
    name="Neptune",
    color="deepskyblue",
    mass=constants.Neptune_m,
    mean_diameter=constants.Neptune_rm,
    orbit=Orbit(constants.Neptune_perihelion, constants.Neptune_T, constants.Neptune_e)
)

Pluto = Body(
    name="Pluto",
    color="lightgrey",
    mass=constants.Pluto_m,
    mean_diameter=constants.Pluto_rm,
    orbit=Orbit(constants.Pluto_perihelion, constants.Pluto_T, constants.Pluto_e)
)

Eris = Body(
    name="Eris",
    color="white",
    mass=constants.Eris_m,
    mean_diameter=constants.Eris_rm,
    orbit=Orbit(constants.Eris_perihelion, constants.Eris_T, constants.Eris_e)
)

Solar_system = System(center = Sun, orbiting = [
    Mercury,
    Venus,
    Earth,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    Eris])