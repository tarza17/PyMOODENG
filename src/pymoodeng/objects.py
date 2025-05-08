
from . import constants
from . import anomaly
import numpy as np

def current_radius(r_p, e, theta):
  """
  Calculate current radius from perihelion radius, eccentricity, and true anomaly.
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

  def update(self, t, c_x = 0.0 , c_y = 0.0, logarithmic = False, depth = 1):
        # Check if orbit exists before trying to use it
        if self.orbit is None:
             self.x = c_x
             self.y = c_y
             return

        mean_anomaly = t * 2 * np.pi / self.orbit.tp # Mean anomaly is always 0 at t=0 (-> no addition of it here)
        true_anomaly = anomaly.Anomaly(mean_anomaly, type = 'M').theta(self.orbit.e)
        r = current_radius(self.orbit.peri_r,self.orbit.e, true_anomaly)
        if logarithmic:
            r = np.log10(r) / depth

        self.x = r * np.cos(true_anomaly) # Local x coordinate
        self.x = self.x + c_x # Global x coordinate
        self.y = r * np.sin(true_anomaly)    # Local y coordinate
        self.y = self.y + c_y # Global y coordinate

  def get_bodies(self):
        return [self]

class System:
    _instances = []
    def __init__(self, name, center, orbiting=None):
        """
        Initializes an orbital System.

        Args:
            center (Body): The central body of the system.
            bodies (list of Body, optional): List of bodies orbiting the center. Defaults to None.
        """
        self.__class__._instances.append(self)
        self.name = name
        self.center = center
        # self.orbiting_objects = orbiting if orbiting is not None else []
        if orbiting is not None:
            if isinstance(orbiting, (System, Body)):
                self.orbiting_objects= [orbiting]
            if isinstance(orbiting, (list, tuple)) and all(isinstance(item, (System, Body)) for item in orbiting):
                self.orbiting_objects = orbiting
        else:
            self.orbiting_objects = []

#        self.all_bodies = [center]

        # for body in self.orbiting_objects:
        #     self.all_bodies.append(body.get_bodies())
    @classmethod
    def get_all_instances(cls):
        # Return the list of all instances
        return cls._instances

    def add_body_or_system(self, body):
        self.orbiting_objects.append(body)

    # def add_system(self, body):
    #     self.bodies.append(body)

    def give_center_body(self, center_body):
        self.center = center_body

    def set_center_pos(self, coords):
        self.center.x = coords[0]
        self.center.y = coords[1]
    '''
    def update(self, t, c_x, c_y, logarithmic = False, depth = 0):
        if self.center.orbit:
          self.center.update(t, c_x, c_y, logarithmic)
        else:
          self.center.x = c_x
          self.center.y = c_y
        depth = depth + 1
        for obj in self.orbiting_objects:
          obj.update(t, self.center.x, self.center.y, logarithmic, depth)
    '''

    def update(self, t, c_x, c_y, logarithmic=False, depth=0):
        # --- MODIFIED LOGIC for positioning the center ---
        # Case 1: This system is a TOP-LEVEL system called by Universe (c_x, c_y are 0,0)
        #         AND its center has an orbit defined (e.g., Earth in EarthMoon sim).
        #         Force the center to stay at (0,0) for this context.
        if depth == 0 and self.center.orbit:
            self.center.x = 0.0
            self.center.y = 0.0
            # Do NOT call self.center.update() in this specific case.

        # Case 2: The center body has an orbit defined AND it's NOT the specific case above
        #         (i.e., this system is orbiting an outer center where c_x, c_y != 0,0).
        #         Update the center's position based on its orbit relative to c_x, c_y.
        elif self.center.orbit:
            # Update center relative to the outer center c_x, c_y
            # NOTE: Pass depth to Body.update if it uses it (it doesn't currently, but good practice)
            self.center.update(t, c_x, c_y, logarithmic, depth)

        # Case 3: The center body has NO orbit defined.
        #         Place it exactly at the reference point c_x, c_y.
        else:
            self.center.x = c_x
            self.center.y = c_y
        # --- END MODIFIED LOGIC ---

        # Now update all orbiting objects relative to the *final* position of the center
        orbiting_center_x = self.center.x
        orbiting_center_y = self.center.y
        depth = depth + 1  # Increment depth for subsystems
        for obj in self.orbiting_objects:
            # Pass the calculated center position and incremented depth
            obj.update(t, orbiting_center_x, orbiting_center_y, logarithmic, depth)

    def get_bodies(self):
       # return [self.center].append(self.orbiting_objects)
       bodies = [self.center]
       for obj in self.orbiting_objects:
           bodies.extend(obj.get_bodies())
       return bodies

class Universe:
    def __init__(self, systems=None):
        """
        Initializes simulation space.

        Args:
            systems : List of systems in the simulation.
        """

        # if systems is None:
        #     systems = []
        if systems is not None:
            if isinstance(systems, (System, Body)):
                systems = [systems]
            if isinstance(systems, (list, tuple)) and all(isinstance(item, (System, Body)) for item in systems):
                systems = systems
        else:
            systems = []

        self.all_elements = systems
        self.all_bodies = []
        for elem in self.all_elements:
            self.all_bodies.extend(elem.get_bodies())

    def update(self, time, logarithmic = False):
        for elem in self.all_elements:
          elem.update(time, 0.0, 0.0, logarithmic)

    def get_bodies(self):
        return self.all_bodies

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

Moon2 = Body(
    name = "Moon2",
    color = "grey",
    mass = constants.Moon_m,
    mean_diameter = constants.Moon_rm,
    orbit = Orbit(constants.Moon_perihelion * 150, constants.Moon_T *2, constants.Moon_e * 0.2))

Moon3 = Body(
    name = "Moon3",
    color = "grey",
    mass = constants.Moon_m,
    mean_diameter = constants.Moon_rm * 0.6,
    orbit = Orbit(constants.Moon_perihelion * 50, constants.Moon_T, constants.Moon_e))

EarthMoon = System(name = "EarthMoon", center = Earth, orbiting = Moon)

SunEarthMoon = System(name="SunEarthMoon", center = Sun, orbiting = EarthMoon)

Solar_system = System(name="Solar_system", center = Sun, orbiting = [
    Mercury,
    Venus,
    EarthMoon,
    Mars,
    Jupiter,
    Saturn,
    Uranus,
    Neptune,
    Pluto,
    Eris
])

Solar_system0 = System(name="Solar_system0", center = Sun, orbiting = [
    Mercury,
    Venus,
    EarthMoon,
    Mars,
])

Dwarfs = System(name="Dwarfs", center = Sun, orbiting = [Pluto, Eris])