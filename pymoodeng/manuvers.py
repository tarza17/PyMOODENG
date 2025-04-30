import warnings
from copy import copy

from numpy import allclose as almost_equal , divmod
from scipy.constants import pi

from . import orbital_conversions as oc
from .orbital_conversions import (
    elements_4_apsides, radius_from_alt, alt_from_radius, impulse_from_finite, saved_state)

class ReprHelperMixin:
    """Mixin class to assist with object representation."""
    def __repr__(self):
        return f"<{self.__class__.__name__}>"

from . import anomaly as anom
from .anomaly import (
     mean_anomaly_from_eccentric, mean_anomaly_from_true)   


########################################################################
#ez most máshonnan veszem 
########################################################################
_copy = copy

class Operation(object):

    def plot(self, orbit, plotter, next_operation=None):
        """Plot the operation on the given orbit using the provided plotter."""
        if hasattr(self, '__plot__') and callable(getattr(self, '__plot__')):
            self.__plot__(orbit, plotter, next_operation)

    def __add__(self, other):
        """Add two operations together."""
        if isinstance(other, Operation):
            return Maneuver(self, other)
        else:
            return NotImplementedError()
        
    def __radd__(self, other):
        """Right add two operations together."""
        if isinstance(other, Operation):
            return Maneuver(other, self)
        else:
            return NotImplementedError()

########################################################################
########################################################################
        
class ImpulseOperation(Operation):
    """Class representing an impulse operation."""

    def __init__(self):
        super(ImpulseOperation, self).__init__()

    def velocity(self):
        """Return the velocity of the impulse operation."""
        raise NotImplementedError()
    
########################################################################
########################################################################
    
class TimeOperation(Operation):
    """Class representing a time operation."""

    def __init__(self, time):
        super(TimeOperation, self).__init__()
        self.time = time

    def time(self):
        """Return the time of the operation."""
        raise NotImplementedError()
    
########################################################################
#APOCENTER
########################################################################


class SetApocenterRadiusTo(ReprHelperMixin, ImpulseOperation): 
    """Class representing an operation to set the apocenter radius."""

    def __init__(self, apocenter_radius, **kwargs):
        super(SetApocenterRadiusTo, self).__init__()
        self.apocenter_radius = apocenter_radius
        #self.kwargs = kwargs

    def __repr__(self):
        return f"SetApocenterRadiusTo(radius={self.radius}, kwargs={self.kwargs})"
    
    def __apply__(self, orbit):
        a , e = elements_4_apsides(self.apocenter_radius, orbit.pericenter_radius)

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        if orbit.apocenter_radius > self.apocenter_radius:
            label = 'Apocenter radius decreased'
        else:
            label = 'Apocenter radius increased'
        self.__apply__(orbit)
       
        """Plot the operation on the given orbit using the provided plotter."""
        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter.plot_apsides(orbit ,f1=0,  f2=f2 , label=label)


    def velocity(self, orbit):
        with saved_state(orbit):
            # get veloctiy and pericenter
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(self.apocenter_radius, orbit.pericenter_radius)

            orbit.semi_major_axis = a
            orbit.eccentricity = e

            new_v = orbit.v

        return new_v - old_v
      

    def represnetation_helper(self , r):
        r.positional_from_attr('apocenter_atlitude')

########################################################################
########################################################################

class SetApocenterAltitudeTo(ReprHelperMixin, ImpulseOperation):
    '''Set  apocenter alitude to a given value, when used must be at orbital pericenter'''

    def __init__(self, apocenter_altitude):
        super(SetApocenterAltitudeTo, self).__init__()
        self.apocenter_radius = apocenter_altitude

    def __apply__(self, orbit):
        apocenter_radius = orbit.body.mean_radius + self.apocenter_radius
        a, e = elements_4_apsides(radius_from_alt(self.apocenter_radius), orbit.pericenter_radius)

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        radius = radius_from_alt(self.apocenter_radius, orbit.body)
        if orbit.apocenter_radius > radius:
            label = 'Apocenter radius decreased'
        else:
            label = 'Apocenter radius increased'
        self.__apply__(orbit)

        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter._plot_apsides(orbit ,f1=0,  f2=f2 , label=label)

    def velocity(self, orbit):
        with saved_state(orbit):
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(radius_from_alt(self.apocenter_radius), orbit.pericenter_radius)

            orbit.a = a
            orbit.e = e

            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('apocenter_altitude')

########################################################################
########################################################################

class ChangeApocenterBy(ReprHelperMixin, ImpulseOperation): 
    '''Operation for changing apocenter radius, when applied must be at orpbital pericenter.'''
    
    def __init__(self, delta):
        super(ChangeApocenterBy, self).__init__()
        self.delta = delta

    def __apply__(self, orbit):
        a, e = elements_4_apsides(orbit.apocenter_radius + self.delta, orbit.pericenter_radius)
    
        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        if self.delta > 0:
            label = 'Apocenter radius increased'
        else:
            label = 'Apocenter radius decreased'
        self.__apply__(orbit)

        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter._plot_apsides(orbit ,f1=0,  f2=f2 , label=label)

    def velocity(self, orbit):
        with saved_state(orbit):
            # get veloctiy and pericenter
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(orbit.apocenter_radius + self.delta, orbit.pericenter_radius)

            orbit.a = a
            orbit.e = e

            new_v = orbit.v

        return new_v - old_v

    def __repr__(self, r):
        r.positional_from_attr('delta')

########################################################################
#PERICENTER
########################################################################

class SetPericenterRadiusTo(ReprHelperMixin, ImpulseOperation):
    '''Set pericenter radius to a given value, when applied must be at orbital apocenter.'''


    def __init__(self, pericenter_radius):
        super(SetPericenterRadiusTo, self).__init__()
        self.pericenter_radius = pericenter_radius

    def __apply__(self, orbit):
        a, e = elements_4_apsides(orbit.apocenter_radius, self.pericenter_radius)

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        if orbit.pericenter_radius > self.pericenter_radius:
            label = 'Pericenter radius decreased'
        else:
            label = 'Pericenter radius increased'
        self.__apply__(orbit)

        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter._plot_apsides(orbit ,f1=0,  f2=f2 , label=label)


    def velocity(self, orbit):
        with saved_state(orbit):
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(self.pericenter_radius, orbit.apocenter_radius)

            orbit.a = a
            orbit.e = e

            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('pericenter_radius')

########################################################################
########################################################################

class SetPericenterAltitudeTo(ReprHelperMixin, ImpulseOperation):
    '''Set pericenter alitude to a given value, when used must be at orbital apocenter'''

    def __init__(self, pericenter_altitude):
        super(SetPericenterAltitudeTo, self).__init__()
        self.pericenter_radius = pericenter_altitude

    def __apply__(self, orbit):
        pericenter_radius = orbit.body.mean_radius + self.pericenter_altitude
        a, e = elements_4_apsides(orbit.apocenter_radius, radius_from_alt(self.pericenter_radius))

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        radius = radius_from_alt(self.pericenter_radius, orbit.body)
        if orbit.pericenter_radius > radius:
            label = 'Pericenter radius decreased'
        else:
            label = 'Pericenter radius increased'
        self.__apply__(orbit)

        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter._plot_apsides(orbit ,f1=0,  f2=f2 , label=label)

    def velocity(self, orbit):
        with saved_state(orbit):
            # get veloctiy and pericenter
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(radius_from_alt(self.pericenter_radius), orbit.apocenter_radius)

            orbit.a = a
            orbit.e = e

            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('pericenter_altitude')


class ChangePericenterBy(ReprHelperMixin, ImpulseOperation):
    '''Operation for changing pericenter radius, when applied must be at orbital apocenter.'''
    
    def __init__(self, delta):
        super(ChangePericenterBy, self).__init__()
        self.delta = delta

    def __apply__(self, orbit):
        a, e = elements_4_apsides(orbit.apocenter_radius, orbit.pericenter_radius + self.delta)

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        if self.delta > 0:
            label = 'Pericenter radius increased'
        else:
            label = 'Pericenter radius decreased'
        self.__apply__(orbit)

        with saved_state(orbit):
            if (next_operation is not None and
                isinstance(next_operation, TimeOperation)):
                orbit.apply_maneuver(next_operation)
                f2 = orbit.theta
                if f2  == 0:
                    f2 = 2 * pi
                else:
                    f2 = 2 * pi

        plotter._plot_apsides(orbit ,f1=0,  f2=f2 , label=label)

    def velocity(self, orbit):
        with saved_state(orbit):
            # get veloctiy and pericenter
            orbit.propagate_anomaly_to(M=0)
            old_v = orbit.v

            a, e = elements_4_apsides(orbit.apocenter_radius, orbit.pericenter_radius + self.delta)

            orbit.a = a
            orbit.e = e

            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('delta')

########################################################################
#INCLANATION
########################################################################

class SetInclinationTo(ReprHelperMixin, ImpulseOperation):
    '''Set inclination to a given value, when used must be at ascending or descending node'''

    def __init__(self, inclination):
        super(SetInclinationTo, self).__init__()
        self.inclination = inclination

    def __apply__(self, orbit):
        orbit.i = self.inclination

    def __plot__(self, orbit, plotter, next_operation=None):
        self.__apply__(orbit)
        plotter._plot_apsides(orbit , label="Changed incalination")
    
    def velocity(self, orbit):
        with saved_state(orbit):
            orbit.f = 2 * pi - orbit.arg_pe 
            old_v = orbit.v
            
            self.__apply__(orbit)
            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('inclination')

#########################################################################
#########################################################################

class ChangeInclanationBy(ReprHelperMixin, ImpulseOperation):
    '''Operation for changing inclination, when applied must be at ascending or descending node'''

    def __init__(self, delta):
        super(ChangeInclanationBy, self).__init__()
        self.delta = delta

    def __apply__(self, orbit):
        orbit.i += self.delta

    def __plot__(self, orbit, plotter, next_operation=None):
        self.__apply__(orbit)
        plotter._plot_apsides(orbit , label="Changed incalination")
    
    def velocity(self, orbit):
        with saved_state(orbit):
            orbit.f = 2 * pi - orbit.arg_pe 
            old_v = orbit.v
            
            self.__apply__(orbit)
            new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.positional_from_attr('delta')

#########################################################################
#ANOMALY
#########################################################################

class PropagateAnomalyTo(ReprHelperMixin, TimeOperation):
    '''Propagate to a given time where the anomaly is equal to the given value
    ONLY 1 VALUE TO BE PASSED '''
    def __init__(self, **kwargs):
        super(PropagateAnomalyTo, self).__init__()

        valid_args = ['M', 'E', 'f']
        extra_args = set(kwargs.keys()) - valid_args

        # Check if any extra arguments are provided
        if extra_args:
            raise ValueError(f"Invalid arguments: {extra_args}. Valid arguments are: {valid_args}")
        
        # Check if there is a value at all 
        if not kwargs:
            raise ValueError("At least one anomaly value must be passed.")
        
        # Check if only 1 value is passed
        if sum([1 for anomaly in kwargs.values() if anomaly is not None]) != 1:
            raise ValueError("Only one anomaly value must be passed.")
        

#########################################################################
#########################################################################

def time(self, orbit):

    if self.key == 'M':
        M = self.anomaly
    if self.key == 'E':
        M = mean_anomaly_from_eccentric(orbit.e, self.anomaly)
    if self.key == 'theta':
        M = mean_anomaly_from_true(orbit.e, self.anomaly)

    # Propagate one orbit if dest is "passed"

    if M < orbit.M:
        M += 2 * pi

    return (M - orbit.M) / orbit.n

def __plot__(self, orbit, plotter, next_operation=None):
    f1 = orbit.f1
    orbit.t += self.time_delta(orbit)
    f2 = orbit.f

    plotter._plot_apsides(orbit , f2=f2 , propageted =True )

def __repr__(self, r):
    r.keyword_with_value(self.key. self.anomaly)

#########################################################################
#########################################################################

class PropagateAnomalyBy(ReprHelperMixin, TimeOperation):
    '''Propagate to by the given  value'''

    def __init__(self, **kwargs):
        super(PropagateAnomalyBy, self).__init__()

        valid_args = ['M', 'E', 'f']
        extra_args = set(kwargs.keys()) - valid_args

        # Check if any extra arguments are provided
        if extra_args:
            raise ValueError(f"Invalid arguments: {extra_args}. Valid arguments are: {valid_args}")
        
        # Check if there is a value at all 
        if not kwargs:
            raise ValueError("At least one anomaly value must be passed.")
        
        # Check if only 1 value is passed
        if sum([1 for anomaly in kwargs.values() if anomaly is not None]) != 1:
            raise ValueError("Only one anomaly value must be passed.")
        

    def time(self, orbit):
        if self.key == 'f':
            orbits, f = divmod(self.anomaly, 2 * pi)
            M = mean_anomaly_from_true(orbit.e, f)
            return orbits * orbit.T + M / orbit.n
        elif self.key == 'E':
            orbits, E = divmod(self.anomaly, 2 * pi)
            M = mean_anomaly_from_eccentric(orbit.e, E)
            return orbits * orbit.T + M / orbit.n
        elif self.key == 'M':
            return self.anomaly / orbit.n
        
    def __plot__(self, orbit, plotter, next_operation=None):
        f1 = orbit.f
        orbit.t += self.time_delta(orbit)
        f2 = orbit.f

        plotter._plot_apsides(orbit , f2=f2 , propageted = True )

    def __repr__(self, r):
        r.keyword_with_value(self.key, self.anomaly)

#########################################################################
#ETC
#########################################################################

class Circularise(ReprHelperMixin, ImpulseOperation):
    '''Circularise the orbit, must be at apsides'''

    def __init__(self, raise_pericenter=True):
        super(Circularise, self).__init__()
        self.raise_pericenter = raise_pericenter

    def __apply__(self, orbit):
        if self.raise_pericenter:
            radius = orbit.apocenter_radius
        else:
            radius = orbit.pericenter_radius

        a, e = elements_4_apsides(radius, radius)

        orbit.a = a
        orbit.e = e

        orbit.v = orbit.v

    def __plot__(self, orbit, plotter, next_operation=None):
        self.__apply__(orbit)
        plotter._plot_apsides(orbit , label="Circularised orbit")
    
    def velocity(self, orbit):
        with saved_state(orbit):
          if self.raise_pericenter:
            orbit.propagate_anomaly_to(M=pi)
            radius = orbit.apocenter_radius
          else:
            orbit.propagate_anomaly_to(M=0)
            radius = orbit.pericenter_radius

        old_v = orbit.v

        a, e = elements_4_apsides(radius, radius)
        orbit.a = a
        orbit.e = e
        new_v = orbit.v

        return new_v - old_v
    
    def __repr__(self, r):
        r.keyword_with_attr('raise_pericenter')

#########################################################################
#########################################################################

class SetPerincenterHere(ReprHelperMixin, ImpulseOperation):
    '''Set pericenter here to the exact location, prep for eliptical oribt, must be from circular orbit'''

    def __init__(self, raise_pericenter=True):
        super(SetPerincenterHere, self).__init__()

    def __apply__(self, orbit):
        orbit.arg_pe = orbit.f
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            orbit.f = 0.0
    
    def __plot__(self, orbit, plotter, next_operation=None):
        self.__apply__(orbit)
        plotter._plot_apsides(orbit , label="Set pericenter here")
    
    def __repr__(self, r):
        pass


#########################################################################
#MANEUVER (Manover? Maneuver? ManRollover with a french accent?)
#########################################################################

class Maneuver(ReprHelperMixin, object):
    ''''Collection of operations '''

    def __init__(self, *operations):
        if not isinstance(operations, (list)):
            operations = [operations]

        self.operations = operations
    
#APOCENTER MANROLLOVERS

    @classmethod
    def set_apocenter_radius_to(cls, apocenter_radius):
        """Maneuver to set the apocenter radius."""
        operations = [
            PropagateAnomalyTo(M=pi),
            SetApocenterRadiusTo(apocenter_radius),
        ]
        return cls(operations)
    
    @classmethod
    def set_apocenter_altitude_to(cls, apocenter_altitude):
        """Maneuver to set the apocenter altitude."""
        operations = [
            PropagateAnomalyTo(M=pi),
            SetApocenterAltitudeTo(apocenter_altitude),
        ]
        return cls(operations)
    
    @classmethod
    def change_apocenter_by(cls, delta):
        """Maneuver to change the apocenter radius."""
        operations = [
            PropagateAnomalyTo(M=pi),
            ChangeApocenterBy(delta),
        ]
        return cls(operations)
    
#PERICENTER MANROLLOVERS

    @classmethod
    def set_pericenter_radius_to(cls, pericenter_radius):
        """Maneuver to set the pericenter radius."""
        operations = [
            PropagateAnomalyTo(M=0),
            SetPericenterRadiusTo(pericenter_radius),
        ]
        return cls(operations)
    
    @classmethod
    def set_pericenter_altitude_to(cls, pericenter_altitude):
        """Maneuver to set the pericenter altitude."""
        operations = [
            PropagateAnomalyTo(M=0),
            SetPericenterAltitudeTo(pericenter_altitude),
        ]
        return cls(operations)
    
    @classmethod
    def change_pericenter_by(cls, delta):
        """Maneuver to change the pericenter radius."""
        operations = [
            PropagateAnomalyTo(M=0),
            ChangePericenterBy(delta),
        ]
        return cls(operations)
    
#INCLINATION MANROLLOVERS

    @classmethod
    def set_inclination_to(cls, inclination):
        """Maneuver to set the inclination."""
        operations = [
            lambda orbit:
            PropagateAnomalyTo(f=2 *pi - orbit.arg_pe),
            SetInclinationTo(inclination),
        ]
        return cls(operations)

    @classmethod
    def change_inclination_by(cls, delta):
        """Maneuver to change the inclination."""
        operations = [
            lambda orbit:
            PropagateAnomalyTo(f=2 *pi - orbit.arg_pe),
            ChangeInclanationBy(delta),
        ]
        return cls(operations)
    
#HOHHMANN MANROLLOVERS
    
    @classmethod
    def hohmann_to_radius(cls, radius):
        """Maneuver to perform a Hohmann transfer to change the apocenter from pos."""
        operations = [
            SetPerincenterHere(),
            SetApocenterRadiusTo(radius),
            PropagateAnomalyTo(M=pi),
            Circularise(),
        ]
        return cls(operations)
    
    @classmethod
    def hohmann_to_altitude(cls, altitude):
        """Maneuver to perform a Hohmann transfer to change the apocenter from pos."""
        operations = [
            SetPerincenterHere(),
            SetApocenterAltitudeTo(altitude),
            PropagateAnomalyTo(M=pi),
            Circularise(),
        ]
        return cls(operations)
    


    def __apply__(self, orbit):
        for operation in self.operations:
            if callable(operation):
                operation = operation(orbit)
            if hasattr(operation, '__apply__') and callable(getattr(operation, '__apply__')):
                operation.__apply__(orbit)
            elif isinstance(operation, ImpulseOperation):
                orbit.v += operation.velocity(orbit)
            elif isinstance(operation, TimeOperation):
                orbit.t += operation.time(orbit)

    def __iapply__(self,orbit, copy=False):
        for operation in self.operations:
            if callable(operation):
                operation = operation(orbit)
            
            if copy:
                yield _copy(orbit), operation
            else:
                yield orbit, operation
            
            if hasattr(operation, '__apply__') and callable(getattr(operation, '__apply__')):
                operation.__apply__(orbit)
            elif isinstance(operation, ImpulseOperation):
                orbit.v += operation.velocity(orbit)
            elif isinstance(operation, TimeOperation):
                orbit.t += operation.time(orbit)

    def __repr__(self, r):
        r.positional_from_attr('operations')

    def __add__(self, other):
        """Add two maneuvers together."""
        if isinstance(other, Maneuver):
            return Maneuver(self.operations + other.operations)
        else:
            return NotImplementedError()