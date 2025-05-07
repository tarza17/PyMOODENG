"""
Orbital and Physical Constants for Solar System Bodies

This module defines key astronomical and planetary constants for the Sun, Moon,
and major Solar System objects, based on authoritative sources such as the IAU
2009 resolutions and NASA planetary fact sheets.

Constants include:
- Speed of light and universal gravitational constant
- Masses and gravitational parameters (μ) for the Sun, planets, and dwarf planets
- Equatorial, polar, and mean radii (in kilometers)
- Orbital eccentricities
- Perihelion distances (in kilometers)
- Orbital periods (in seconds)

Sources:
- Masses and gravitational constants: IAU 2009 recommendations  
https://iau-a3.gitlab.io/NSFA/IAU2009_consts.html
- Radii and orbital elements:  
https://orbital-mechanics.space/reference/planetary-parameters.html  
https://nssdc.gsfc.nasa.gov/planetary/factsheet/  
https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html

All distances are in kilometers, times in seconds, and masses in kilograms unless otherwise noted.

Intended for use in orbital mechanics simulations, ephemeris generation, and astrodynamics analysis.

Example:
    from solar_system_constants import Earth_mu, Mars_e, Jupiter_T
"""

#IAU constants (source for mass: https://iau-a3.gitlab.io/NSFA/IAU2009_consts.html)
#Equitorial radii (sources: https://orbital-mechanics.space/reference/planetary-parameters.html, https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html)

c = 299792458 #speed of light (m/s)

G = 6.67428e-11  #gravitational constant (m^3/kg*s^2)

S_M_P = 1.32712440041e20  #Solar mass parameter (m^3/s^2)
Sun_m = S_M_P / G #mass of Sun (kg)

Earth_mu = 3.986004415e14  #Geocentric gravitational constant (m^3/s^2)
Earth_m  = Earth_mu / G #mass of Earth (kg)
Earth_re = 6378.1366 #equatorial radius of Earth (km)
Earth_rp = 6356.8 #polar radius of Earth (km)
Earth_rm = 6371 #mean radius of Earth (km)
Earth_perihelion   = 147098290  #km perihelion distance from Sun
Earth_T   = 365.256 * 24 * 3600     # 31_558_118.4 seconds

Moon_rm = 1737.4           # mean radius (km)
Moon_m = 7.342e22       # mass in kilograms
Moon_e = 0.0549            # orbital eccentricity
Moon_perihelion = 363300 #km
Moon_T    = 27.3217 * 24 * 3600     # 2_360_591.68 seconds

#masses (kg)
Mercury_m = Sun_m / 6.02365733e6
Venus_m   = Sun_m / 4.08523719e5
Mars_m    = Sun_m / 3.09870359e6
Jupiter_m = Sun_m / 1.047348644e3
Saturn_m  = Sun_m / 3.4979018e3
Uranus_m  = Sun_m / 2.2902951e4
Neptune_m = Sun_m / 1.941226e4
Pluto_m   = Sun_m / 1.3605e8
Eris_m    = Sun_m / 1.191e8

#gravitational constants (m^3/s^2)
Mercury_mu = Mercury_m * G
Venus_mu   = Venus_m   * G
Mars_mu    = Mars_m    * G
Jupiter_mu = Jupiter_m * G
Saturn_mu  = Saturn_m  * G
Uranus_mu  = Uranus_m  * G
Neptune_mu = Neptune_m * G
Pluto_mu   = Pluto_m * G
Eris_mu    = Eris_m * G

#equatorial radii (km)
Sun_re = 695700
Mercury_re = 2440.53
Venus_re   = 6051.8
Mars_re    = 3396.19
Jupiter_re = 71492
Saturn_re  = 60268
Uranus_re  = 25559
Neptune_re = 24764
Pluto_re   = 1188.3
Eris_re    = 1163

#polar radii (km)
Sun_rp = Sun_re
Mercury_rp = 2438.26
Venus_rp   = Venus_re
Mars_rp    = 3376.2
Jupiter_rp = 66854
Saturn_rp  = 54364
Uranus_rp  = 24973
Neptune_rp = 24341
Pluto_rp   = Pluto_re
Eris_rp    = Eris_re

#mean radii (km)
Sun_rm = Sun_rp
Mercury_rm = 2439.4
Venus_rm   = Venus_rp
Mars_rm    = 3389.5
Jupiter_rm = 69911
Saturn_rm  = 58232
Uranus_rm  = 25362
Neptune_rm = 24622
Pluto_rm   = Pluto_rp
Eris_rm    = Eris_rp

#perihelion radius
Mercury_perihelion = 46001200     # km
Venus_perihelion   = 107476000    # km
Mars_perihelion    = 206700000    # km
Jupiter_perihelion = 740573600    # km
Saturn_perihelion  = 1353572956   # km
Uranus_perihelion  = 2734998229   # km
Neptune_perihelion = 4452940833   # km
Pluto_perihelion   = 4436750000   # km
Eris_perihelion    = 5765740000   # km

#Orbital eccentricity
Mercury_e = 0.2056
Venus_e   = 0.0068
Earth_e   = 0.0167
Mars_e    = 0.0934
Jupiter_e = 0.0489
Saturn_e  = 0.0565
Uranus_e  = 0.0457
Neptune_e = 0.0113
Pluto_e   = 0.2488
Eris_e    = 0.4418

#Periodic time (s)
Mercury_T = 87.97 * 24 * 3600       # 7_600_128 seconds
Venus_T   = 224.70 * 24 * 3600      # 19_415_680 seconds
Mars_T    = 686.98 * 24 * 3600      # 59_355_072 seconds
Jupiter_T = 4332.59 * 24 * 3600     # 374_335_776 seconds
Saturn_T  = 10759.22 * 24 * 3600    # 929_596_608 seconds
Uranus_T  = 30688.5 * 24 * 3600     # 2_652_836_800 seconds
Neptune_T = 60182 * 24 * 3600       # 5_199_244_800 seconds
Pluto_T   = 90560 * 24 * 3600       # 7_824_384_000 seconds
Eris_T    = 203830 * 24 * 3600      # 17_613_312_000 seconds (approx)
