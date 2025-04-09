#IAU constants (source for mass: https://iau-a3.gitlab.io/NSFA/IAU2009_consts.html)
#Equitorial radii (sources: https://orbital-mechanics.space/reference/planetary-parameters.html, https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html)

c = 299792458 #speed of light (m/s)

G = 6.67428e-11  #gravitational constant (m^3/kg*s^2)

S_M_P = 1.32712440041e20  #Solar mass parameter (m^3/s^2)
Sun_m = S_M_P / G #mass of Sun (kg)

Earth_mu = 3.986004415e14  #Geocentric gravitational constant (m^3/s^2)
Earth_m  = Earth_mu / G #mass of Earth (kg)
Eart_re = 6378.1366 #equatorial radius of Earth (km)
Earth_rp = 6356.8 #polar radius of Earth (km)
Earth_rm = 6371 #mean radius of Earth (km)

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