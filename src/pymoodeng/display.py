import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
from . import objects
matplotlib.use("TkAgg")


"""
Solar System Orbital Animation and Simulation Script

This script initializes a planetary system and runs a basic update of orbital states 
using a custom `Universe` class and object definitions from the `objects` module. 
It includes optional animated visualization using `matplotlib.animation`, though it is 
commented out in the current version.

Features:
- Imports planetary data from the `objects` module (e.g., Solar_system).
- Optionally animates a test plot of circular motion for visualization/debugging.
- Initializes a Universe simulation and updates planetary positions over a 1-day timestep.
- Prints updated x, y positions of all celestial bodies after simulation.

Modules:
- `numpy`: numerical operations
- `matplotlib`: plotting and animation
- `objects`: user-defined module containing planetary body and system definitions

Note:
- The animation section is currently commented out but can be enabled for basic orbital plotting.
- The script uses `TkAgg` backend for interactive plotting—ensure Tkinter is installed.

Example:
    $ python simulate.py
    Earth   1.4709829e+08   0.0
    Mars    2.067e+08       0.0
    ...

Requirements:
- matplotlib
- numpy
- A valid `objects.py` module with `Universe` and `Solar_system` definitions.

"""

# fig, ax = plt.subplots(figsize=(10, 10))
# ax.set_xlim(-1.5, 1.5)
# ax.set_ylim(-1.5, 1.5)
#
# test_plot, = ax.plot(0, 0, 'o', color = 'red')

# def init():
#     test_plot.set_data([0], [0])
#     return test_plot,
#
# def animate(i):
#     test_plot.set_data([np.cos(i / 360 * np.pi)], [np.sin(i / 360 * np.pi)])
#     return test_plot,
#
# anim = animation.FuncAnimation(fig, animate, init_func=init, frames = 720, interval = 1)

plt.show()

sim = objects.Universe([objects.Solar_system])

sim.update(1*24*3600)

for o in sim.get_bodies():
    print(o.name,' ', o.x, ' ', o.y)
        #print(type(ob))
        #print(ob.x, ob.y)

# for i in range (0,3650):
#     sim.update(i)
#     for ob in sim.get_bodies():
#         print(i)
#         for l in ob:
#             print(l.name,' ', l.x, ' ', l.y)
#         #print(type(ob))
#         #print(ob.x, ob.y)