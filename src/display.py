import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib
import objects
matplotlib.use("TkAgg")

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