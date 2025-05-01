import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import matplotlib
import objects
matplotlib.use("TkAgg")

# Animation and simulation parameters
SECONDS_PER_DAY = 24 * 3600
START_DAY = 0
END_DAY = 730

# How many simulation days pass per animation frame
DAYS_PER_FRAME = 1

# Animation speed control (milliseconds between frames)
INTERVAL_MS = 20


# Creating the universe containing the desired system(s)
sim = objects.Universe([objects.Solar_system])

# Plotting setup
fig, ax = plt.subplots(figsize=(10,10))

# Background color (both figure and axes)
fig.set_facecolor('black')
ax.set_facecolor('black')

# Aspect ratio set to equal to avoid distortion
ax.set_aspect('equal', adjustable='box')

# Labels and title
ax.set_xlabel("X (m)", color="white")
ax.set_ylabel("Y (m)", color="white")
plot_title = ax.set_title("", color="white")
#ax.grid(True, linestyle='--', alpha=0.6)

# Adjust tick parameters
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')

# Adjust spine colors
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white')
ax.spines['left'].set_color('white')
ax.spines['right'].set_color('white')

# Scaling factor for body sizes
scale = 5000

# Global list to hold plot elements
plot_elements = []
bodies_to_plot = sim.get_bodies()

c_display_radius = 0.0
c_x, c_y = 0.0, 0.0

# Find the center body of the plotted system
center_obj = next((b for b in bodies_to_plot if b.name == objects.Solar_system.center.name), None)
if center_obj:
    # Calculate the radius the center body will be displayed with
    c_display_radius = (center_obj.mean_diameter / 2.0) * scale
    c_x = center_obj.x # Get the center's actual calculated position
    c_y = center_obj.y
else:
    print("Warning: Center body not found. Cannot calculate shift.")

def init():
    for patch in plot_elements:
        try:
            patch.remove()
        except ValueError:
            pass
    plot_elements.clear()

    # --- Run simulation for the very first frame (t=0) ---
    sim.update(START_DAY * SECONDS_PER_DAY, logarithmic=False)

    # --- Get initial positions and limits for t=0 ---
    min_x_init, max_x_init = np.inf, -np.inf
    min_y_init, max_y_init = np.inf, -np.inf

    # Get center's position at t=0
    c_x_init, c_y_init = (0.0, 0.0)
    if center_obj:
        c_x_init = center_obj.x
        c_y_init = center_obj.y

    # --- Create circle patches for t=0 and calculate initial limits ---
    for body in bodies_to_plot:
        x_orig = body.x # Position at t=0
        y_orig = body.y # Position at t=0
        color = body.color

        # Determine radius and initial plot position (with shift at t=0)
        if center_obj and body.name == center_obj.name:
            display_radius = c_display_radius
            z_order = 10
            plot_x, plot_y = x_orig, y_orig
        else:
            display_radius = (body.mean_diameter / 2.0) * scale
            z_order = 20
            # Shift calculation relative to center's t=0 position
            dx = x_orig - c_x_init
            dy = y_orig - c_y_init
            distance_from_center = np.sqrt(dx ** 2 + dy ** 2)
            if distance_from_center > 1e-9:
                ux = dx / distance_from_center
                uy = dy / distance_from_center
                plot_x = x_orig + ux * c_display_radius
                plot_y = y_orig + uy * c_display_radius
            else:
                plot_x, plot_y = x_orig, y_orig

        # Update limits based on initial positions/sizes
        min_x_init = min(min_x_init, plot_x - display_radius)
        max_x_init = max(max_x_init, plot_x + display_radius)
        min_y_init = min(min_y_init, plot_y - display_radius)
        max_y_init = max(max_y_init, plot_y + display_radius)

        # Create the circle at t=0 POSITION and add it
        circle = patches.Circle((plot_x, plot_y), radius=display_radius, color=color, zorder=z_order)
        ax.add_patch(circle)
        plot_elements.append(circle) # Store the circle object

    # --- Set initial axes limits for t=0 ---
    if np.isinf(min_x_init): min_x_init, max_x_init, min_y_init, max_y_init = -1e10, 1e10, -1e10, 1e10 # Fallback
    x_range = max_x_init - min_x_init
    y_range = max_y_init - min_y_init
    if x_range < 1e-9: x_range = abs(min_x_init) * 0.2 + 1e-9 if abs(min_x_init) > 1e-9 else 1e9
    if y_range < 1e-9: y_range = abs(min_y_init) * 0.2 + 1e-9 if abs(min_y_init) > 1e-9 else 1e9
    padding = 0.15
    ax.set_xlim(min_x_init - x_range * padding, max_x_init + x_range * padding)
    ax.set_ylim(min_y_init - y_range * padding, max_y_init + y_range * padding)

    # Set initial title
    plot_title.set_text(f"Solar System at t = {START_DAY:.1f} days")
    return plot_elements + [plot_title]


def animate(frame):
    # Calculate current simulation time based on frame number and step
    current_day = START_DAY + frame * DAYS_PER_FRAME
    current_time_seconds = current_day * SECONDS_PER_DAY

    # Update simulation state
    sim.update(current_time_seconds, logarithmic=False)

    # Variables to track plot extents for this frame
    min_x, max_x = np.inf, -np.inf
    min_y, max_y = np.inf, -np.inf

    # Get center's current position for shifting calculations in this frame
    c_x_now, c_y_now = (0.0,0.0)
    if center_obj:
        c_x_now = center_obj.x
        c_y_now = center_obj.y

    # Get the updates positions and calculate new axes limits
    for i, body in enumerate(bodies_to_plot):
        circle = plot_elements[i]  # Get the corresponding circle patch
        x_orig = body.x
        y_orig = body.y
        display_radius = circle.get_radius()

        if body.name == objects.Solar_system.center.name:
            plot_x, plot_y = x_orig, y_orig
        else:
            # Shift calculation relative to center body's current position
            dx = x_orig - c_x_now
            dy = y_orig - c_y_now
            distance_from_center = np.sqrt(dx ** 2 + dy ** 2)
            if distance_from_center > 1e-9:
                ux = dx / distance_from_center
                uy = dy / distance_from_center
                plot_x = x_orig + ux * c_display_radius
                plot_y = y_orig + uy * c_display_radius
            else:
                plot_x, plot_y = x_orig, y_orig

        # Update the circle's center position
        circle.set_center((plot_x, plot_y))

        # Update limits based on new positions
        min_x = min(min_x, plot_x - display_radius)
        max_x = max(max_x, plot_x + display_radius)
        min_y = min(min_y, plot_y - display_radius)
        max_y = max(max_y, plot_y + display_radius)

        # Update dynamic limits with padding
        if np.isinf(min_x): min_x, max_x, min_y, max_y = ax.get_xlim()[0], ax.get_xlim()[1], ax.get_ylim()[0], \
        ax.get_ylim()[1]  # Use previous if error
        x_range = max_x - min_x
        y_range = max_y - min_y
        # Prevent zero range if all objects are at one point
        if x_range < 1e-9: x_range = abs(min_x) * 0.2 + 1e-9
        if y_range < 1e-9: y_range = abs(min_y) * 0.2 + 1e-9
        padding = 0.15
        ax.set_xlim(min_x - x_range * padding, max_x + x_range * padding)
        ax.set_ylim(min_y - y_range * padding, max_y + y_range * padding)

        # Update the title
        plot_title.set_text(f"Solar System at t = {current_day:.1f} days")

        # Return iterable of plot elements that have been updated
    return plot_elements + [plot_title]

# Running the animation

# Calculate the number of frames needed
num_frames = int((END_DAY - START_DAY) / DAYS_PER_FRAME) + 1

print(f"Starting animation: {num_frames} frames, {DAYS_PER_FRAME} day(s) per frame, interval={INTERVAL_MS}ms")

# Create the animation object
ani = animation.FuncAnimation(fig, animate, frames=num_frames,
                              init_func=init, blit=False, interval=INTERVAL_MS,
                              repeat=True) # repeat=False stops after one loop

# Display the animation
plt.show()

print("Animation finished or window closed.")

#ani.save('solar_system_animation.mp4', fps=1000/INTERVAL_MS, dpi=150)
#print("Animation saved.")