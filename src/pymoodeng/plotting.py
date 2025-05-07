import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.widgets import Slider, Button, RadioButtons
import matplotlib
import inspect
from objects import System
import objects
from collections import deque

matplotlib.use("TkAgg")

def get_system_instances():
    # Getting all members from the objects module
    all_members = inspect.getmembers(objects)
    # Filtering for instances of System
    system_instances = [member for name, member in all_members if isinstance(member, System)]

    return system_instances

available_systems = get_system_instances()
system_names = [sys.name for sys in available_systems]


# Creating the universe containing the desired system(s)
sim = objects.Universe([objects.Solar_system])
selected_system = objects.Solar_system

# Interaction State
user_interacting = False
pan_info = {'button': None, 'start_x': None, 'start_y': None, 'current_xlim': None, 'current_ylim': None}

# Animation and simulation parameters
SECONDS_PER_DAY = 24 * 3600
START_DAY = 0
END_DAY = 365*300

# Animation speed control parameters
PAUSE = False
INTERVAL_MS = 20
DAYS_PER_FRAME = 1
CURRENT_FRAME = 0

# Orbital Trail Settings
TRAIL_VISIBLE = True
MAX_TRAIL_POINTS = None
TRAIL_COLOR = 'white'
TRAIL_LINEWIDTH = 0.5
TRAIL_ALPHA = 0.7

# Plotting setup
fig, ax = plt.subplots(figsize=(15,15))
NAMES_VISIBLE = False

# Background color (both figure and axes)
fig.set_facecolor('black')
ax.set_facecolor('black')

# Aspect ratio set to equal to avoid distortion
ax.set_aspect('equal', adjustable='box')

# Labels and title
ax.set_xlabel("X (m)", color="white")
ax.set_ylabel("Y (m)", color="white")
plot_title = ax.set_title("", color="white")

# Adjust tick parameters
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')

# Adjust spine colors
ax.spines['bottom'].set_color('white')
ax.spines['top'].set_color('white')
ax.spines['left'].set_color('white')
ax.spines['right'].set_color('white')

# Scaling factors for center body and other body sizes
c_scale = 1
scale = 1


# Global list to hold plot elements
plot_elements = []
trail_lines = []
trajectory_history = {}
plot_positions_this_frame = []

bodies_to_plot = sim.get_bodies()

# Create radio buttons for system selection
ax_systems = plt.axes([0.8, 0.3, 0.15, 0.15])  # Position for the system selector
radio_systems = RadioButtons(ax_systems, system_names)

# Variables for sticky axes
global_min_x, global_max_x = np.inf, -np.inf
global_min_y, global_max_y = np.inf, -np.inf
limits_initialized = False # Flag to track if initial limits are set
padding = 0.15

c_display_radius = 0.0
c_x, c_y = 0.0, 0.0
center_obj = None

# Heuristic to find the "main" center if multiple systems exist
# Prefers the center of the first system in the universe list
if sim.all_elements and isinstance(sim.all_elements[0], objects.System):
    center_obj_candidate = sim.all_elements[0].center
    # Verify it's actually in the list of bodies to plot
    center_obj = next((b for b in bodies_to_plot if b.name == center_obj_candidate.name), None)

if center_obj:
    # Calculate the radius the center body will be displayed with
    c_display_radius = (center_obj.mean_diameter / 2.0) * c_scale
    # We'll get its actual position during the t=0 update in init(), this is just needed for shifting
else:
    print("Warning: Center body for shifting not definitively identified. Shifting might be inconsistent.")

# Callback functions for sliders
def update_c_scale(val):
    global c_scale
    c_scale = val
    slider_c_scale.label.set_text(f'Center Body Scale  ({c_scale:.2f}): ')  # Update label text

    # Clear trail lines and history
    for line in trail_lines:
        line.set_data([], [])  # Clear the line data
    for i in range(len(trajectory_history)):
        trajectory_history[i].clear()  # Clear the history for each body

def update_scale(val):
    global scale
    scale = val
    slider_scale.label.set_text(f'Orbiting Bodies Scale ({scale:.2f}): ')  # Update label text


# Callback functions for buttons
def toggle_trail(event):
    global TRAIL_VISIBLE
    TRAIL_VISIBLE = not TRAIL_VISIBLE

def toggle_names(event):
    global NAMES_VISIBLE
    NAMES_VISIBLE = not NAMES_VISIBLE


# Create sliders
ax_c_scale = plt.axes([0.8, 0.14, 0.15, 0.03])  # Position for c_scale slider
slider_c_scale = Slider(ax_c_scale, ' ', 1, 10000.0, valinit=c_scale, color="white")
slider_c_scale.label.set_color('white')
slider_c_scale.label.set_text(f'Center Body Scale  ({c_scale:.2f}): ')
slider_c_scale.on_changed(update_c_scale)

ax_scale = plt.axes([0.8, 0.1, 0.15, 0.03])  # Position for scale slider
slider_scale = Slider(ax_scale, ' ', 1, 10000.0, valinit=scale, color="white")
slider_scale.label.set_color('white')
slider_scale.label.set_text(f'Orbiting Bodies Scale ({scale:.2f}): ')
slider_scale.on_changed(update_scale)

# Create buttons
ax_toggle_trail = plt.axes([0.8, 0.199, 0.15, 0.04])  # Position for toggle trail button
button_trail = Button(ax_toggle_trail, 'Trails Visibility')
button_trail.on_clicked(toggle_trail)

ax_toggle_names = plt.axes([0.8, 0.249, 0.15, 0.04])  # Position for toggle names button
button_names = Button(ax_toggle_names, 'Names Visibility')
button_names.on_clicked(toggle_names)

def apply_global_limits_with_padding():
    # Applies the current global limits (with padding) to the axes.
    global global_min_x, global_max_x, global_min_y, global_max_y
    if not limits_initialized:
        return # Don't apply if not initialized

    # Calculate range based on the *global* min/max data points seen so far
    x_range = global_max_x - global_min_x
    y_range = global_max_y - global_min_y

    # Add a small absolute value or fraction of the coordinate if non-zero
    min_abs_range = 1e9 # Minimum range to prevent excessive zoom on single points
    if x_range < 1e-9:
        x_range = max(abs(global_min_x) * 0.2, min_abs_range) if abs(global_min_x) > 1e-9 else min_abs_range
    if y_range < 1e-9:
        y_range = max(abs(global_min_y) * 0.2, min_abs_range) if abs(global_min_y) > 1e-9 else min_abs_range

    # Calculate desired padded limits based on global range
    target_xmin = global_min_x - x_range * padding
    target_xmax = global_max_x + x_range * padding
    target_ymin = global_min_y - y_range * padding
    target_ymax = global_max_y + y_range * padding

    # Set the limits on the axes using the calculated global targets
    ax.set_xlim(target_xmin, target_xmax)
    ax.set_ylim(target_ymin, target_ymax)

# Event Handlers

def on_scroll(event):
    """Handle mouse scroll for zooming."""
    global user_interacting
    if event.inaxes != ax: return

    user_interacting = True # User has taken control
    base_scale = 1.1 # Zoom factor
    cur_xlim = ax.get_xlim()
    cur_ylim = ax.get_ylim()
    xdata = event.xdata # Mouse x position
    ydata = event.ydata # Mouse y position

    if event.button == 'up': # Zoom in
        scale_factor = 1 / base_scale
    elif event.button == 'down': # Zoom out
        scale_factor = base_scale
    else: # Unexpected event
        scale_factor = 1

    new_width = (cur_xlim[1] - cur_xlim[0]) * scale_factor
    new_height = (cur_ylim[1] - cur_ylim[0]) * scale_factor

    relx = (cur_xlim[1] - xdata) / (cur_xlim[1] - cur_xlim[0])
    rely = (cur_ylim[1] - ydata) / (cur_ylim[1] - cur_ylim[0])

    ax.set_xlim([xdata - new_width * (1 - relx), xdata + new_width * relx])
    ax.set_ylim([ydata - new_height * (1 - rely), ydata + new_height * rely])
    ax.figure.canvas.draw_idle() # Request redraw

def on_press(event):
    """Handle mouse button press for panning."""
    global user_interacting, pan_info
    if event.inaxes != ax: return
    # Use middle mouse button for panning (button 2), or left (button 1) if preferred
    if event.button != 1: return

    user_interacting = True # User has taken control
    pan_info['button'] = event.button
    pan_info['start_x'] = event.xdata
    pan_info['start_y'] = event.ydata
    pan_info['current_xlim'] = ax.get_xlim()
    pan_info['current_ylim'] = ax.get_ylim()

def on_motion(event):
    """Handle mouse motion while button is pressed for panning."""
    global pan_info
    if event.inaxes != ax or pan_info['button'] != event.button: return
    if pan_info['start_x'] is None or event.xdata is None: return # Avoid errors if mouse leaves axes

    dx = event.xdata - pan_info['start_x']
    dy = event.ydata - pan_info['start_y']

    new_xlim = (pan_info['current_xlim'][0] - dx, pan_info['current_xlim'][1] - dx)
    new_ylim = (pan_info['current_ylim'][0] - dy, pan_info['current_ylim'][1] - dy)

    ax.set_xlim(new_xlim)
    ax.set_ylim(new_ylim)
    ax.figure.canvas.draw_idle() # Request redraw

def on_release(event):
    """Handle mouse button release to stop panning."""
    global pan_info
    if event.button == pan_info['button']:
        pan_info = {'button': None, 'start_x': None, 'start_y': None, 'current_xlim': None, 'current_ylim': None}

def on_key(event):
    """Handle key press events, e.g., reset view."""
    global user_interacting, global_min_x, global_max_x, global_min_y, global_max_y, PAUSE, c_scale
    if event.key == 'r': # 'r' for reset
        user_interacting = False
        # Recalculate limits based on *current* positions to avoid jumping
        # This requires getting current positions, similar to animate()
        frame_min_x, frame_max_x = np.inf, -np.inf
        frame_min_y, frame_max_y = np.inf, -np.inf

        # Need current center position if shifting is enabled
        c_x_now, c_y_now = (0.0, 0.0)
        c_display_radius = (center_obj.mean_diameter / 2.0) * c_scale

        if center_obj:
            c_x_now = center_obj.x
            c_y_now = center_obj.y

        for i, body in enumerate(bodies_to_plot):
             x_orig = body.x
             y_orig = body.y
             # Use stored radius if possible, otherwise recalculate
             try:
                 display_radius = plot_elements[i].get_radius()
             except (IndexError, AttributeError):
                 display_radius = (body.mean_diameter / 2.0) * scale # Fallback

             plot_x, plot_y = x_orig, y_orig
             if center_obj and body.name != center_obj.name:
                 dx = x_orig - c_x_now
                 dy = y_orig - c_y_now
                 distance_from_center = np.sqrt(dx**2 + dy**2)
                 if distance_from_center > 1e-9 and c_display_radius > 1e-9:
                     ux = dx / distance_from_center
                     uy = dy / distance_from_center
                     plot_x = x_orig + ux * c_display_radius
                     plot_y = y_orig + uy * c_display_radius

             frame_min_x = min(frame_min_x, plot_x - display_radius)
             frame_max_x = max(frame_max_x, plot_x + display_radius)
             frame_min_y = min(frame_min_y, plot_y - display_radius)
             frame_max_y = max(frame_max_y, plot_y + display_radius)

        if not np.isinf(frame_min_x):
            global_min_x, global_max_x = frame_min_x, frame_max_x
            global_min_y, global_max_y = frame_min_y, frame_max_y
            apply_global_limits_with_padding() # Apply the reset limits
            ax.figure.canvas.draw_idle() # Request redraw
        else:
            print("Warning: Could not determine current limits for reset.")

    if event.key == ' ':
        PAUSE = not PAUSE
        if PAUSE:
            ani.event_source.stop()
        else:
            ani.event_source.start()

def init():

    global global_min_x, global_max_x, global_min_y, global_max_y, limits_initialized, user_interacting, c_scale, bodies_to_plot

    bodies_to_plot = sim.get_bodies()

    user_interacting = False

    for patch in plot_elements:
        try:
            patch.remove()
        except ValueError:
            pass
    plot_elements.clear()

    for line in trail_lines:
        try: line.remove()
        except ValueError: pass
    trail_lines.clear()

    trajectory_history.clear() # Clear history

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
        c_display_radius = (center_obj.mean_diameter / 2.0) * c_scale

    # --- Create circle patches for t=0 and calculate initial limits ---
    for i, body in enumerate(bodies_to_plot):
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

        # Initialize Trail
        trajectory_history[i] = deque(maxlen=MAX_TRAIL_POINTS)
        trajectory_history[i].append((plot_x, plot_y))
        line, = ax.plot([], [], color=TRAIL_COLOR, linewidth=TRAIL_LINEWIDTH, alpha=TRAIL_ALPHA,
                        zorder=z_order - 1)  # Trail behind body
        trail_lines.append(line)  # Store the line object

    # Set initial axes limits for t=0
    if np.isinf(min_x_init): min_x_init, max_x_init, min_y_init, max_y_init = -1e10, 1e10, -1e10, 1e10 # Fallback

    global_min_x, global_max_x = min_x_init, max_x_init
    global_min_y, global_max_y = min_y_init, max_y_init
    limits_initialized = True

    apply_global_limits_with_padding()

    '''
    if x_range < 1e-9: x_range = abs(min_x_init) * 0.2 + 1e-9 if abs(min_x_init) > 1e-9 else 1e9
    if y_range < 1e-9: y_range = abs(min_y_init) * 0.2 + 1e-9 if abs(min_y_init) > 1e-9 else 1e9
    padding = 0.15
    ax.set_xlim(min_x_init - x_range * padding, max_x_init + x_range * padding)
    ax.set_ylim(min_y_init - y_range * padding, max_y_init + y_range * padding)
    '''
    # Set initial title
    plot_title.set_text(f"{selected_system.name} centric system at t = {START_DAY:.1f} days ({START_DAY/365.2425:.2f} years)")
    return plot_elements + trail_lines + [plot_title]


def animate(frame):
    global CURRENT_FRAME, global_min_x, global_max_x, global_min_y, global_max_y, plot_positions_this_frame, scale, c_scale, NAMES_VISIBLE, DAYS_PER_FRAME

    # Update the frame counter
    CURRENT_FRAME = frame
    # Calculate current simulation time based on frame number and step
    current_day = START_DAY + frame * DAYS_PER_FRAME
    current_time_seconds = current_day * SECONDS_PER_DAY

    # Update simulation state
    sim.update(current_time_seconds, logarithmic=False)

    # Get center's current position for shifting calculations in this frame
    c_x_now, c_y_now = (0.0,0.0)
    if center_obj:
        c_x_now = center_obj.x
        c_y_now = center_obj.y
        c_display_radius = (center_obj.mean_diameter / 2.0) * c_scale

    # Clear previous frame's positions
    plot_positions_this_frame.clear()

    # Get the updates positions and calculate new axes limits
    for i, body in enumerate(bodies_to_plot):
        if i >= len(plot_elements):  # Safety check
            print(
                f"Warning: Mismatch between bodies_to_plot ({len(bodies_to_plot)}) and plot_elements ({len(plot_elements)}) at index {i}")
            continue
        circle = plot_elements[i]  # Get the corresponding circle patch
        line = trail_lines[i]

        x_orig = body.x
        y_orig = body.y

        if body.name == center_obj.name:
            display_radius = (body.mean_diameter / 2.0) * c_scale  # Center body
        else:
            display_radius = (body.mean_diameter / 2.0) * scale  # Orbiting bodies

        plot_x, plot_y = x_orig, y_orig
        if center_obj and body.name != center_obj.name:
            dx = x_orig - c_x_now
            dy = y_orig - c_y_now
            distance_from_center = np.sqrt(dx**2 + dy**2)
            if distance_from_center > 1e-9 and c_display_radius > 1e-9:
                ux = dx / distance_from_center
                uy = dy / distance_from_center
                plot_x = x_orig + ux * c_display_radius
                plot_y = y_orig + uy * c_display_radius

        # Update the circle's center position
        circle.set_center((plot_x, plot_y))
        circle.set_radius(display_radius)
        plot_positions_this_frame.append((plot_x, plot_y, display_radius)) # Store for limits/reset

        # --- Update Trail ---
        # Append current *plotted* position to this body's history
        trajectory_history[i].append((plot_x, plot_y))
        # Extract history for plotting (handle case of < 2 points)
        hist = list(trajectory_history[i])  # Convert deque to list for slicing/zipping
        if len(hist) >= 2:
            hist_x, hist_y = zip(*hist)
            line.set_data(hist_x, hist_y)
        elif len(hist) == 1:  # Only one point
            line.set_data([hist[0][0]], [hist[0][1]])  # Plot single point
        else:  # No history yet
            line.set_data([], [])

        # Update trail visibility
        if TRAIL_VISIBLE:
            line.set_visible(True)
        else:
            line.set_visible(False)

        # Update Names
        if NAMES_VISIBLE:
            # Check if the body is within axes limits
            if ax.get_xlim()[0] <= plot_x <= ax.get_xlim()[1] and ax.get_ylim()[0] <= plot_y <= ax.get_ylim()[1]:
                # Create or update the text label for the body
                if not hasattr(circle, 'name_text'):
                    # Create a new text object if it doesn't exist
                    circle.name_text = ax.text(plot_x, plot_y, body.name, color='white', fontsize=8, ha='center',
                                               va='bottom', zorder=30)
                else:
                    # Update the position of the existing text object
                    circle.name_text.set_position((plot_x, plot_y + display_radius))
                    circle.name_text.set_text(body.name)  # Update the name text
            else:
                # If names are not visible, remove any existing text
                if hasattr(circle, 'name_text'):
                    circle.name_text.remove()
                    del circle.name_text
        else:
            # If names are not visible, remove any existing text
            if hasattr(circle, 'name_text'):
                circle.name_text.remove()
                del circle.name_text

        # Handle Axes Limits
    if not user_interacting:
            # Calculate frame limits ONLY if not user interacting
        frame_min_x, frame_max_x = np.inf, -np.inf
        frame_min_y, frame_max_y = np.inf, -np.inf
        for px, py, radius in plot_positions_this_frame:
            frame_min_x = min(frame_min_x, px - radius)
            frame_max_x = max(frame_max_x, px + radius)
            frame_min_y = min(frame_min_y, py - radius)
            frame_max_y = max(frame_max_y, py + radius)

        if not np.isinf(frame_min_x):  # Ensure we got valid numbers this frame
            global_min_x = min(global_min_x, frame_min_x)
            global_max_x = max(global_max_x, frame_max_x)
            global_min_y = min(global_min_y, frame_min_y)
            global_max_y = max(global_max_y, frame_max_y)

            # Update dynamic limits with padding
            '''
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
            '''

        # Apply the potentially updated global limits with padding
        # This ensures the axes only expand or stay the same, never shrink
        apply_global_limits_with_padding()

    # Update the title
    plot_title.set_text(f"{selected_system.name} at t = {current_day:.1f} days ({current_day/365.2425:.2f} years)")

        # Return iterable of plot elements that have been updated
    return plot_elements + trail_lines + [plot_title]

# Callback function for changing systems
def switch_system(label):
    global sim, bodies_to_plot, plot_elements, trail_lines, trajectory_history, center_obj, selected_system, CURRENT_FRAME

    # Reset the frame counter
    CURRENT_FRAME = 0

    # Find the selected system by name
    selected_system = available_systems[system_names.index(label)]
    sim = objects.Universe([selected_system])  # Update the simulation with the new system

    # Clear all existing name texts from the axes
    for text in ax.texts:
        try:
            text.remove()
        except:
            pass

    # Reset plot elements and trajectory history
    for patch in plot_elements:
        try:
            # Remove any name text attached to the patch
            if hasattr(patch, 'name_text'):
                try:
                    patch.name_text.remove()
                except:
                    pass
            patch.remove()
        except ValueError:
            pass
    plot_elements.clear()

    for line in trail_lines:
        try:
            line.remove()
        except ValueError:
            pass
    trail_lines.clear()

    trajectory_history.clear()

    # Update bodies to plot based on the new system
    bodies_to_plot = sim.get_bodies()

    # Update center object for the new system
    if sim.all_elements and isinstance(sim.all_elements[0], objects.System):
        center_obj_candidate = sim.all_elements[0].center
        center_obj = next((b for b in bodies_to_plot if b.name == center_obj_candidate.name), None)
    else:
        center_obj = None

    # Reinitialize the simulation
    init()  # Call the init function to set up the new system

    # Restart the animation from frame 0
    ani.frame_seq = ani.new_frame_seq()
    print(f"Switched to system: {label}")


# Connect the radio button to the callback function
radio_systems.on_clicked(switch_system)

# Connect Event Handlers
fig.canvas.mpl_connect('scroll_event', on_scroll)
fig.canvas.mpl_connect('button_press_event', on_press)
fig.canvas.mpl_connect('motion_notify_event', on_motion)
fig.canvas.mpl_connect('button_release_event', on_release)
fig.canvas.mpl_connect('key_press_event', on_key)

# Running the animation

# Calculate the number of frames needed
num_frames = int((END_DAY - START_DAY) / DAYS_PER_FRAME) + 1

print(f"Starting animation: {num_frames} frames, {DAYS_PER_FRAME} day(s) per frame, interval={INTERVAL_MS}ms")
print("\n--- Interactive Controls ---")
print("Mouse Wheel: Zoom in/out")
print("Left Mouse Button + Drag: Pan")
print("Press 'space': Pause/Unpause animation")
print("Press 'r': Reset view")
print("--------------------------\n")
# Create the animation object
ani = animation.FuncAnimation(fig, animate, frames=num_frames,
                              init_func=init, blit=False, interval=INTERVAL_MS,
                              repeat=True) # repeat=False stops after one loop

# Display the animation
plt.show()

print("Animation finished or window closed.")