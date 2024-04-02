import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def read_data(filename):
    data = np.loadtxt(filename, delimiter=',', skiprows=1)

    time = data[:, 0]
    u = data[:, 1]
    v = data[:, 2]

    return time, u, v


def plot_data(filename):
    time, u, v = read_data(filename)

    # Plot 2D graphs for u(t), v(t)
    plt.figure(figsize=(10, 5))
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    ax_u = plt.subplot(2, 1, 1)
    plt.plot(time, u, label='u(t)', color='r')
    plt.title('u(t)')
    plt.xlabel('t')
    plt.ylabel('u')
    plt.legend()

    ax_v = plt.subplot(2, 1, 2)
    plt.plot(time, v, label='v(t)', color='g')
    plt.title('v(t)')
    plt.xlabel('t')
    plt.ylabel('v')
    plt.legend()

    t_limits = ax_u.get_xlim()
    u_limits = ax_u.get_ylim()
    v_limits = ax_v.get_ylim()

    fig2d = plt.figure(figsize=(10, 5))
    ax2d = fig2d.add_subplot(111)
    ax2d.plot(time, u, label='u(t)', color='r')
    ax2d.plot(time, v, label='v(t)', color='g')
    ax2d.set_xlabel('t')
    ax2d.set_ylabel('u,v')
    ax2d.set_title('u(t), v(t)')
    ax2d.legend()

    # Animated graph
    fig_animated = plt.figure(figsize=(10, 5))
    plt.subplots_adjust(hspace=0.5, wspace=0.5)

    u_animated = fig_animated.add_subplot(211)
    u_animated.set_xlim(t_limits)
    u_animated.set_ylim(u_limits)
    v_animated = fig_animated.add_subplot(212)
    v_animated.set_xlim(t_limits)
    v_animated.set_ylim(v_limits)

    # Initialize first subplot with the first frame
    line_u, = u_animated.plot(time[0], u[0], label='u(t)', color='r')
    u_animated.set_xlabel('t')
    u_animated.set_ylabel('u')
    u_animated.set_title('u(t)')
    u_animated.legend()

    # Initialize second subplot with the first frame
    line_v, = v_animated.plot(time[0], v[0], label='v(t)', color='g')
    v_animated.set_xlabel('t')
    v_animated.set_ylabel('v')
    v_animated.set_title('v(t)')
    v_animated.legend()

    current_point_u, = u_animated.plot(time[0], u[0], 'ro')
    current_point_v, = v_animated.plot(time[0], v[0], 'go')

    def update(frame):
        nonlocal current_point_u, current_point_v
        current_point_u.remove()
        current_point_v.remove()
        line_u.set_data(time[:frame+1], u[:frame+1])
        line_v.set_data(time[:frame+1], v[:frame+1])
        current_point_u, = u_animated.plot(time[frame], u[frame], 'ro')
        current_point_v, = v_animated.plot(time[frame], v[frame], 'go')
        return line_u, current_point_u, line_v, current_point_v,

    ani = FuncAnimation(fig_animated, update, frames=len(time), interval=1, blit=True)
    # ani.save('animation.gif', writer=PillowWriter(fps=15))

    # Show plots
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <data_file>")
        sys.exit(1)

    file = sys.argv[1]
    plot_data(file)
