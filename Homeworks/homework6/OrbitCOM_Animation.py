#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 11:43:39 2025

@author: petershea
"""
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from OrbitCOM import readOrbit

# reads in positional data for each galaxy, can subtract to get orbitaal motion 
# relative to an object
x1,y1,z1 = readOrbit('Orbit_MW.txt')[1] #- readOrbit('Orbit_MW.txt')[1]
x2,y2,z2 = readOrbit('Orbit_M31.txt')[1] #- readOrbit('Orbit_MW.txt')[1]
x3,y3,z3 = readOrbit('Orbit_M33.txt')[1] #- readOrbit('Orbit_MW.txt')[1]

# Combine data into lists
x_data = [x1, x2, x3]
y_data = [y1, y2, y3]
z_data = [z1, z2, z3]

num_objects = len(x_data)
num_frames = min(len(x1), len(x2), len(x3))  # Ensure all have the same length

# Initialize figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set axis limits dynamically based on data range
ax.set_xlim([min(map(min, x_data)), max(map(max, x_data))])
ax.set_ylim([min(map(min, y_data)), max(map(max, y_data))])
ax.set_zlim([min(map(min, z_data)), max(map(max, z_data))])

# Sets axis Labels and Titles
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title("COM Orbital Paths")

# lables for each object
labels = ["MW", "M31", "M33"]

# Create line objects for each object
lines = [ax.plot([], [], [], label=labels[i])[0] for i in range(num_objects)]
points = [ax.scatter([], [], [], marker='o', s=50) for _ in range(num_objects)]  # 's' controls marker size

# Initialize function
def init():
    for line, point in zip(lines, points):
        line.set_data([], [])
        line.set_3d_properties([])
        point._offsets3d = ([], [], [])  # Initialize scatter points
    return lines + points

# Update function
def update(frame):
    for i, (line, point) in enumerate(zip(lines, points)):
        # Update trajectory line
        line.set_data(x_data[i][:frame], y_data[i][:frame])
        line.set_3d_properties(z_data[i][:frame])

        # Update scatter point (current position)
        point._offsets3d = ([x_data[i][frame]], [y_data[i][frame]], [z_data[i][frame]])

    return lines + points

# Animate
ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, interval=30, blit=False)

# Show animation
plt.legend()
plt.show()

ani.save('COM_Orbits.gif', writer='pillow', fps=30)

