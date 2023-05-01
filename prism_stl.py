#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Charlotte Bodenmüller

Description: 
    This program calculates the vertices of a customized prism.
Input: Parameters that define a prism.
Output: STL file of the prism

"""

import numpy as np
import math as m

# Plot libs
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Reading input
print("Welcome to the prism calculator, please enter the following parameters\n")
faces = int(input("Number of side faces (>2)\n"))
facewidth = int(input("Length of side faces\n"))
height = int(input("Prism height\n"))
location_x = int(input("x-axis of the center location\n"))
location_y = int(input("y-axis of the center location\n"))
location_z = int(input("z-axis of the center location\n"))
angle_x = int(input("Rotation angle for the x-axis in degrees\n"))
angle_y = int(input("Rotation angle for the y-axis in degrees\n"))
angle_z = int(input("Rotation angle for the z-axis in degrees\n"))



# Calculates vertices of prism
def prism(height, facewidth, n):
    coordinates_vertices = []
    
    if n < 3:
        print("Error: Number of side faces should be >= 3")
        return coordinates_vertices
    
    # Calc lengths
    angle = 360/n
    center_to_vertice = (facewidth/2)/np.sin(np.deg2rad(angle/2))
    center_to_middle_of_face = center_to_vertice*np.cos(np.deg2rad(angle/2))
    
    ## Calc vertices of basement
    # Calc first two vertices
    coordinates_vertices.append(np.array([center_to_middle_of_face, -facewidth/2, 0]))
    coordinates_vertices.append(np.array([center_to_middle_of_face, facewidth/2, 0]))
    
    # Calc residual vertices
    vector = coordinates_vertices[1] - coordinates_vertices[0]
    rotation_matrix = rot_z(angle)
    
    for i in range(2, n): 
        vector = np.dot(rotation_matrix, vector)
        coordinates_vertices.append(coordinates_vertices[i-1] + vector)
        
    ## Calc vertices of top face
    for i in range(n):
        coordinates_vertices.append(coordinates_vertices[i] + np.array([0,0,height]))
    
    return coordinates_vertices


def rotating_prism(coordinates_vertices, angle_x, angle_y, angle_z):
    for i in range(len(coordinates_vertices)):
        coordinates_vertices[i] = np.dot(rot_x(angle_x), coordinates_vertices[i])
        coordinates_vertices[i] = np.dot(rot_y(angle_y), coordinates_vertices[i])
        coordinates_vertices[i] = np.dot(rot_z(angle_z), coordinates_vertices[i])  
        
    return coordinates_vertices
        
def rot_x(theta):
    theta = np.deg2rad(theta)
    return np.array([[ 1, 0           , 0           ],
                      [ 0, m.cos(theta),-m.sin(theta)],
                      [ 0, m.sin(theta), m.cos(theta)]])
  
def rot_y(theta):
    theta = np.deg2rad(theta)
    return np.array([[ m.cos(theta), 0, m.sin(theta)],
                      [ 0           , 1, 0           ],
                      [-m.sin(theta), 0, m.cos(theta)]])
  
def rot_z(theta):
    theta = np.deg2rad(theta)
    return np.array([[ m.cos(theta), -m.sin(theta), 0 ],
                      [ m.sin(theta), m.cos(theta) , 0 ],
                      [ 0           , 0            , 1 ]])

# Plots vertices of a prism
def plot_prism_vertices(coordinates_vertices):   
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cube = coordinates_vertices
    hull = ConvexHull(cube)
    # draw the polygons of the convex hull
    for s in hull.simplices:
        tri = Poly3DCollection([cube[s]])
        tri.set_color("blue")
        tri.set_alpha(0.5)
        ax.add_collection3d(tri)
    # draw the vertices
    ax.scatter(cube[:, 0], cube[:, 1], cube[:, 2], marker='o', color='red')
    plt.show()


coordinates_vertices = prism(height, facewidth, faces)
coordinates_vertices = np.array(coordinates_vertices)
coordinates_vertices = rotating_prism(coordinates_vertices, angle_x, angle_y, angle_z)

# Change location
for i in range(len(coordinates_vertices)):
    coordinates_vertices[i][0] = coordinates_vertices[i][0] + location_x
    coordinates_vertices[i][1] = coordinates_vertices[i][1] + location_y
    coordinates_vertices[i][2] = coordinates_vertices[i][2] + location_z
    
plot_prism_vertices(coordinates_vertices)


# Saving to file

coordinates_vertices_basement = coordinates_vertices[:int(len(coordinates_vertices)/2)]
coordinates_vertices_top = coordinates_vertices[int(len(coordinates_vertices)/2):]

def vertex_to_str(arr):
    return str(arr[0]) + " " + str(arr[1]) + " " + str(arr[2]) 

def triangle_to_stl (v1, v2, v3,normal_vector):
    f.write("facet normal " + vertex_to_str(normal_vector) + "\n")
    f.write("   outer loop\n")
    f.write("       vertex " + vertex_to_str(v1) + "\n")
    f.write("       vertex " + vertex_to_str(v2) + "\n")
    f.write("       vertex " + vertex_to_str(v3) + "\n")
    f.write("   end loop\n")
    f.write("endfacet\n")

with open('prism_stl.txt', 'w') as f:
    f.write("solid prism\n")
    
    ## Top and basement triangles
    # Calc norm vector 
    a = coordinates_vertices_top[1] - coordinates_vertices_top[0]
    b = coordinates_vertices_top[2] - coordinates_vertices_top[0]
    normal_vector_top = np.cross(a,b)
    
    # Triangles of top face 
    for i in range(1, len(coordinates_vertices_top)-1):
        triangle_to_stl(coordinates_vertices_top[0],coordinates_vertices_top[i],coordinates_vertices_top[i+1],normal_vector_top)

    # Triangles of basement face 
    for i in range(1, len(coordinates_vertices_basement)-1):       
        triangle_to_stl(coordinates_vertices_basement[0],coordinates_vertices_basement[i],coordinates_vertices_basement[i+1],normal_vector_top*(-1))
        

        
    ## Triangles of side faces
    for i in range(len(coordinates_vertices_top)):
        
        # Calc index of third vertex
        j = -1
        if i == len(coordinates_vertices_top) -1:
            j = 0
        else:
            j = i +1
            
        # Calc norm vector
        a = coordinates_vertices_top[i] - coordinates_vertices_basement[i]
        b = coordinates_vertices_top[i] - coordinates_vertices_top[j]
        normal_vector_side = np.cross(a,b)
        
        # Side face 1
        triangle_to_stl(coordinates_vertices_basement[i],coordinates_vertices_top[i],coordinates_vertices_top[j],normal_vector_side)
        
        # Side face 2
        triangle_to_stl(coordinates_vertices_basement[i],coordinates_vertices_basement[j],coordinates_vertices_top[j],normal_vector_side)

    f.write("endsolid prism")
    print("Prism saved in stl file")





