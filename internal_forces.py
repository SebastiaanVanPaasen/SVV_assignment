import numpy as np
import matplotlib.pyplot as plt
import copy
#from SVV_assignment.SVV_assignment.geometry import *
from geometry import *

# aileron parameters
l_a = 1.691
h_a = 0.173

# force locations along the beam
x_begin = 0
x_1 = 0.149
x_2 = 0.554
x_3 = 1.541
x_end = 1.691

# actuator distances
theta = 28  # degrees
x_a = 0.272
x_a1 = x_2 - x_a / 2
y_a1 = np.cos(np.radians(h_a / 2 - np.tan(theta) * h_a / 2))
z_a1 = np.cos(theta) * h_a / 2 + np.sin(theta) * h_a / 2
x_a2 = x_2 + x_a / 2

# load parameters
q = 2710
p = 37900
E = 73.1 * (10 ** 9)
G = 28 * (10 ** 9)
delta_1 = 0.00681
delta_3 = 0.0203
izz, iyy = get_moi(20)


class Force:
        
    def __init__(self, magnitude, direction, position):
        self.position = position
        self.direction = direction
        self.magnitude = magnitude
        self.d = np.linalg.norm(self.direction)
        self.x = self.magnitude * self.direction[0]/self.d
        self.y = self.magnitude * self.direction[1]/self.d
        self.z = self.magnitude * self.direction[2]/self.d
        
    def modify(self, magnitude, direction = 1):
        self.magnitude = magnitude
        self.direction *= direction
        self.x = self.magnitude * self.direction[0]/self.d
        self.y = self.magnitude * self.direction[1]/self.d
        self.z = self.magnitude * self.direction[2]/self.d
        
    def resultant(self, force_2):
        position = self.position
        direction= np.zeros(3)
        x = self.x + force_2.x 
        y = self.y + force_2.y
        z = self.z + force_2.z
        vector = [x,y,z]
        magnitude = np.linalg.norm(vector)
        for i in range(3):
            direction[i] = vector[i]/magnitude
        res_force = Force(magnitude, direction, position)
        return res_force
        
    def determine_force(self, direction):
        if direction == 'y':
            return self.y
        elif direction == 'z':
            return self.z
        elif direction == 'x':
            return self.x
        else:
            return 0
    
    def determine_moment(self, position):
        distance = self.position - position
#        moments = np.zeros(3)
#        moments[0] = (self.y * distance[2] - self.z * distance[1])
#        moments[1] = (self.x * distance[2] - self.z * distance[0])
#        moments[2] = (-self.x * distance[1] + self.y * distance[0])
        moments = np.cross(distance, [self.x,self.y,self.z])
        
        return moments

    # rotation specifically in yz axis
    def rotate(self, angle):
        rotated_force = copy.deepcopy(self)
        
        rotated_force.y = self.y * np.cos(np.radians(angle)) - self.z * np.sin(np.radians(angle))
        rotated_force.z = self.z * np.cos(np.radians(angle)) + self.y * np.sin(np.radians(angle))
        rotated_force.direction[1] = rotated_force.y/rotated_force.magnitude
        rotated_force.direction[2] = rotated_force.z/rotated_force.magnitude
        
        return rotated_force


# def calc_reaction_forces():  THIS ONE
y_1 = Force(1, np.array([0, 1, 0]), np.array([x_1, 0, 0]))
y_2 = Force(1, np.array([0, 1, 0]), np.array([x_2, 0, 0]))
y_3 = Force(1, np.array([0, 1, 0]), np.array([x_3, 0, 0]))
z_1 = Force(1, np.array([0, 0, 1]), np.array([x_1, 0, 0]))
z_2 = Force(1, np.array([0, 0, 1]), np.array([x_2, 0, 0]))
z_3 = Force(1, np.array([0, 0, 1]), np.array([x_3, 0, 0]))
a_1y = Force(1, np.array([0, 1, 0]), np.array([x_a1, y_a1, z_a1]))
a_1z = Force(1, np.array([0, 0, 1]), np.array([x_a1, y_a1, z_a1]))
# TODO: decompose the forces below into the corresponding tilted axis
p_y = Force(abs(p * np.sin(theta)), np.array([0, -1, 0]), np.array([x_a2, y_a1, z_a1]))
p_z = Force(abs(p * np.cos(theta)), np.array([0, 0, -1]), np.array([x_a2, y_a1, z_a1]))
q_y = Force(abs(q * np.cos(theta) * l_a), np.array([0, -1, 0]), np.array([l_a / 2, 0, 0]))
q_z = Force(abs(q * np.sin(theta) * l_a), np.array([0, 0, 1]), np.array([l_a / 2, 0, 0]))

forces = [y_1, y_2, y_3, z_1, z_2, z_3, a_1y, a_1z, p_y, p_z, q_y, q_z]
        
#R1 = Force(np.sqrt(2), np.array([0,1,1]), np.array([x_1,0,0]))
#R2 = Force(np.sqrt(2), np.array([0,1,1]), np.array([x_2,0,0]))
#R3 = Force(np.sqrt(2), np.array([0,1,1]), np.array([x_3,0,0]))
#A1 = Force(np.sqrt(2), np.array([0,1,1]), np.array([x_a1,y_a1,z_a1]))
#P = Force(p, np.array([0,-1*np.sin(theta),-1*np.cos(theta)]), np.array([x_a2,y_a1,z_a1]))
#Q = Force(q, np.array([0,-1*np.cos(theta),1*np.sin(theta)]), np.array([l_a/2, 0, 0]))

#forces = [R1, R2, R3, A1, P, Q]

sum_forces_y = []
sum_forces_z = []
sum_moments_x = []
sum_moments_y = []
sum_moments_z = []

for force in forces:
    sum_forces_y.append(force.determine_force('y'))
    sum_forces_z.append(force.determine_force('z'))
    sum_moments_x.append(force.determine_moment([0, 0, 0])[0])
    sum_moments_y.append(force.determine_moment([0, 0, 0])[1])
    sum_moments_z.append(force.determine_moment([0, 0, 0])[2])

# print(sum_forces_y)
# print(sum_forces_z)
# print(sum_moments)

ce_eq_z_1 = (1 / (E * izz)) * np.array([(l_a - x_1) ** 3 / 3,
                                        (l_a - x_2) ** 3 / 3 - (x_1 + x_2) * (l_a - x_2) ** 2 / 2 + x_1 * x_2 * (
                                                l_a - x_2),
                                        (l_a - x_3) ** 3 / 3 - (x_1 + x_3) * (l_a - x_3) ** 2 / 2 + x_1 * x_3 * (
                                                l_a - x_3),
                                        0., 0., 0.,
                                        (l_a - x_a1) ** 3 / 3 - (x_1 + x_a1) * (l_a - x_a1) ** 2 / 2 + x_1 * x_a1 * (
                                                l_a - x_a1),
                                        0.,
                                        (l_a - x_a2) ** 3 / 3 - (x_1 + x_a2) * (l_a - x_a2) ** 2 / 2 + x_1 * x_a2 * (
                                                l_a - x_a2),
                                        0.,
                                        -0.5 * (l_a ** 3 / 3 - x_1 * l_a * l_a / 2), delta_1])

ce_eq_z_3 = (1 / (E * izz)) * np.array(
    [(l_a - x_1) ** 3 / 3 - (x_3 + x_1) * (l_a - x_1) ** 2 / 2 + x_3 * x_1 * (l_a - x_1),
     (l_a - x_2) ** 3 / 3 - (x_3 + x_2) * (l_a - x_2) ** 2 / 2 + x_3 * x_2 * (l_a - x_2),
     (l_a - x_3) ** 3 / 3,
     0., 0., 0.,
     (l_a - x_a1) ** 3 / 3 - (x_3 + x_a1) * (l_a - x_a1) ** 2 / 2 + x_3 * x_a1 * (l_a - x_a1),
     0.,
     (l_a - x_a2) ** 3 / 3 - (x_3 + x_a2) * (l_a - x_a2) ** 2 / 2 + x_3 * x_a2 * (l_a - x_a2),
     0.,
     -0.5 * (l_a ** 3 / 3 - x_3 * l_a * l_a / 2), delta_3])

ce_eq_y = (1 / (E * izz)) * np.array([0., 0., 0., (l_a - x_1) ** 3 / 3,
                                      (l_a - x_2) ** 3 / 3 - (x_1 + x_2) * (l_a - x_2) ** 2 / 2 + x_1 * x_2 * (
                                              l_a - x_2),
                                      (l_a - x_3) ** 3 / 3 - (x_1 + x_3) * (l_a - x_3) ** 2 / 2 + x_1 * x_3 * (
                                              l_a - x_3),
                                      0.,
                                      (l_a - x_a1) ** 3 / 3 - (x_1 + x_a1) * (l_a - x_a1) ** 2 / 2 + x_1 * x_a1 * (
                                              l_a - x_a1),
                                      0.,
                                      (l_a - x_a2) ** 3 / 3 - (x_1 + x_a2) * (l_a - x_a2) ** 2 / 2 + x_1 * x_a2 * (
                                              l_a - x_a2),
                                      0.,
                                      -0.5 * (l_a ** 3 / 3 - x_1 * l_a * l_a / 2), 0.])

system = [sum_forces_y, sum_forces_z, sum_moments_x, sum_moments_y, sum_moments_z, ce_eq_y, ce_eq_z_1, ce_eq_z_3]

# solving reaction forces
sys_mat = np.zeros((8, 8))
sys_vec = np.zeros(8)
for i in range(len(system)):
    for j in range(len(system[0]) - 4):
        sys_mat[i][j] = system[i][j]
    sys_vec[i] = -1 * np.sum(system[i][len(system):])

unk = np.linalg.solve(sys_mat, sys_vec)
for i in range(len(forces) - 4):
    forces[i].modify(abs(unk[i]),int(np.sign(unk[i])))  # TODO: make magnitude positive and the direction according to the correct sign
    

# return forces  ENDS HERE

class Slice:
        
    def __init__(self, position, dx):
        self.position = position
        self.dx = dx
        self.vx = 0.  # Force(1., np.array([1,0,0]), np.array([x, 0, 0]))
        self.vy = 0.  # Force(1., np.array([0,1,0]), np.array([x, 0, 0]))
        self.vz = 0.  # Force(1., np.array([0,1,0]), np.array([x, 0, 0]))
        self.mx = 0.
        self.my = 0.
        self.mz = 0.
        self.dy = 0.

#    def int_dist(self, applied_forces, l_a):
#        ext_forces = []
#        app_forces = copy.deepcopy(applied_forces)
#        for i in [-2, -1]:
#            app_forces[i].magnitude *= self.position[0] * 1 / l_a
#            app_forces[i].position[0] *= self.position[0] * 1 / l_a
#
#        for i in range(len(app_forces)):
#            if app_forces[i].position[0] < self.position[0]:
#                ext_forces.append(app_forces[i])
#        for i in range(len(ext_forces)):
#            self.vx += -1 * ext_forces[i].determine_force('x')
#            self.vy += -1 * ext_forces[i].determine_force('y')
#            self.vz += -1 * ext_forces[i].determine_force('z')
#            self.mx += -1 * ext_forces[i].determine_moment(self.position)[0]
#            self.my += -1 * ext_forces[i].determine_moment(self.position)[1]
#            self.mz += -1 * ext_forces[i].determine_moment(self.position)[2]
#        return app_forces, ext_forces  # Did this to just check that the acquired forces are correct. Saving this data is unnecessary.
        
    def int_dist(self, applied_forces, l_a):
        ext_forces = []
        app_forces = copy.deepcopy(applied_forces)
        for i in [-2, -1]:
            app_forces[i].modify(app_forces[i].magnitude*self.position[0] * 1 / l_a)
            app_forces[i].position[0] *= self.position[0] * 1 / l_a

        for i in range(len(app_forces)):
            if app_forces[i].position[0] < self.position[0]:
                ext_forces.append(app_forces[i])
        for i in range(len(ext_forces)):
            self.vx += -1 * ext_forces[i].determine_force('x')
            self.vy += -1 * ext_forces[i].determine_force('y')
            self.vz += -1 * ext_forces[i].determine_force('z')
            self.mx += -1 * ext_forces[i].determine_moment(self.position)[0]
            self.my += -1 * ext_forces[i].determine_moment(self.position)[1]
            self.mz += -1 * ext_forces[i].determine_moment(self.position)[2]
        return app_forces, ext_forces  # Did this to just check that the acquired forces are correct. Saving this data is unnecessary.

    # def macauley(self):


# def distribution(forces, bc1, bc2, l_a, dx):
# x_bound = []
# for i in range(len(forces)):
#    if forces[i].position[0] != l_a/2:
#        x_bound.append(forces[i].position[0])
# x_bound = (list(set(x_bound)))
# x_bound.sort()

dx = 0.0001
x_slice = np.arange(0., l_a + dx, dx)
total_n = len(x_slice)
v_y = np.zeros(total_n)
v_z = np.zeros(total_n)
m_z = np.zeros(total_n)
m_y = np.zeros(total_n)
d_y = np.zeros(total_n)

slice_list = []
slice_app_force = []  # unnecessary
for i in range(total_n):
    slice_list.append(Slice([x_slice[i],0,0], dx))
    slice_app_force.append(slice_list[i].int_dist(forces, l_a))
    v_y[i] = slice_list[i].vy
    v_z[i] = slice_list[i].vz
    m_z[i] = slice_list[i].mz
    m_y[i] = slice_list[i].my
for i in range(total_n-1):
    d_y[i] = np.trapz([m_z[i],m_z[i+1]],[x_slice[i],x_slice[i+1]])


f1 = plt.figure()
plt.plot(x_slice, m_z)
f2 = plt.figure()
plt.plot(x_slice, m_y)
plt.show()