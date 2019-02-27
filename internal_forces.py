import numpy as np
import matplotlib.pyplot as plt
import copy

from SVV_assignment.SVV_assignment.shear import get_shear_flow
from SVV_assignment.SVV_assignment.torque import get_torque
from SVV_assignment.SVV_assignment.shear_center import *

# from shear_center import *
# from torque import get_torque
# from shear import get_shear_flow

# aileron parameters
l_a = 1.691
h_a = 0.173
c_a = 0.484
t_skin = 1.1 / 1000
t_spar = 2.5 / 1000

# force locations along the beam
x_begin = 0
x_1 = 0.149
x_2 = 0.554
x_3 = 1.541
x_end = 1.691

# actuator distances
theta = 26  # degrees
x_a = 0.272
x_a1 = x_2 - x_a / 2
y_a1 = h_a / 2  # np.cos(np.radians(h_a / 2 - np.tan(np.radians(theta)) * h_a / 2))
z_a1 = h_a / 2  # np.cos(np.radians(theta)) * h_a / 2 + np.sin(np.radians(theta)) * h_a / 2
x_a2 = x_2 + x_a / 2

# load parameters
q = 2710
qy = -q * np.cos(np.radians(theta))
qz = q * np.sin(np.radians(theta))
p = 37900
E = 73.1 * 10 ** 9
G = 28 * 10 ** 9
delta_1 = 0.00681
delta_1y = delta_1 * np.cos(np.radians(theta))
delta_1z = -delta_1 * np.sin(np.radians(theta))
delta_3 = 0.0203
delta_3y = delta_3 * np.cos(np.radians(theta))
delta_3z = -delta_3 * np.sin(np.radians(theta))

izz, iyy = 5.697e-6, 6.947e-5  # get_moi(20)


class Force:

    def __init__(self, magnitude, direction, position):
        self.position = position
        self.direction = direction
        self.magnitude = magnitude
        self.d = np.linalg.norm(self.direction)
        self.x = self.magnitude * self.direction[0] / self.d
        self.y = self.magnitude * self.direction[1] / self.d
        self.z = self.magnitude * self.direction[2] / self.d

    def modify(self, magnitude, direction=1):
        self.magnitude = magnitude
        self.direction *= direction
        self.x = self.magnitude * self.direction[0] / self.d
        self.y = self.magnitude * self.direction[1] / self.d
        self.z = self.magnitude * self.direction[2] / self.d

    def resultant(self, force_2):
        position = self.position
        direction = np.zeros(3)
        x = self.x + force_2.x
        y = self.y + force_2.y
        z = self.z + force_2.z
        vector = [x, y, z]
        magnitude = np.linalg.norm(vector)
        for i in range(3):
            direction[i] = vector[i] / magnitude
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
        moments = np.cross(distance, [self.x, self.y, self.z])

        return moments

    # rotation specifically in yz axis
    def rotate(self, angle):
        rotated_force = copy.deepcopy(self)

        rotated_force.y = self.y * np.cos(np.radians(angle)) - self.z * np.sin(np.radians(angle))
        rotated_force.z = self.z * np.cos(np.radians(angle)) + self.y * np.sin(np.radians(angle))
        rotated_force.direction[1] = rotated_force.y / rotated_force.magnitude
        rotated_force.direction[2] = rotated_force.z / rotated_force.magnitude

        return rotated_force

    def macauley(self, deflection, direction, x):
        term = 1 / 6 * np.array([self.x, self.y, self.z]) * direction * ((x - self.position[0]) ** 3) * np.heaviside(
            x - self.position[0], 1)
        if direction[2] != 0:
            term *= -1
        return term


# def calc_reaction_forces():  THIS ONE
y_1 = Force(1, np.array([0, 1, 0]), np.array([x_1, 0, 0]))
# y_1 = Force(1, np.array([0, 1, 0]), np.array([x_1, delta_1y, delta_1z]))
y_2 = Force(1, np.array([0, 1, 0]), np.array([x_2, 0, 0]))
y_3 = Force(1, np.array([0, 1, 0]), np.array([x_3, 0, 0]))
# y_3 = Force(1, np.array([0, 1, 0]), np.array([x_3, delta_3y, delta_3z]))

z_1 = Force(1, np.array([0, 0, 1]), np.array([x_1, 0, 0]))
# z_1 = Force(1, np.array([0, 0, 1]), np.array([x_1, delta_1y, delta_1z]))
z_2 = Force(1, np.array([0, 0, 1]), np.array([x_2, 0, 0]))
z_3 = Force(1, np.array([0, 0, 1]), np.array([x_3, 0, 0]))
# z_3 = Force(1, np.array([0, 0, 1]), np.array([x_3, delta_3y, delta_3z]))

a_1y = Force(1, np.array([0, 1, 0]), np.array([x_a1, y_a1, z_a1]))
a_1z = Force(1, np.array([0, 0, 1]), np.array([x_a1, y_a1, z_a1]))

p_y = Force(abs(p * np.sin(np.radians(theta))), np.array([0, -1, 0]), np.array([x_a2, y_a1, z_a1]))
p_z = Force(abs(p * np.cos(np.radians(theta))), np.array([0, 0, -1]), np.array([x_a2, y_a1, z_a1]))
q_y = Force(abs(q * np.cos(np.radians(theta)) * l_a), np.array([0, -1, 0]),
            np.array([l_a / 2, 0, h_a / 2 - 0.25 * c_a]))
q_z = Force(abs(q * np.sin(np.radians(theta)) * l_a), np.array([0, 0, 1]), np.array([l_a / 2, 0, h_a / 2 - 0.25 * c_a]))

forces = [y_1, y_2, y_3, z_1, z_2, z_3, a_1y, a_1z, p_y, p_z, q_y, q_z]

sum_forces_y = []
sum_forces_z = []
sum_moments_x = []
sum_moments_y = []
sum_moments_z = []
# defl_y1 = []
# defl_y2 = []
# defl_y3 = []
# defl_z1 = []
# defl_z2 = []
# defl_z3 = []

for force in forces:
    sum_forces_y.append(force.determine_force('y'))
    sum_forces_z.append(force.determine_force('z'))
    sum_moments_x.append(force.determine_moment([0, 0, 0])[0])
    sum_moments_y.append(force.determine_moment([0, 0, 0])[1])
    sum_moments_z.append(force.determine_moment([0, 0, 0])[2])

# for i in range(len(forces)):
#    sum_forces_y.append(forces[i].determine_force('y'))
#    sum_forces_z.append(forces[i].determine_force('z'))
#    sum_moments_x.append(forces[i].determine_moment([0, 0, 0])[0])
#    sum_moments_y.append(forces[i].determine_moment([0, 0, 0])[1])
#    sum_moments_z.append(forces[i].determine_moment([0, 0, 0])[2])
#    if i < len(forces)-2:
#        defl_y1.append(forces[i].macauley(delta_1y, [0,1,0], x_1))
#        defl_y2.append(forces[i].macauley(0., [0,1,0], x_2))
#        defl_y3.append(forces[i].macauley(delta_3y, [0,1,0], x_3))
#        defl_z1.append(forces[i].macauley(delta_1y, [0,0,1], x_1))
#        defl_z2.append(forces[i].macauley(0, [0,0,1], x_2))
#        defl_z3.append(forces[i].macauley(delta_3y, [0,0,1], x_3))     #wanted to do this but how do you insert elements from a list to another list


defl_y1 = 1 / 6 * np.array([(x_1 - x_1) ** 3 * np.heaviside(x_1 - x_1, 1),
                            (x_1 - x_2) ** 3 * np.heaviside(x_1 - x_2, 1),
                            (x_1 - x_3) ** 3 * np.heaviside(x_1 - x_3, 1),
                            0., 0., 0.,
                            (x_1 - x_a1) ** 3 * np.heaviside(x_1 - x_a1, 1),
                            0.,
                            6 * x_1, 6 * 1., 0., 0.,
                            forces[8].y * ((x_1 - x_a2) ** 3) * np.heaviside(x_1 - x_a2, 1),
                            0.,
                            6 * qy / 24 * (x_1 ** 4),
                            0.,
                            -E * izz * 6 * delta_1y])

defl_y2 = 1 / 6 * np.array([(x_2 - x_1) ** 3 * np.heaviside(x_2 - x_1, 1),
                            (x_2 - x_2) ** 3 * np.heaviside(x_2 - x_2, 1),
                            (x_2 - x_3) ** 3 * np.heaviside(x_2 - x_3, 1),
                            0., 0., 0.,
                            (x_2 - x_a1) ** 3 * np.heaviside(x_2 - x_a1, 1),
                            0.,
                            6 * x_2, 6 * 1., 0., 0.,
                            forces[8].y * ((x_2 - x_a2) ** 3) * np.heaviside(x_2 - x_a2, 1),
                            0.,
                            6 * qy / 24 * (x_2 ** 4),
                            0., 0.])

defl_y3 = 1 / 6 * np.array([(x_3 - x_1) ** 3 * np.heaviside(x_3 - x_1, 1),
                            (x_3 - x_2) ** 3 * np.heaviside(x_3 - x_2, 1),
                            (x_3 - x_3) ** 3 * np.heaviside(x_3 - x_3, 1),
                            0., 0., 0.,
                            (x_3 - x_a1) ** 3 * np.heaviside(x_3 - x_a1, 1),
                            0.,
                            6 * x_3, 6 * 1., 0., 0.,
                            forces[8].y * ((x_3 - x_a2) ** 3) * np.heaviside(x_3 - x_a2, 1),
                            0.,
                            6 * qy / 24 * (x_3 ** 4),
                            0.,
                            -E * izz * 6 * delta_3y])

defl_z1 = 1 / 6 * np.array([0., 0., 0.,
                            -(x_1 - x_1) ** 3 * np.heaviside(x_1 - x_1, 1),
                            -(x_1 - x_2) ** 3 * np.heaviside(x_1 - x_2, 1),
                            -(x_1 - x_3) ** 3 * np.heaviside(x_1 - x_3, 1),
                            0.,
                            -(x_1 - x_a1) ** 3 * np.heaviside(x_1 - x_a1, 1),
                            0., 0., 6 * x_1, 6 * 1.,
                            0.,
                            -forces[9].z * ((x_1 - x_a2) ** 3) * np.heaviside(x_1 - x_a2, 1),
                            0.,
                            -6 * qz / 24 * (x_1 ** 4),
                            -E * iyy * 6 * delta_1z])

defl_z2 = 1 / 6 * np.array([0., 0., 0.,
                            -(x_2 - x_1) ** 3 * np.heaviside(x_2 - x_1, 1),
                            -(x_2 - x_2) ** 3 * np.heaviside(x_2 - x_2, 1),
                            -(x_2 - x_3) ** 3 * np.heaviside(x_2 - x_3, 1),
                            0.,
                            -(x_2 - x_a1) ** 3 * np.heaviside(x_2 - x_a1, 1),
                            0., 0., 6 * x_2, 6 * 1.,
                            0.,
                            -forces[9].z * ((x_2 - x_a2) ** 3) * np.heaviside(x_2 - x_a2, 1),
                            0.,
                            -6 * qz / 24 * (x_2 ** 4),
                            0.])

defl_z3 = 1 / 6 * np.array([0., 0., 0.,
                            -(x_3 - x_1) ** 3 * np.heaviside(x_3 - x_1, 1),
                            -(x_3 - x_2) ** 3 * np.heaviside(x_3 - x_2, 1),
                            -(x_3 - x_3) ** 3 * np.heaviside(x_3 - x_3, 1),
                            0.,
                            -(x_3 - x_a1) ** 3 * np.heaviside(x_3 - x_a1, 1),
                            0., 0., 6 * x_3, 6 * 1.,
                            0.,
                            -forces[9].z * ((x_3 - x_a2) ** 3) * np.heaviside(x_3 - x_a2, 1),
                            0.,
                            -6 * qz / 24 * (x_3 ** 4),
                            -E * iyy * 6 * delta_3z])

trig = np.zeros(16)
trig[6] = np.cos(np.radians(theta))
trig[7] = -np.sin(np.radians(theta))

system = [sum_forces_y, sum_forces_z, sum_moments_x, sum_moments_y, sum_moments_z, defl_z1, defl_z2, defl_z3,
          defl_y1, defl_y2, defl_y3, trig]

for i in range(5):
    for j in range(4):
        system[i].insert((j + 8), 0)

# solving reaction forces
sys_mat = np.zeros((12, 12))
sys_vec = np.zeros(12)
for i in range(len(system)):
    for j in range(len(system[0]) - 4):
        sys_mat[i][j] = system[i][j]
    sys_vec[i] = -1 * np.sum(system[i][len(system):])

unk = np.linalg.solve(sys_mat, sys_vec)
for i in range(len(forces) - 4):
    forces[i].modify(abs(unk[i]), int(np.sign(unk[i])))

int_constant_z = [unk[10], unk[11]]
int_constant_y = [unk[8], unk[9]]

# return forces  ENDS HERE
forces_resultant = []
forces_global = []
forces_resultant_global = []

for i in range(len(forces)):
    forces_global.append(forces[i].rotate(28))
for i in range(len(forces)):
    if i % 2 == 0:
        forces_resultant.append(forces[i].resultant(forces[i + 1]))
        forces_resultant_global.append(forces_global[i].resultant(forces_global[i + 1]))


class Slice:

    def __init__(self, position, dx):
        self.position = position
        self.dx = dx
        self.vx = 0.
        self.vy = 0.
        self.vz = 0.
        self.mx = 0.
        self.my = 0.
        self.mz = 0.
        self.dy = 0.

    def int_dist(self, applied_forces, l_a):
        ext_forces = []
        app_forces = copy.deepcopy(applied_forces)
        for i in [-2, -1]:
            app_forces[i].modify(app_forces[i].magnitude * self.position[0] * 1 / l_a)
            app_forces[i].position[0] *= self.position[0] * 1 / l_a

        for i in range(len(app_forces)):
            if app_forces[i].position[0] < self.position[0]:
                ext_forces.append(app_forces[i])
        for i in range(len(ext_forces)):
            self.vx += 1 * ext_forces[i].determine_force('x')
            self.vy += 1 * ext_forces[i].determine_force('y')
            self.vz += 1 * ext_forces[i].determine_force('z')
            self.mx += 1 * ext_forces[i].determine_moment(self.position)[0]
            self.my += 1 * ext_forces[i].determine_moment(self.position)[1]
            self.mz += 1 * ext_forces[i].determine_moment(self.position)[2]
        return app_forces, ext_forces  # Did this to just check that the acquired forces are correct. Saving this data is unnecessary.


# def distribution(forces, bc1, bc2, l_a, dx):
d_x = 0.0001
x_slice = np.arange(0., l_a + d_x, d_x)
total_n = len(x_slice)
v_y = np.zeros(total_n)
v_z = np.zeros(total_n)
m_z = np.zeros(total_n)
m_y = np.zeros(total_n)
m_x = np.zeros(total_n)

# ----------------------------------------------------------------------------------------------------------------------
n_booms = 260
cg_z, boom_locations = get_cg(n_booms)
moi_zz, moi_yy = get_moi(n_booms)
boom_pos, distances, boom_areas, top, bottom = get_boom_information(n_booms)
sc_z = get_shear_center()

print("sc : " + str(sc_z))
print("izz : " + str(moi_zz))
print("iyy : " + str(moi_yy))
# ----------------------------------------------------------------------------------------------------------------------

slice_list = []
slice_app_force = []  # unnecessary
twist = []
shear_flows = []
shear_stress_yz = []
shear_stress_xy = []
for i in range(total_n):
    slice_list.append(Slice([x_slice[i], 0, sc_z], d_x))
    slice_app_force.append(slice_list[i].int_dist(forces, l_a))
    v_y[i] = slice_list[i].vy
    v_z[i] = slice_list[i].vz
    m_z[i] = slice_list[i].mz
    m_y[i] = slice_list[i].my
    m_x[i] = slice_list[i].mx

    flow_i_shear, flow_ii_shear, twist_shear = get_shear_flow(n_booms, v_y[i], v_z[i], [x_slice[i], 0, sc_z])
    flow_i_torque, flow_ii_torque, twist_torque = get_torque(n_booms, m_x[i])

    twist.append(twist_shear + twist_torque)
    shear_flows.append([np.add(flow_i_shear, flow_i_torque), np.add(flow_ii_shear, flow_ii_torque)])
    shear_stress_yz.append([np.divide(np.add(flow_i_shear, flow_i_torque), t_skin),
                         np.divide(np.add(flow_ii_shear, flow_ii_torque), t_skin)])

    shear_stress_yz[i][0][-1] = shear_stress_yz[i][0][-1] * t_skin / t_spar
    shear_stress_yz[i][1][-1] = shear_stress_yz[i][1][-1] * t_skin / t_spar


# print("the twist : " + str(twist))
# print("shear flows : " + str(shear_flows))
# print("shear stress : " + str(shear_stress))


# deflection in y

def eq_def_y(x, constant):
    return 1 / (E * izz) * (1 / 6 * (forces[0].y * ((x - x_1) ** 3) * np.heaviside(x - x_1, 1) +
                                      forces[1].y * ((x - x_2) ** 3) * np.heaviside(x - x_2, 1) +
                                      forces[2].y * ((x - x_3) ** 3) * np.heaviside(x - x_3, 1) +
                                      forces[6].y * ((x - x_a1) ** 3) * np.heaviside(x - x_a1, 1) +
                                      forces[8].y * ((x - x_a2) ** 3) * np.heaviside(x - x_a2, 1)) + qy / 24 * (
                                     x ** 4) + constant[0] * x + constant[1])


# deflection in z


def eq_def_z(x, constant):
    return 1 / (E * iyy) * (-1 / 6 * (forces[3].z * ((x - x_1) ** 3) * np.heaviside(x - x_1, 1) +
                                       forces[4].z * ((x - x_2) ** 3) * np.heaviside(x - x_2, 1) +
                                       forces[5].z * ((x - x_3) ** 3) * np.heaviside(x - x_3, 1) +
                                       forces[7].z * ((x - x_a1) ** 3) * np.heaviside(x - x_a1, 1) +
                                       forces[9].z * ((x - x_a2) ** 3) * np.heaviside(x - x_a2, 1)) - qz / 24 * (
                                     x ** 4) + constant[0] * x + constant[1])


d_y = eq_def_y(x_slice, int_constant_y)
d_z = eq_def_z(x_slice, int_constant_z)
plt.plot(x_slice, d_z)
plt.plot(x_slice, d_y)


def absolute_def(span_defy, span_defz, twist, shear_center):
    dy_le = span_defy - h_a/2*np.sin(np.radians(theta)) - (h_a/2-shear_center[2])*np.sin(twist)
    dy_te = span_defy + (c_a-h_a/2)*np.sin(np.radians(theta)) + (c_a-h_a/2+shear_center[2])*np.sin(twist)
    dz_le = span_defz 
    dz_te = span_defz
    return dy_le, dy_te, dz_le, dz_te

def normal_stress(momenty, momentz, moi, booms_geometry):
    sigma_x = np.zeros((len(momenty),len(booms_geometry)))
    for i in range(len(momenty)):
        for j in range(len(booms_geometry)):
            sigma_x[i][j] = momenty[i]/moi[0]*booms_geometry[1] + momentz[i]/moi[1]*booms_geometry[0]
    return sigma_x


def von_mises_stres(normal_stress, shear_stress_1):
    vonmis_stress = np.zeros(np.shape(normal_stress))
    for i in range(len(normal_stress[0])):
        for j in range(len(normal_stress[1])):
            vonmis_stress[i][j] = np.sqrt(0.5*(normal_stress[i][j]*normal_stress[i][j])+
                         3*shear_stress_1[i][j]*shear_stress_1[i][j])
    return vonmis_stress

    




