import math
import numpy as np
from SVV_assignment.SVV_assignment.geometry import *


class Torque:

    def __init__(self):
        self.cell_i = []
        self.cell_ii = []
        self.h_a = 0.173
        self.c_a = 0.484
        self.t_skin = 1.1 / 1000
        self.t_spar = 2.5 / 1000
        self.shear_mod = 23 * 10 ** 9

    def create_cell_properties(self, s, n_top, n_bottom):
        bottom_locations = s[n_bottom:]
        top_locations = s[:n_top + 1]

        s_cell_i_bottom = np.zeros(len(bottom_locations) - 1)
        for i in range(len(bottom_locations) - 1):
            s_cell_i_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
        s_cell_i_bottom[-1] += top_locations[0]

        s_cell_i_top = np.zeros(len(top_locations) - 1)
        for i in range(len(top_locations) - 1):
            s_cell_i_top[i] = top_locations[i + 1] - top_locations[i]

        s_cell_i = np.append(np.append(s_cell_i_bottom, s_cell_i_top), np.array([self.h_a]))
        t_cell_i = np.append(np.array((len(s_cell_i) - 1) * [self.t_skin]), np.array([self.t_spar]))
        self.cell_i.append(s_cell_i)
        self.cell_i.append(t_cell_i)

        s_cell_ii_locations = s[n_top:n_bottom + 1]
        s_cell_ii = np.zeros(len(s_cell_ii_locations) - 1)
        for i in range(len(s_cell_ii_locations) - 1):
            s_cell_ii[i] = s_cell_ii_locations[i + 1] - s_cell_ii_locations[i]

        s_cell_ii = np.append(s_cell_ii, np.array([self.h_a]))
        self.cell_ii.append(s_cell_ii)
        t_cell_ii = np.append(np.array((len(s_cell_ii) - 1) * [self.t_skin]), np.array([self.t_spar]))
        self.cell_ii.append(t_cell_ii)

        area_i = math.pi * ((self.h_a / 2) ** 2) / 2
        self.cell_i.append(area_i)
        area_ii = self.h_a * (self.c_a - (self.h_a / 2))
        self.cell_ii.append(area_ii)

    def q_t(self, cell):
        x = np.sum(np.divide(cell[0], cell[1]))
        y = -1 * (cell[0][-1] / cell[1][-1])
        z = -2 * cell[2] * self.shear_mod
        num = 0

        return np.array([x, y, z, num])

    def shear_torque(self, torque):
        torque_eq = np.array([2 * self.cell_i[2], 2 * self.cell_ii[2], 0, torque])

        return torque_eq

    def solve_system(self, torque):
        torque_eq = self.shear_torque(torque)
        d_theta_i = self.q_t(self.cell_i)
        d_theta_ii = self.q_t(self.cell_ii)

        sys = np.array([torque_eq[:3], d_theta_i[:3], d_theta_ii[:3]])
        sys_vec = np.array([[torque_eq[3]], [d_theta_i[3]], [d_theta_ii[3]]])

        values = np.linalg.solve(sys, sys_vec)

        return values


def get_torque(n_booms, m_x):
    torque = Torque()
    booms, distances, boom_areas, top, bottom = get_boom_information(n_booms)

    torque.create_cell_properties(distances, top, bottom)
    final_torque_values = torque.solve_system(m_x)

    q_flow_i = final_torque_values[0][0]
    q_flow_ii = final_torque_values[0][1]
    d_theta = final_torque_values[0][2]

    return q_flow_i, q_flow_ii, d_theta
