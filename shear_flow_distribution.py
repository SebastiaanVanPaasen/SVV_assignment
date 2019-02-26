import math
import numpy as np


class Torque:

    def __int__(self):
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

        return self.cell_i, self.cell_ii

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

        sys = np.array(torque_eq[:3], d_theta_i[:3], d_theta_ii[:3])
        sys_vec = np.array(torque_eq[3], d_theta_i[3], d_theta_ii[3])

        values = np.linalg.solve(sys, sys_vec)

        return values


class Shear:

    def __int__(self):
        self.cell_i = []
        self.cell_ii = []
        self.h_a = 0.173
        self.c_a = 0.484
        self.t_skin = 1.1 / 1000
        self.t_spar = 2.5 / 1000
        self.shear_mod = 23 * 10 ** 9

    def create_cell_properties(self, s, b, n_top, n_bottom, base_q):
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

        cell_i_coordinates = b[n_bottom:] + b[:n_top + 1]
        self.cell_i.append(cell_i_coordinates)

        s_cell_ii_locations = s[n_top:n_bottom + 1]
        s_cell_ii = np.zeros(len(s_cell_ii_locations) - 1)
        for i in range(len(s_cell_ii_locations) - 1):
            s_cell_ii[i] = s_cell_ii_locations[i + 1] - s_cell_ii_locations[i]

        s_cell_ii = np.append(s_cell_ii, np.array([self.h_a]))
        self.cell_ii.append(s_cell_ii)
        t_cell_ii = np.append(np.array((len(s_cell_ii) - 1) * [self.t_skin]), np.array([self.t_spar]))
        self.cell_ii.append(t_cell_ii)

        cell_ii_coordinates = b[n_top:n_bottom + 1]
        self.cell_ii.append(cell_ii_coordinates)

        area_i = math.pi * ((self.h_a / 2) ** 2) / 2
        self.cell_i.append(area_i)
        area_ii = self.h_a * (self.c_a - (self.h_a / 2))
        self.cell_ii.append(area_ii)

        base_q_cell_i = np.append(base_q[n_bottom:], base_q[:n_top])
        base_q_cell_i = np.append(base_q_cell_i, np.array([0]))
        base_q_cell_ii = np.append(base_q[n_top:n_bottom], np.array([0]))

        self.cell_i.append(base_q_cell_i)
        self.cell_ii.append(base_q_cell_ii)

    @staticmethod
    def q_b(booms, areas, i_zz, i_yy, v_y, v_z):
        diff_q = np.zeros(len(booms))

        alpha = -v_y / i_zz
        beta = -v_z / i_yy
        value = 0

        for i in range(len(booms)):
            delta_q = alpha * booms[i][1] * areas[i] + beta * booms[i][0] * areas[i]
            value += delta_q
            diff_q[i] = value

        base_q = np.zeros(len(booms))
        start = 0
        for i in range(1, len(diff_q)):
            base_q[i] = diff_q[i] - start
            start = diff_q[i]

        return base_q

    def d_theta(self, cell):
        flow = np.multiply(cell[4], cell[0])
        flow = np.divide(flow, cell[1])

        num = -np.sum(flow)
        x = np.sum(np.divide(cell[0], cell[1]))
        y = -1 * (cell[0][-1] / cell[1][-1])
        z = -2 * cell[4] * self.shear_mod

        return np.array([x, y, z, num])

    @staticmethod
    def area_triangle(pos1, pos2, pos3):
        len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
        len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
        len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)

        s = (len1 + len2 + len3) / 2
        area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))

        return area

    @staticmethod
    def area_circle(distance, h_a):
        ratio = distance / (h_a * np.pi / 2)
        area = ratio * np.pi * ((h_a / 2) ** 2) / 2

        return area

    def sum_moments_i(self):
        moment = 0

        for i in range(len(self.cell_i[0]) - 1):
            area = self.area_circle(self.cell_i[0], self.h_a)
            moment += 2 * area * self.cell_i[4][i]

        return moment

    def sum_moments_ii(self, summation_point):
        moment = 0
        flow_areas = []

        for i in range(len(self.cell_ii) - 1):
            flow_area = self.area_triangle(summation_point, self.cell_ii[2][i], self.cell_ii[2][i + 1])
            flow_areas.append(flow_area)

        midpoint = len(flow_areas) // 2

        flow_areas[midpoint] += self.area_triangle(self.cell_ii[2][midpoint], self.cell_ii[2][midpoint + 1],
                                                   [self.c_a, 0])

        for i in range(len(self.cell_ii[4]) - 1):
            moment += 2 * flow_areas[i] * self.cell_ii[4][i]

        return moment

    def solve_system(self, summation_point):
        moments = np.array(
            [2 * self.cell_i[3], 2 * self.cell_ii[3], 0, -self.sum_moments_i() - self.sum_moments_ii(summation_point)])
        d_theta_i = self.d_theta(self.cell_i)
        d_theta_ii = self.d_theta(self.cell_ii)

        sys = np.array([moments[:3], d_theta_i[:3], d_theta_ii[:3]])
        sys_vec = np.array([[moments[3]], [d_theta_i[3]], [d_theta_ii[3]]])

        values = np.linalg.solve(sys, sys_vec)

        return values
