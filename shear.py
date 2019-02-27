from SVV_assignment.SVV_assignment.geometry import *

class Shear:

    def __init__(self):
        self.cell_i = []
        self.cell_ii = []
        self.h_a = 0.173
        self.c_a = 0.484
        self.t_skin = 1.1 / 1000
        self.t_spar = 2.5 / 1000
        self.shear_mod = 23 * 10 ** 9

    def create_cell_properties(self, s, b, n_top, n_bottom, base_q):
        """

        :param s:
        :param b:
        :param n_top:
        :param n_bottom:
        :param base_q:
        :return:
        """
        bottom_locations = s[n_bottom:]
        top_locations = s[:n_top + 1]

        s_cell_i_bottom = np.zeros(len(bottom_locations))
        for i in range(len(bottom_locations) - 1):
            s_cell_i_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
        s_cell_i_bottom[-1] += top_locations[0]*2

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

        base_q_cell_i = np.append(base_q[n_bottom + 1:], base_q[:n_top + 1])
        base_q_cell_i = np.append(base_q_cell_i, np.array([0]))
        base_q_cell_ii = np.append(base_q[n_top + 1:n_bottom + 1], np.array([0]))

        self.cell_i.append(base_q_cell_i)
        self.cell_ii.append(base_q_cell_ii)

        return self.cell_i, self.cell_ii

    @staticmethod
    def q_b(booms, areas, i_zz, i_yy, v_y, v_z):
        diff_q = np.zeros(len(booms))
        alpha = -v_y / i_zz
        beta = -v_z / i_yy

        for i in range(len(booms)):
            delta_q = alpha * booms[i][1] * areas[i] + beta * booms[i][0] * areas[i]
            diff_q[i] = delta_q

        base_q = np.zeros(len(booms))
        previous = diff_q[0]
        for i in range(1, len(diff_q)):
            base_q[i] = previous + diff_q[i]
            previous = diff_q[i]

        return base_q

    def d_theta(self, cell):
        flow = np.multiply(cell[4], cell[0])
        flow = np.divide(flow, cell[1])

        num = -np.sum(flow)
        x = np.sum(np.divide(cell[0], cell[1]))
        y = -1 * (cell[0][-1] / cell[1][-1])
        z = -2 * cell[3] * self.shear_mod

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
            area = self.area_circle(self.cell_i[0][i], self.h_a)
            moment += 2 * area * self.cell_i[4][i]

        return moment

    def sum_moments_ii(self, summation_point):
        moment = 0
        flow_areas = []

        for i in range(len(self.cell_ii[2]) - 1):
            flow_area = self.area_triangle(summation_point, self.cell_ii[2][i], self.cell_ii[2][i + 1])
            flow_areas.append(flow_area)

        flow_areas.append(0)
        midpoint = len(flow_areas) // 2

        flow_areas[midpoint] += self.area_triangle(self.cell_ii[2][midpoint], self.cell_ii[2][midpoint + 1],
                                                   [self.c_a, 0])

        for i in range(len(self.cell_ii[4])):
            moment += 2 * flow_areas[i] * self.cell_ii[4][i]

        return moment

    def solve_system(self, summation_point):
        moments_i = -self.sum_moments_i()
        moments_ii = - self.sum_moments_ii(summation_point)
        total_moment = moments_i + moments_ii

        moments = np.array([2 * self.cell_i[3], 2 * self.cell_ii[3], 0, total_moment])
        d_theta_i = self.d_theta(self.cell_i)
        d_theta_ii = self.d_theta(self.cell_ii)

        sys = np.array([moments[:3], d_theta_i[:3], d_theta_ii[:3]])
        sys_vec = np.array([[moments[3]], [d_theta_i[3]], [d_theta_ii[3]]])

        values = np.linalg.solve(sys, sys_vec)

        return values


def get_shear_flow(n_booms, v_y, v_z, position):
    shear = Shear()
    cg_z, boom_locations = get_cg(n_booms)
    moi_zz, moi_yy = get_moi(n_booms)
    booms, distances, boom_areas, top, bottom = get_boom_information(n_booms)

    q_b = shear.q_b(boom_locations, boom_areas, moi_zz, moi_yy, v_y, v_z)
    cell_i, cell_ii = shear.create_cell_properties(distances, boom_locations, top, bottom, q_b)

    final_shear_values = shear.solve_system(position)
    q_base_i = np.array(cell_i[4])
    q_base_ii = np.array(cell_ii[4])

    q_flow_i = np.add(q_base_i, final_shear_values[0])
    q_flow_ii = np.add(q_base_ii, final_shear_values[1])

    return q_flow_i, q_flow_ii, final_shear_values[2]
