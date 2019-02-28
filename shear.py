from geometry import *
from shear_center import *

from geometry import *
from shear_center import *


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
        Creates the cell properties of the two closed cells. The circumferential distance parameter s starts at the
        bottom of the spar for the first cell and at the top of the spar for the second cell. Both going in counter
        clockwise direction.

        :param s: List of floats, containing the positions of all the booms around the cross-section.
        :param b: List of float lists, containing all the positions of the booms around the cross-section.
        :param n_top: Integer, containing the boom number of the boom at the top of the spar.
        :param n_bottom: Integer, containing the boom number of the boom at the bottom of the spar.
        :param base_q: List of floats, containing the base shear flows between all the booms.
        :return: Two lists, containing all the required parameters of each cell.
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
        """
        Calculates the base shear flow due to shear forces in both y- and z-direction.

        :param booms: List of float lists, containing the locations of all the booms in the cross-section.
        :param areas: List of floats, containing the areas of all the booms in the cross-section.
        :param i_zz: Float, value of the moment of inertia around the z-axis.
        :param i_yy: Float, value of the moment of inertia around the y-axis.
        :param v_y: Float, containing the value of the shear force in y-direction.
        :param v_z: Float, containing the value of the shear force in z-direction.
        :return: List of floats containing the base shear flows between all the booms.
        """
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
        """
        Calculates the angle of twist equation of a cell

        :param cell: List of parameters, containing the circumferential distances between booms, the thickness of each
        part, the base shear flows and the total area of the cell.
        :return: The values of the constant shear flows of the two cells in one of the cells, the angle of twist and the
        twist created by the base shear flow.
        """
        flow = np.multiply(cell[4], cell[0])
        flow = np.divide(flow, cell[1])

        num = -np.sum(flow)
        x = np.sum(np.divide(cell[0], cell[1]))
        y = -1 * (cell[0][-1] / cell[1][-1])
        z = -2 * cell[3] * self.shear_mod

        return np.array([x, y, z, num])

    @staticmethod
    def area_triangle(pos1, pos2, pos3):
        """
        Calculates the enclosed area between three positions in the yz-coordinate frame.

        :param pos1: First corner of the triangle.
        :param pos2: Second corner of the triangle.
        :param pos3: Third corner of the triangle.
        :return: A float representing the area enclosed by the triangle.
        """
        len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
        len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
        len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)

        s = (len1 + len2 + len3) / 2
        area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))

        return area

    @staticmethod
    def area_circle(distance, h_a):
        """
        Calculates the enclosed area between two points on the semi-circle at the front of the cross-section and the
        middle of the spar.

        :param distance: Float, representing the circumferential distance along the semi-circle between two booms.
        :param h_a: Float, representing the height of the spar
        :return: Float, representing the area enclosed by the circumferential distance and the middle of the spar.
        """
        ratio = distance / (h_a * np.pi / 2)
        area = ratio * np.pi * ((h_a / 2) ** 2) / 2

        return area

    def sum_moments_i(self):
        """
        Calculates the moment generated by all the base shear flows in the nose-cell.

        :return: Float, representing the total moment created by all the base shear flows.
        """
        moment = 0

        for i in range(len(self.cell_i[0]) - 1):
            area = self.area_circle(self.cell_i[0][i], self.h_a)
            moment += 2 * area * self.cell_i[4][i]

        return moment

    def sum_moments_ii(self, summation_point):
        """
        Calculates the moment generated by all the base shear flows in the second cell.

        :param summation_point: The yz-point around which the moments are summed.
        :return: Float, representing the total moment generated by the base shear flows in the second cell.
        """
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

    def solve_system(self, summation_point, v_y, sc_z):
        """
        Solves for the angle of twist and the constant shear flows in both cells.

        :param summation_point: The yz-point around which the moments in the second cell will be summed.
        :return: A 1D array containing the constant shear flows and the angle of twist.
        """
        moments_i = -self.sum_moments_i()
        moments_ii = - self.sum_moments_ii(summation_point)
        moments_shear_force = abs(summation_point[2] - sc_z)*v_y
        total_moment = moments_i + moments_ii + moments_shear_force

        moments = np.array([2 * self.cell_i[3], 2 * self.cell_ii[3], 0, total_moment])
        d_theta_i = self.d_theta(self.cell_i)
        d_theta_ii = self.d_theta(self.cell_ii)

        sys = np.array([moments[:3], d_theta_i[:3], d_theta_ii[:3]])
        sys_vec = np.array([[moments[3]], [d_theta_i[3]], [d_theta_ii[3]]])

        values = np.linalg.solve(sys, sys_vec)

        return values


# Simple get function to obtain the final shear flow values and the angle of twist

def get_shear_flow(n_booms, v_y, v_z, position):
    shear = Shear()
    cg_z, boom_locations = get_cg(n_booms)
    moi_zz, moi_yy = get_moi(n_booms)
    booms, distances, boom_areas, top, bottom = get_boom_information(n_booms)
    shear_center = get_shear_center()

    q_b = shear.q_b(boom_locations, boom_areas, moi_zz, moi_yy, v_y, v_z)
    cell_i, cell_ii = shear.create_cell_properties(distances, boom_locations, top, bottom, q_b)

    final_shear_values = shear.solve_system(position, v_y, shear_center)
    q_base_i = np.array(cell_i[4])
    q_base_ii = np.array(cell_ii[4])

    q_flow_i = np.add(q_base_i, final_shear_values[0])
    q_flow_ii = np.add(q_base_ii, final_shear_values[1])

    q_flow_i[-1] = q_flow_i[-1] - final_shear_values[1]
    q_flow_ii[-1] = q_flow_ii[-1] - final_shear_values[0]

    return q_flow_i, q_flow_ii, final_shear_values[2], cell_i[2], cell_ii[2]
