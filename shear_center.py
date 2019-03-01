from geometry import *
import numpy as np

from geometry import *


def q_b(boom_locations, area, moi):
    """
    Calculates the base shear flow due to a dummy force of 1 N.

    :param boom_locations: List of float lists, containing all the locations of the booms around the cross-section.
    :param area: List of floats, containing the areas of all the booms around the cross-section.
    :param moi: Float, representing the moment of inertia around the z-axis.
    :return: List of floats, containing the base shear flows around the cross-section.
    """
    diff_q = np.zeros(len(boom_locations))

    alpha = -1 / moi
    for i in range(len(boom_locations)):
        delta_q = alpha * boom_locations[i][1] * area[i]
        diff_q[i] = delta_q

    base_q = np.zeros(len(boom_locations))
    previous = diff_q[0]
    for i in range(1, len(diff_q)):
        base_q[i] = previous + diff_q[i]
        previous = diff_q[i]

    # print(base_q)
    # print(len(base_q))
    return base_q


def q_s0(base_q, t, s):
    """
    Calculates the constant shear flow equations leaving them as unknowns in the equations.

    :param base_q: List of floats, containing the base shear flows between all the booms in one cell.
    :param t: List of floats, containing the thickness between all the booms in one cell.
    :param s: List of floats, containing the distances between the booms in one cell.
    :return: values for the two unknown shear flows and the numerical value due to the base shear flow.
    """
    flow = np.multiply(base_q, s)
    flow = np.divide(flow, t)

    num = np.sum(flow)
    x = np.sum(np.divide(s, t))
    y = -1 * (s[-1] / t[-1])

    return x, y, num


def cell_separation(skin, base_q, s, n_top, n_bottom, d_spar, spar, h_a, c_a, booms):
    """
    Creates the cell parameters for both cells.

    :param skin: Float, representing the thickness of the skin.
    :param base_q: List of floats, representing all the base shear flows around the section.
    :param s: List of floats, representing the circumferential positions of the booms around the section.
    :param n_top: Integer, boom number of the boom at the top of the spar.
    :param n_bottom: Integer, boom number of the boom at the bottom of the spar.
    :param d_spar: Float, circumferential distance to the boom at the top of the spar.
    :param spar: Float, representing the thickness of the spar.
    :param h_a: Float, representing the height of the spar.
    :param c_a: Float, representing the length of the chord.
    :param booms: List of float lists, containing the locations of all the booms in the yz-coordinate frame.
    :return: Two lists, containing the base shear flow, circumferential distances between booms, thickness of the
    distances between the booms, the area of the cell and the coordinates of all the booms in the cell.
    """
    base_q_cell_i = np.append(base_q[n_bottom + 1:], base_q[:n_top + 1])
    base_q_cell_i = np.append(base_q_cell_i, np.array([0]))
    base_q_cell_ii = np.append(base_q[n_top + 1:n_bottom + 1], np.array([0]))

    bottom_locations = s[n_bottom:]
    top_locations = s[:n_top + 1]

    s_cell_i_bottom = np.zeros(len(bottom_locations))
    for i in range(len(bottom_locations) - 1):
        s_cell_i_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
    s_cell_i_bottom[-1] += top_locations[0] * 2

    s_cell_i_top = np.zeros(len(top_locations) - 1)
    for i in range(len(top_locations) - 1):
        s_cell_i_top[i] = top_locations[i + 1] - top_locations[i]

    s_cell_i = np.append(np.append(s_cell_i_bottom, s_cell_i_top), np.array([d_spar]))

    s_cell_ii_locations = s[n_top:n_bottom + 1]
    s_cell_ii = np.zeros(len(s_cell_ii_locations) - 1)
    for i in range(len(s_cell_ii_locations) - 1):
        s_cell_ii[i] = s_cell_ii_locations[i + 1] - s_cell_ii_locations[i]

    s_cell_ii = np.append(s_cell_ii, np.array([d_spar]))

    t_cell_i = np.append(np.array((len(s_cell_i) - 1) * [skin]), np.array([spar]))
    t_cell_ii = np.append(np.array((len(s_cell_ii) - 1) * [skin]), np.array([spar]))

    cell_i_coordinates = booms[n_bottom:] + booms[:n_top + 1]
    cell_ii_coordinates = booms[n_top:n_bottom + 1]

    area_i = np.pi * ((h_a / 2) ** 2) / 2
    area_ii = h_a * (c_a - (h_a / 2))

    return [base_q_cell_i, s_cell_i, t_cell_i, cell_i_coordinates, area_i], \
           [base_q_cell_ii, s_cell_ii, t_cell_ii, cell_ii_coordinates, area_ii]


def area_triangle(pos1, pos2, pos3):
    len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
    len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
    len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)

    s = (len1 + len2 + len3) / 2
    area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))

    return area


def area_circle(distance, h_a):
    ratio = distance / (h_a * np.pi / 2)

    area = ratio * np.pi * ((h_a / 2) ** 2) / 2

    return area


def sum_moments_i(base_flow, constant_flow, cell_area, s, h_a):
    moment = 0

    for i in range(len(s) - 1):
        area = area_circle(s[i], h_a)
        moment += 2 * area * base_flow[i]

    moment += 2 * cell_area * constant_flow

    return moment


def sum_moments_ii(base_flow, constant_flow, cell_area, pos, summation_point, c_a):
    moment = 0
    flow_areas = []

    for i in range(len(pos) - 1):
        flow_area = area_triangle(summation_point, pos[i], pos[i + 1])
        flow_areas.append(flow_area)

    midpoint = len(flow_areas) // 2

    flow_areas[midpoint] += area_triangle(pos[midpoint], pos[midpoint + 1], [c_a, 0])

    for i in range(len(base_flow) - 1):
        moment += 2 * flow_areas[i] * base_flow[i]

    moment += 2 * constant_flow * cell_area

    return moment


def get_shear_flows(b, areas, i_zz, t_skin, d, top, bottom, h_a, t_spar, c_a):
    q_diff = q_b(b, areas, i_zz)

    cell_i, cell_ii = cell_separation(t_skin, q_diff, d, top, bottom, h_a, t_spar, h_a, c_a, b)

    qs_i = np.array([q_s0(cell_i[0], cell_i[2], cell_i[1])])
    qs_ii = np.array([q_s0(cell_ii[0], cell_ii[2], cell_ii[1])])

    sys = np.array([[qs_i[0][0], qs_i[0][1]], [qs_ii[0][0], qs_ii[0][1]]])
    sys_vec = np.array([[-qs_i[0][2]], [-qs_ii[0][2]]])

    values = np.linalg.solve(sys, sys_vec)

    return values, cell_i, cell_ii


def summing_moments(constants, cell_i, cell_ii, z_pos):
    moments_i = sum_moments_i(cell_i[0], constants[0][0], cell_i[4], cell_i[1], height)
    moments_ii = sum_moments_ii(cell_ii[0], constants[1][0], cell_ii[4], cell_ii[3], [height / 2, 0], chord)

    z_pos_sc = ((moments_i + moments_ii) / -1)

    return z_pos_sc


# pre shear-flow calculations
n_booms = 260
b, boom_distances, boom_areas, top_spar, bottom_spar = get_boom_information(n_booms)
z_pos, boom_loc = get_cg(n_booms)

izz, iyy = get_moi(n_booms)
# required parameters
height = 0.173
chord = 0.484
skin_t = 1.1 / 1000
spar_t = 2.5 / 1000


def get_shear_center():
    constant, first_cell, second_cell = get_shear_flows(boom_loc, boom_areas, izz, skin_t, boom_distances,
                                                        top_spar,
                                                        bottom_spar, height, spar_t, chord)

    shear_center = summing_moments(constant, first_cell, second_cell, z_pos)

    return shear_center


sc_z = get_shear_center()
# print(sc_z)
