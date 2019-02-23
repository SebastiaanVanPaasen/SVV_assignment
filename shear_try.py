import numpy as np
from SVV_assignment.SVV_assignment.geometry import *


def q_b(booms, a, moi):
    diff_q = np.zeros(len(booms))

    alpha = -1 / moi
    value = 0
    for i in range(len(booms)):
        delta_q = alpha * booms[i][1] * a[i]
        value += delta_q
        diff_q[i] = value

    base_q = np.zeros(len(booms))
    start = 0
    for i in range(1, len(diff_q)):
        base_q[i] = diff_q[i] - start
        start = diff_q[i]

    return base_q


def q_s0(base_q, t, s):
    flow = np.multiply(base_q, s)
    flow = np.divide(flow, t)

    num = np.sum(flow)
    x = np.sum(np.divide(s, t))
    y = -1 * (s[-1] / t[-1])

    return x, y, num


def cell_separation(skin, base_q, s, n_top, n_bottom, d_spar, spar, h_a, c_a):
    base_q_celli = np.append(base_q[n_bottom:], base_q[:n_top])
    base_q_celli = np.append(base_q_celli, np.array([0]))

    bottom_locations = s[n_bottom:]
    top_locations = s[:n_top + 1]

    s_celli_bottom = np.zeros(len(bottom_locations) - 1)
    for i in range(len(bottom_locations) - 1):
        s_celli_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
    s_celli_bottom[-1] += top_locations[0]

    s_celli_top = np.zeros(len(top_locations) - 1)
    for i in range(len(top_locations) - 1):
        s_celli_top[i] = top_locations[i + 1] - top_locations[i]

    s_celli = np.append(np.append(s_celli_bottom, s_celli_top), np.array([d_spar]))
    t_celli = np.append(np.array((len(s_celli) - 1) * [skin]), np.array([spar]))

    celli_coordinates = booms[n_bottom:] + booms[:n_top + 1]

    base_q_cellii = np.append(base_q[n_top:n_bottom], np.array([0]))

    s_cellii_locations = s[n_top:n_bottom + 1]
    s_cellii = np.zeros(len(s_cellii_locations) - 1)
    for i in range(len(s_cellii_locations) - 1):
        s_cellii[i] = s_cellii_locations[i + 1] - s_cellii_locations[i]

    s_cellii = np.append(s_cellii, np.array([d_spar]))
    t_cellii = np.append(np.array((len(s_cellii) - 1) * [skin]), np.array([spar]))

    cellii_coordinates = booms[n_top:n_bottom + 1]

    areai = np.pi * ((h_a / 2) ** 2) / 2
    areaii = h_a * (c_a - (h_a / 2))

    return [base_q_celli, s_celli, t_celli, celli_coordinates, areai], \
           [base_q_cellii, s_cellii, t_cellii, cellii_coordinates, areaii]


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
    
    celli, cellii = cell_separation(t_skin, q_diff, d, top, bottom, h_a, t_spar, b, h_a, c_a)
    
    qsi = np.array([q_s0(celli[0], celli[2], celli[1])])
    qsii = np.array([q_s0(cellii[0], cellii[2], cellii[1])])
    
    sys = np.array([[qsi[0][0], qsi[0][1]], [qsii[0][0], qsii[0][1]]])
    sys_vec = np.array([[qsi[0][2]], [qsii[0][2]]])
    
    values = np.linalg.solve(sys, sys_vec)
    
    return values, celli, cellii


def summing_moments(constants, celli, cellii, z_pos):
    moments_i = sum_moments_i(celli[0], constants[0][0], celli[4], celli[1], height)
    moments_ii = sum_moments_ii(cellii[0], constants[1][0], cellii[4], cellii[3], [height / 2, 0], chord)

    z_pos_sc = ((moments_i + moments_ii) / -1 + (height / 2)) - z_pos

    return z_pos_sc


# pre shear-flow calculations
n_booms = 20
boom_distances, boom_areas, top_spar, bottom_spar = get_boom_information(n_booms)
z_pos, booms = get_cg(n_booms)
izz, iyy = get_moi(n_booms)

# required parameters
height = 0.173
chord = 0.484
skin_t = 1.1 / 1000
spar_t = 2.5 / 1000


def get_shear_center():
    constant, first_cell, second_cell = get_shear_flows(booms, boom_areas, izz, skin_t, boom_distances, top_spar,
                                                        bottom_spar, height, spar_t, chord)

    shear_center = summing_moments(constant, first_cell, second_cell, z_pos)

    return shear_center
