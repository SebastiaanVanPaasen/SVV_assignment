import numpy as np
from SVV_assignment.SVV_assignment.geometry import Idealization as Idealization


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


def cell_separation(skin, base_q, s, n_top, n_bottom, d_spar, spar, booms):
    base_q_cellI = np.append(base_q[n_bottom:], base_q[:n_top])
    base_q_cellI = np.append(base_q_cellI, np.array([0]))

    bottom_locations = s[n_bottom:]
    top_locations = s[:n_top + 1]

    s_cellI_bottom = np.zeros(len(bottom_locations) - 1)
    for i in range(len(bottom_locations) - 1):
        s_cellI_bottom[i] = bottom_locations[i + 1] - bottom_locations[i]
    s_cellI_bottom[-1] += top_locations[0]

    s_cellI_top = np.zeros(len(top_locations) - 1)
    for i in range(len(top_locations) - 1):
        s_cellI_top[i] = top_locations[i + 1] - top_locations[i]

    s_cellI = np.append(np.append(s_cellI_bottom, s_cellI_top), np.array([d_spar]))
    t_cellI = np.append(np.array((len(s_cellI) - 1) * [skin]), np.array([spar]))

    cellI_coordinates = booms[n_bottom:] + booms[:n_top + 1]

    base_q_cellII = np.append(base_q[n_top:n_bottom], np.array([0]))

    s_cellII_locations = s[n_top:n_bottom + 1]
    s_cellII = np.zeros(len(s_cellII_locations) - 1)
    for i in range(len(s_cellII_locations) - 1):
        s_cellII[i] = s_cellII_locations[i + 1] - s_cellII_locations[i]

    s_cellII = np.append(s_cellII, np.array([d_spar]))
    t_cellII = np.append(np.array((len(s_cellII) - 1) * [skin]), np.array([spar]))

    cellII_coordinates = booms[n_top:n_bottom + 1]

    return [base_q_cellI, s_cellI, t_cellI, cellI_coordinates], [base_q_cellII, s_cellII, t_cellII, cellII_coordinates]


def area_triangle(pos1, pos2, pos3):
    len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
    len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
    len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)

    s = (len1 + len2 + len3) / 2
    area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))

    return area


def area_circle(distance):
    ratio = distance / (h_a * np.pi / 2)

    area = ratio * np.pi * ((h_a / 2) ** 2) / 2

    return area


def sum_moments_I(base_flow, constant_flow, cell_area, s):
    moment = 0

    for i in range(len(s) - 1):
        area = area_circle(s[i])
        moment += 2 * area * base_flow[i]

    moment += 2 * cell_area * constant_flow

    return moment


def sum_moments_II(base_flow, constant_flow, cell_area, pos, summation_point):
    moment = 0
    flow_areas = []

    for i in range(len(pos) - 1):
        flow_area = area_triangle(summation_point, pos[i], pos[i + 1])
        flow_areas.append(flow_area)

    midpoint = len(areas) // 2

    flow_areas[midpoint] += area_triangle(pos[midpoint], pos[midpoint + 1], [c_a, 0])

    for i in range(len(base_flow) - 1):
        moment += 2 * flow_areas[i] * base_flow[i]

    moment += 2 * constant_flow * cell_area

    return moment


# pre shear-flow calculations
ideal = Idealization(20)
b, d = ideal.set_boom_locations()
areas = ideal.calculate_boom_area(b)
b, d, areas, top, bottom = ideal.add_spar_booms(b, d, areas)
z_pos, b = ideal.calculate_cg(b, areas)
i_zz, i_yy = ideal.calculate_moi(b, areas)

# required parameters
h_a = 0.173
c_a = 0.484
t_skin = 1.1 / 1000
t_spar = 2.5 / 1000

# shear-flows
q_diff = q_b(b, areas, i_zz)

cellI, cellII = cell_separation(t_skin, q_diff, d, top, bottom, h_a, t_spar, b)

qsI = np.array([q_s0(cellI[0], cellI[2], cellI[1])])
qsII = np.array([q_s0(cellII[0], cellII[2], cellII[1])])

sys = np.array([[qsI[0][0], qsI[0][1]], [qsII[0][0], qsII[0][1]]])
sys_vec = np.array([[qsI[0][2]], [qsII[0][2]]])

values = np.linalg.solve(sys, sys_vec)

# summing the moments
areaI = np.pi * ((h_a / 2) ** 2) / 2
areaII = h_a * (c_a - (h_a / 2))

moments_I = sum_moments_I(cellI[0], values[0][0], areaI, cellI[1])
moments_II = sum_moments_II(cellII[0], values[1][0], areaII, cellII[3], [h_a / 2, 0])

z_pos_sc = ((moments_I + moments_II) / -1 + (h_a / 2)) - z_pos

print(z_pos)
print(z_pos_sc)
