import numpy as np


def shear_torque(torque, area):
    shear_flow = torque / (2 * area)

    return shear_flow


def q_b(booms, a, i_zz, i_yy, v_y, v_z):
    diff_q = np.zeros(len(booms))

    alpha = -v_y / i_zz
    beta = -v_z / i_yy
    value = 0

    for i in range(len(booms)):
        delta_q = alpha * booms[i][1] * a[i] + beta * booms[i][0] * a[i]
        value += delta_q
        diff_q[i] = value

    base_q = np.zeros(len(booms))
    start = 0
    for i in range(1, len(diff_q)):
        base_q[i] = diff_q[i] - start
        start = diff_q[i]

    return base_q


def q_s0(base_q, t, s, shear_mod, d_theta, area):
    flow = np.multiply(base_q, s)
    flow = np.divide(flow, t)

    num = np.sum(flow) + 2*area*shear_mod*d_theta
    x = np.sum(np.divide(s, t))
    y = -1 * (s[-1] / t[-1])

    return x, y, num


def cell_separation(skin, base_q, s, n_top, n_bottom, d_spar, spar, h_a, c_a, booms):
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


def get_shear_flows(b, areas, i_zz, t_skin, d, top, bottom, h_a, t_spar, c_a, i_yy, v_y, v_z, shear_mod, d_theta):
    q_diff = q_b(b, areas, i_zz, i_yy, v_y, v_z)

    celli, cellii = cell_separation(t_skin, q_diff, d, top, bottom, h_a, t_spar, h_a, c_a, b)

    qsi = np.array([q_s0(celli[0], celli[2], celli[1], shear_mod, d_theta, celli[4])])
    qsii = np.array([q_s0(cellii[0], cellii[2], cellii[1], shear_mod, d_theta, cellii[4])])

    sys = np.array([[qsi[0][0], qsi[0][1]], [qsii[0][0], qsii[0][1]]])
    sys_vec = np.array([[-qsi[0][2]], [-qsii[0][2]]])

    values = np.linalg.solve(sys, sys_vec)

    return values, celli, cellii







