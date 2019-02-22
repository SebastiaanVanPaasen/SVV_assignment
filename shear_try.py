import numpy as np
from SVV_assignment.SVV_assignment.geometry import Idealization as Idealization

ideal = Idealization(200)

b, d = ideal.set_boom_locations()

areas = ideal.calculate_boom_area(b)

b, d, areas = ideal.add_spar_booms(b, d, areas)

z_pos, b = ideal.calculate_cg(b, areas)

i_zz, i_yy = ideal.calculate_moi(b, areas)


h_a = 0.173
c_a = 0.484

n_I = 149
n_spar = 75
n_II = 453

t_skin = 1.1 / 1000
t_spar = 2.5 / 1000

s_I = np.array(n_I * [s[1] - s[0]])
s_I[-1] = h_a

s_II = np.array(n_II * [s[1] - s[0]])
s_II[-1] = h_a

t_I = np.array(n_I * [t_skin])
t_II = np.array(n_II * [t_skin])
t_I[-1] = t_spar
t_II[-1] = t_spar


def q_b(booms, a, moi):
    diff_q = np.zeros(len(booms))

    alpha = -1 / moi
    value = 0
    for i in range(len(booms)):
        delta_q = alpha * booms[i][1] * a[i]
        value += delta_q
        diff_q[i] = value

    base_q = np.zeros(len(booms) - 1)
    start = 0
    for i in range(len(diff_q) - 1):
        base_q[i] = diff_q[i] - start
        start = diff_q[i]

    return base_q


def q_s0(base_q, t, s):
    flow = np.multiply(base_q, s)
    flow = np.divide(flow, t)

    num = np.sum(flow)
    x = np.sum(np.add(flow, np.divide(s, t)))
    y = s[-1] / t[-1]

    return x, y, num


def area_tiangle(pos1, pos2, pos3):
    len1 = np.sqrt((pos1[0] - pos2[0]) ** 2 + (pos1[1] - pos2[1]) ** 2)
    len2 = np.sqrt((pos2[0] - pos3[0]) ** 2 + (pos2[1] - pos3[1]) ** 2)
    len3 = np.sqrt((pos1[0] - pos3[0]) ** 2 + (pos1[1] - pos3[1]) ** 2)

    s = (len1 + len2 + len3) / 3
    area = np.sqrt((s * (s - len1) * (s - len2) * (s - len3)))

    return area


def area_circle(distance):
    ratio = distance / (h_a * np.pi / 2)

    area = ratio * np.pi * (h_a / 2) ** 2 / 2

    return area


qdiff = q_b(B, y, Izz)

qdiffI = np.append(qdiff[len(qdiff) - n_spar + 1:], qdiff[:n_spar])
qdiffI[-1] = 0

qdiffII = qdiff[n_spar - 1:len(qdiff) - n_spar + 2]
qdiffII[-1] = 0

qsI = np.array([q_s0(qdiffI, t_I, s_I)])
qsII = np.array([q_s0(qdiffII, t_II, s_II)])

sys = np.array([[qsI[0][0], qsI[0][1]], [qsII[0][0], qsII[0][1]]])
sys_vec = np.array([[qsI[0][2]], [qsII[0][2]]])

values = np.linalg.solve(sys, sys_vec)

total_area = np.pi * (h_a / 2) ** 2 / 2 + h_a * (c_a - h_a / 2)

q_finalI = qdiffI + values[0][0]
q_finalII = qdiffII + values[1][0]

yI = np.append(y[len(y) - n_spar + 1:], y[:n_spar])

# print(y)
# print(values)
# print(len(yI))
# print(len(qdiffII))

def sum_moments_I(y, z, shear_flow, summation_point, s):

    for i in range(len(y)):
        area = area_circle(s[i])

        moment = 2*area*shear_flow[i]


    

