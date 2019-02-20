import numpy as np
import matplotlib.pyplot as plt


def stiffner_area(t, h, w):  # Ai
    return w * h - ((w - t) * (h - t))


def wingbox_idealize(h_aileron, cord_aileron, n_stiffner):
    contour_length = np.pi * h_aileron / 2 + 2 * np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)
    distance_stiffner = contour_length / n_stiffner
    pos_stiffner = np.array([])
    pos_stiffner = np.append(pos_stiffner, 0)
    for i in range(1, n_stiffner):
        pos_stiffner = np.append(pos_stiffner, pos_stiffner[-1] + distance_stiffner)
    pos_spar = np.array([np.pi * h_aileron / 4, contour_length - np.pi * h_aileron / 4])
    return pos_stiffner, pos_spar, contour_length


def single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin, distance, y_pos):  # Bi
    B = 0
    # print(y_pos)
    for i in [0, 2]:
        if y_pos[i] != 0.:
            B = B + t_skin * distance / 6 * (2 + y_pos[1] / y_pos[i])
        else:
            B = B + t_skin * distance / 6 * (2 - 1)  # assume pure bending moment
    if stiffner:
        B = B + a_stiffner
    if spar:
        B = B + a_spar / 6
    return B


def area_units(stiffner, spar, a_stiffner, a_spar, t_skin, distance):
    a = t_skin * distance  # assume pure bending moment
    if stiffner:
        a = a + a_stiffner
    if spar:
        a = a + a_spar / 2
    return a


def boom_area(aileron, stiffner, a_spar, n_discretize):  # pass aileron, stiffner as tupel, see below
    t_skin, h_aileron, cord_aileron = aileron
    t_stiffner, h_stiffner, w_stiffner = stiffner

    pos_stiffner, pos_spar, contour_length = wingbox_idealize(h_aileron, cord_aileron, 13)
    distance_boom = contour_length / n_discretize

    pos_boom = np.array([])
    pos_boom = np.append(pos_boom, 0)
    for i in range(1, n_discretize):
        pos_boom = np.append(pos_boom, pos_boom[-1] + distance_boom)

    stiffner = False
    spar = False

    a_stiffner = stiffner_area(t_stiffner, h_stiffner, w_stiffner)

    z_pos_booms = np.array([])
    y_pos_booms = np.array([])
    for i in range(n_discretize):
        y, z = map_s_to_yz(i * distance_boom, h_aileron, cord_aileron)
        z_pos_booms = np.append(z_pos_booms, z)
        y_pos_booms = np.append(y_pos_booms, y)

    boom_areas_z = np.array([])
    areas = np.array([])
    for i in range(n_discretize):
        areas = np.append(areas, area_units(stiffner, spar, a_stiffner, a_spar, t_skin, distance_boom))
        stiffner, spar = False, False
        for j in pos_stiffner:
            if pos_boom[i] - distance_boom / 2 < j < pos_boom[i] + distance_boom / 2:
                stiffner = True

        for j in pos_spar:
            if pos_boom[i] - distance_boom / 2 < j < pos_boom[i] + distance_boom / 2:
                spar = True

        if i == n_discretize - 1:  # not sure if still necessary, for z=0 for elements on y axis would be nicer and functionally equivalent
            z_pos_local_booms = np.array([z_pos_booms[i - 1], z_pos_booms[i], z_pos_booms[0]])
        else:
            z_pos_local_booms = np.array([z_pos_booms[i - 1], z_pos_booms[i], z_pos_booms[i + 1]])
        boom_areas_z = np.append(boom_areas_z,
                                 single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin, distance_boom,
                                                  z_pos_local_booms))

    cg_y_coordinate = cg(areas, y_pos_booms)  # gives cg distance from LE
    boom_areas_y = np.array([])
    for i in range(n_discretize):
        stiffner, spar = False, False
        for j in pos_stiffner:
            if pos_boom[i] - distance_boom / 2 < j < pos_boom[i] + distance_boom / 2:
                stiffner = True

        for j in pos_spar:
            if pos_boom[i] - distance_boom / 2 < j < pos_boom[i] + distance_boom / 2:
                spar = True
        if i == n_discretize - 1:  # not sure if still necessary, for z=0 for elements on y axis would be nicer and functionally equivalent
            y_pos_local_booms = np.array([y_pos_booms[i - 1] - cg_y_coordinate, y_pos_booms[i] - cg_y_coordinate,
                                          y_pos_booms[0] - cg_y_coordinate])
        else:
            y_pos_local_booms = np.array([y_pos_booms[i - 1] - cg_y_coordinate, y_pos_booms[i] - cg_y_coordinate,
                                          y_pos_booms[i + 1] - cg_y_coordinate])
        boom_areas_y = np.append(boom_areas_y,
                                 single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin, distance_boom,
                                                  y_pos_local_booms))

    return boom_areas_y, y_pos_booms, z_pos_booms, boom_areas_z


def moi(B, pos_booms):
    I = 0
    for i in range(len(B)):
        I = I + B[i] * pos_booms[i] ** 2
    return I


def cg(A, pos_booms):
    Q = 0
    for i in range(len(A)):
        Q = Q + A[i] * pos_booms[i]
    return Q / sum(A)


def map_s_to_yz(s, h_aileron, cord_aileron):  # not generalized for different wing geometries
    lower_side = False
    if s > 2 * (np.pi * h_aileron / 4 + np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)) or s < 0.:
        print("S is outside bounds. s = " + str(s))
        print(s / np.pi * h_aileron / 2 + 2 * np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2))

    if s >= np.pi * h_aileron / 4 + np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2):
        diff = s - (np.pi * h_aileron / 4 + np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2))
        s = np.pi * h_aileron / 4 + np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2) - diff
        lower_side = True

    if 0. <= s <= np.pi * h_aileron / 4:
        y = h_aileron / 2 * (1 - np.cos(s / (h_aileron / 2)))
        z = h_aileron / 2 * np.sin(s / (h_aileron / 2))
        if lower_side:
            z = -z

    elif np.pi * h_aileron / 4 < s <= np.pi * h_aileron / 4 + np.sqrt(
            (cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2):
        y = h_aileron / 2 + (s - np.pi * h_aileron / 4) * np.cos(
            np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
        z = h_aileron / 2 - (s - np.pi * h_aileron / 4) * np.sin(
            np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
        if lower_side:
            z = -z
    return y, z


def main():
    aileron = (t_skin, h_aileron, cord_aileron) = (0.0011, 0.173, 0.484)
    stiffner = (t_stiffner, h_stiffner, w_stiffner) = (0.0012, 0.014, 0.018)
    a_spar = 0.0025 * 0.173
    n_discretize = 601  # only ODD integer, 13<n<600 [contour_length/w_stiffner]
    boom_areas_y, y_pos_booms, z_pos_booms, boom_areas_z = boom_area(aileron, stiffner, a_spar, n_discretize)
    # for i in range(len(boom_areas_z)):
    #     plt.scatter(y_pos_booms[i], z_pos_booms[i], s=60 * boom_areas_z[i] / max(boom_areas_z), color="blue",
    #                 label="around y axis")
    # for i in range(len(z_pos_booms)):
    #     z_pos_booms[i] = z_pos_booms[i] + h_aileron * 1.2
    # for i in range(len(boom_areas_y)):
    #     plt.scatter(y_pos_booms[i], z_pos_booms[i], s=60 * boom_areas_y[i] / max(boom_areas_y), color="red",
    #                 label="around z' axis")
    # plt.axis('equal')
    Izz = moi(boom_areas_z, z_pos_booms)
    Iyy = moi(boom_areas_y, y_pos_booms)
    print("Izz: " + str(Izz))
    print("Iyy: " + str(Iyy))
    plt.show()
    return

main()
