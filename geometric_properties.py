import numpy as np
import matplotlib.pyplot as plt

def stiffner_area(t, h, w):  # Ai
    return w * h - ((w - t) * (h - t))


def wingbox_idealize(h_aileron, cord_aileron, n_stiffner):
    # ONLY FOR TESTING, DELETE AFTER:
    h_aileron, cord_aileron, n_stiffner = 0.173, 0.484, 13

    contour_length = np.pi * h_aileron / 2 + 2 * np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)
    distance_stiffner = contour_length / n_stiffner
    pos_stiffner = np.array([])
    pos_stiffner = np.append(pos_stiffner, 0)
    for i in range(1, n_stiffner):
        pos_stiffner = np.append(pos_stiffner, pos_stiffner[-1] + distance_stiffner)
    pos_spar = np.array([np.pi * h_aileron / 4, contour_length - np.pi * h_aileron / 4])
    return pos_stiffner, pos_spar, contour_length


def single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin, distance, y_pos):  # Bi
    for i in [-1, 1]:
        B = B + t_skin * distance / 6 * (2 + y_pos[1] / y_pos[i])
    if stiffner:
        B = B + a_stiffner
    if spar:
        B = B + a_spar / 6
    return B


def boom_area(aileron, stiffner, a_spar):  # pass aileron, stiffner as tupel, see below
    t_skin, h_aileron, cord_aileron = aileron
    t_stiffner, h_stiffner, w_stiffner = stiffner

    n_discretize = 20  # integer, 13<n<600 [contour_length/w_stiffner]
    pos_stiffner, pos_spar, contour_length = wingbox_idealize(1, 1, 1)
    distance_boom = contour_length / n_discretize

    pos_boom = np.array([])
    pos_boom = np.append(pos_boom, 0)
    for i in range(1, n_discretize):
        pos_boom = np.append(pos_boom, pos_boom[-1] + distance_boom)

    stiffner = False
    spar = False

    a_stiffner = stiffner_area(t_stiffner, h_stiffner, w_stiffner)

    y_pos_booms = np.array(sorin_func(h, c, range[n_discretize]))
    boom_areas = np.array([])
    for i in range(n_discretize):
        for j in pos_stiffner:
            if pos_boom[i] - distance_boom / 2 < pos_stiffner[j] < pos_boom[i] + distance_boom / 2:
                stiffner = True

        for j in pos_spar:
            if pos_boom[i] - distance_boom / 2 < pos_spar[j] < pos_boom[i] + distance_boom / 2:
                spar = True
        boom_areas = np.append(boom_areas, single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin, distance_boom,
                                                            y_pos_booms[i - 1:i + 1]))
    return


def moi(B, pos_booms):
    I = 0
    for i in range(len(B)):
        I = I + B[i] * pos_booms[i] ** 2
    return I


def cg(B, pos_booms):
    Q = 0
    for i in range(len(B)):
        Q = Q + B[i] * pos_booms[i]
    return Q / sum(B)


def map_s_to_yz(s, h_aileron, cord_aileron):  # not generalized for different wing geometries
    if 0 <= s <= np.pi * h_aileron / 4:
        y = h_aileron / 2 * (1 - np.cos(s / h_aileron / 2))
        z = h_aileron / 2 * np.sin(s / h_aileron / 2)
    elif np.pi * h_aileron / 4 <= s <= np.pi * h_aileron / 4 + np.sqrt(
            (cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2):
        y = h_aileron / 2 + (s - np.pi * h_aileron / 4) * np.cos(
            np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
        z = h_aileron / 2 - (s - np.pi * h_aileron / 4) * np.sin(
            np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
    else:
        s = s - np.pi * h_aileron / 4 + np.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)
        return map_s_to_yz(s, h_aileron, cord_aileron)
    return y, z


def main():
    return


# main()

#STILL TESTING, SORINS PROGRAM DOES NOT WORK COMPLETLY SO IM ADJUSTING IT TOMORROW (map_s_to_yz)
yl=[]
zl=[]
for i in np.linspace(0.,np.pi * 0.173 / 4 + np.sqrt((0.484 - 0.173 / 2) ** 2 + (0.173 / 2) ** 2),100):
    y,z = map_s_to_yz(i, 0.173, 0.484)
    yl.append(y)
    zl.append(z)

plt.plot(yl, zl)
y,z=map_s_to_yz(np.pi * 0.173 / 4,  0.173, 0.484)
plt.scatter(y,z,c="red")
y,z=map_s_to_yz(np.pi * 0.173 / 4 + np.sqrt((0.484 - 0.173 / 2) ** 2 + (0.173 / 2) ** 2),  0.173, 0.484)
plt.scatter(y,z,c="green")
plt.axis('equal')
plt.show()

print(str(map_s_to_yz(np.pi * 0.173 / 4, 0.173, 0.484)))
