import numpy as np
import math
import matplotlib.pyplot as plt


class Idealization:

    def __init__(self, n_boom):
        self.h_a = 0.173
        self.c_a = 0.484
        self.t_skin = 1.1 / 1000
        self.n_boom = n_boom
        self.n_stiff = 13
        self.w_stiff = 1.8 / 1000
        self.h_stiff = 1.4 / 1000
        self.t_stiff = 1.2 / 1000
        self.t_spar = 2.5 / 1000
        self.stiff_area = self.n_stiff * (self.w_stiff * self.h_stiff - ((self.w_stiff - self.t_stiff) *
                                                                         (self.h_stiff - self.t_stiff)))

    def _calculate_boom_distance(self):
        s = np.pi * self.h_a / 2 + 2 * np.sqrt((self.c_a - (self.h_a / 2)) ** 2 + (self.h_a / 2) ** 2)
        d_boom = s / self.n_boom
        # print("the boom distance equals : " + str(d_boom))
        return d_boom

    def set_boom_locations(self):
        distance = self._calculate_boom_distance()
        s = distance / 2
        booms = []
        distances = []

        for i in range(self.n_boom):
            z, y, quadrant = self._map_s_to_yz(s)
            # print("the coordinates are (z, y) : " + str(z), " " + str(y), " " + str(quadrant))
            booms.append([z, y])
            distances.append(s)
            s += distance

        distances.append(distance * self.n_boom)
        return booms, distances

    def _map_s_to_yz(self, s):
        d_spar_top = math.pi * self.h_a / 4
        d_triangle = math.sqrt((self.c_a - (self.h_a / 2)) ** 2 + (self.h_a / 2) ** 2)
        alpha = np.arcsin((self.c_a - (self.h_a / 2)) / d_triangle)  # corner at top of spar
        beta = np.arcsin((self.h_a / 2) / d_triangle)  # corner at trailing edge

        if s < 0. or s > 2 * d_spar_top + 2 * d_triangle:
            print("S is outside bounds. s = " + str(s))
            return None
        # print("the spar top distance is " + str(d_spar_top))
        # print("the value of s is " + str(s))

        if 0. < s <= d_spar_top:
            # print("boom in first quadrant")
            angle = (s / (math.pi * self.h_a / 4)) * math.pi / 2

            # print("the angle is " + str(np.degrees(angle)))

            z = self.h_a / 2 - (np.cos(angle) * self.h_a / 2)
            y = self.h_a / 2 * np.sin(angle)
            return z, y, 1

        elif d_spar_top < s <= d_spar_top + d_triangle:
            # print("boom in second quadrant")
            s -= d_spar_top

            z = self.h_a / 2 + s * np.sin(alpha)
            y = self.h_a / 2 - s * np.cos(alpha)
            return z, y, 2

        elif d_spar_top + d_triangle < s <= d_spar_top + 2 * d_triangle:
            s = s - d_spar_top - d_triangle

            z = self.c_a - s * np.cos(beta)
            y = 0 - s * np.sin(beta)
            return z, y, 3

        else:
            s = s - d_spar_top - 2 * d_triangle
            angle = (s / (math.pi * self.h_a / 4)) * math.pi / 2

            z = self.h_a / 2 - self.h_a / 2 * np.sin(angle)
            y = - self.h_a / 2 + (self.h_a / 2 - np.cos(angle) * self.h_a / 2)
            return z, y, 4

    def add_spar_booms(self, booms, distances, boom_areas):
        distance = self._calculate_boom_distance()

        for i in range(len(booms)):
            if booms[i][0] > self.h_a / 2:
                n_top_spar = i
                print("top spar at " + str(i + 1))
                break

        booms = booms[:n_top_spar] + [[self.h_a / 2, self.h_a / 2]] + booms[n_top_spar:]
        distances = distances[:n_top_spar] + [np.pi * self.h_a / 4] + distances[n_top_spar:]

        area_before_top_spar = self.t_skin * distance / 6 * (2 + booms[n_top_spar - 2][1] / booms[n_top_spar - 1][1]) + \
                               self.t_skin * (distances[n_top_spar] - distances[n_top_spar - 1]) / 6 * \
                               (2 + booms[n_top_spar][1] / booms[n_top_spar - 1][1]) + self.stiff_area / (
                                       self.n_boom + 2)

        area_after_top_spar = self.t_skin * distance / 6 * (2 + booms[n_top_spar + 2][1] / booms[n_top_spar + 1][1]) + \
                              self.t_skin * (distances[n_top_spar + 1] - distances[n_top_spar]) / 6 * \
                              (2 + booms[n_top_spar][1] / booms[n_top_spar + 1][1]) + self.stiff_area / (
                                      self.n_boom + 2)

        area_top_spar = self.t_skin * (distances[n_top_spar + 1] - distances[n_top_spar]) / 6 * \
                        (2 + booms[n_top_spar + 1][1] / booms[n_top_spar][1]) + self.stiff_area / (self.n_boom + 2) + \
                        self.t_skin * (distances[n_top_spar] - distances[n_top_spar - 1]) / 6 * \
                        (2 + booms[n_top_spar - 1][1] / booms[n_top_spar][1]) + \
                        self.t_spar * self.h_a / 6 + self.stiff_area / (self.n_boom + 2)

        boom_areas = boom_areas[:n_top_spar - 1] + [area_before_top_spar] + [area_top_spar] + [
            area_after_top_spar] + boom_areas[n_top_spar + 1:]

        for i in range(len(booms)):
            if booms[i][1] < 0 and booms[i][0] < self.h_a / 2:
                n_bottom_spar = i
                print("bottom spar at " + str(i + 1))
                break

        booms = booms[:n_bottom_spar] + [[self.h_a / 2, -self.h_a / 2]] + booms[n_bottom_spar:]
        distances = distances[:n_bottom_spar] + [
            np.pi * self.h_a / 4 + 2 * np.sqrt((self.h_a / 2) ** 2 + (self.c_a - (self.h_a / 2)) ** 2)] + distances[
                                                                                                          n_bottom_spar:]

        area_before_top_spar = self.t_skin * distance / 6 * (
                2 + booms[n_bottom_spar - 2][1] / booms[n_bottom_spar - 1][1]) + \
                               self.t_skin * (distances[n_bottom_spar] - distances[n_bottom_spar - 1]) / 6 * \
                               (2 + booms[n_bottom_spar][1] / booms[n_bottom_spar - 1][1]) + self.stiff_area / (
                                       self.n_boom + 2)

        area_after_top_spar = self.t_skin * distance / 6 * (
                2 + booms[n_bottom_spar + 2][1] / booms[n_bottom_spar + 1][1]) + \
                              self.t_skin * (distances[n_bottom_spar + 1] - distances[n_bottom_spar]) / 6 * \
                              (2 + booms[n_bottom_spar][1] / booms[n_bottom_spar + 1][1]) + self.stiff_area / (
                                      self.n_boom + 2)

        area_top_spar = self.t_skin * (distances[n_bottom_spar + 1] - distances[n_bottom_spar]) / 6 * \
                        (2 + booms[n_bottom_spar + 1][1] / booms[n_bottom_spar][1]) + self.stiff_area / (
                                self.n_boom + 2) + \
                        self.t_skin * (distances[n_bottom_spar] - distances[n_bottom_spar - 1]) / 6 * \
                        (2 + booms[n_bottom_spar - 1][1] / booms[n_bottom_spar][1]) + \
                        self.t_spar * self.h_a / 6 + self.stiff_area / (self.n_boom + 2)

        boom_areas = boom_areas[:n_bottom_spar - 1] + [area_before_top_spar] + [area_top_spar] + [
            area_after_top_spar] + boom_areas[n_bottom_spar + 1:]

        return booms, distances, boom_areas, n_top_spar, n_bottom_spar

    def calculate_boom_area(self, booms):
        boom_areas = []
        distance = self._calculate_boom_distance()

        for i in range(len(booms)):
            before = self.t_skin * distance / 6 * (2 + booms[-1 + i][1] / booms[i][1])
            if i == len(booms) - 1:
                after = self.t_skin * distance / 6 * (2 + booms[0][1] / booms[i][1])
            else:
                after = self.t_skin * distance / 6 * (2 + booms[i + 1][1] / booms[i][1])
            boom_areas.append(before + after + self.stiff_area / (self.n_boom + 2))

        return boom_areas

    @staticmethod
    def calculate_cg(pos, a):
        nom = 0
        denom = 0

        for i in range(len(pos)):
            nom += a[i] * pos[i][0]
            denom += a[i]

        z_pos = nom / denom

        for i in range(len(pos)):
            pos[i][0] - z_pos

        return z_pos, pos

    @staticmethod
    def calculate_moi(pos, a):
        izz = 0
        iyy = 0

        for i in range(len(pos)):
            izz += a[i] * pos[i][1] ** 2
            iyy += a[i] * pos[i][0] ** 2

        return izz, iyy


ideal = Idealization(20)

b, d = ideal.set_boom_locations()

areas = ideal.calculate_boom_area(b)

b, d, areas, top, bottom = ideal.add_spar_booms(b, d, areas)

z_pos, b = ideal.calculate_cg(b, areas)

i_zz, i_yy = ideal.calculate_moi(b, areas)

# def plotting(booms):
#     for i in range(len(booms)):
#         plt.scatter(booms[i][0], booms[i][1], color= 'blue', label="Cross-section")
#
#     plt.show()
#
# plotting(b)

print(top)
print(bottom)
