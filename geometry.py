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
        self.w_stiff = 1.8 / 100
        self.h_stiff = 1.4 / 100
        self.t_stiff = 1.2 / 1000
        self.t_spar = 2.5 / 1000
        self.stiff_area = self.n_stiff * (self.w_stiff + self.h_stiff) * self.t_stiff

    def _calculate_boom_distance(self):
        """
        Calculates the distances between all the booms if they are spaced equally around the section.

        :return: Float, containing the distance between all the booms
        """
        s = math.pi * self.h_a / 2 + 2 * math.sqrt((self.c_a - (self.h_a / 2)) ** 2 + (self.h_a / 2) ** 2)
        d_boom = s / self.n_boom

        # print("the boom distance equals : " + str(d_boom))
        return d_boom

    def set_boom_locations(self):
        """
        Setting the locations of the booms around the section by calculating the position based on the position along
        the circumference of the section

        :return: Two float lists. The first contains the locations of all the booms. The second contains their position
        along the circumference of the section.
        """
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

        # distances.append(distance * self.n_boom)
        return booms, distances

    def _map_s_to_yz(self, s):
        """
        Calculates the position of a boom based on its position along the circumference of the section.

        :param s: Float, the circumference distance of the boom along the cross-section.
        :return: List of floats, containing the positions in the z-y frame of the cross-section. Together with a number
        representing the quadrant in which the boom is located.
        """
        d_spar_top = math.pi * self.h_a / 4
        d_triangle = math.sqrt((self.c_a - (self.h_a / 2)) ** 2 + (self.h_a / 2) ** 2)
        alpha = np.arcsin((self.c_a - (self.h_a / 2)) / d_triangle)  # corner at top of spar
        beta = np.arcsin((self.h_a / 2) / d_triangle)  # corner at trailing edge

        if s < 0. or s > (2 * d_spar_top + 2 * d_triangle):
            print("S is outside bounds. s = " + str(s))
            return None

        # print("the spar top distance is " + str(d_spar_top))
        # print("the value of s is " + str(s))

        if 0. < s <= d_spar_top:
            # print("boom in first quadrant")
            angle = (s / d_spar_top) * math.pi / 2
            # print("the angle is " + str(np.degrees(angle)))

            z = self.h_a / 2 - (math.cos(angle) * self.h_a / 2)
            y = self.h_a / 2 * math.sin(angle)
            return z, y, 1

        elif d_spar_top < s <= (d_spar_top + d_triangle):
            # print("boom in second quadrant")
            s -= d_spar_top

            z = self.h_a / 2 + s * math.sin(alpha)
            y = self.h_a / 2 - s * math.cos(alpha)
            return z, y, 2

        elif (d_spar_top + d_triangle) < s <= (d_spar_top + 2 * d_triangle):
            s = s - d_spar_top - d_triangle

            z = self.c_a - s * math.cos(beta)
            y = 0 - s * math.sin(beta)
            return z, y, 3

        else:
            s = s - d_spar_top - 2 * d_triangle
            angle = (s / d_spar_top) * math.pi / 2

            z = self.h_a / 2 - self.h_a / 2 * math.sin(angle)
            y = - self.h_a / 2 + (self.h_a / 2 - math.cos(angle) * self.h_a / 2)
            return z, y, 4

    def add_spar_booms(self, booms, distances, boom_areas):
        """
        Adds two more booms that are on top of and at the bottom of the spar. Furthermore, all parameters affected by
        this are changed as well.

        :param booms: List, containing all the positions of the booms in the cross-section.
        :param distances: List, containing floats that represent the circumferential position of all the booms in the
        cross-section.
        :param boom_areas: List, containing floats that represent the area of all the booms in the cross-section.
        :return: Returns all the lists containing the booms at the spars as well.
        """
        distance = self._calculate_boom_distance()

        for i in range(len(booms)):
            if booms[i][0] > self.h_a / 2:
                n_top_spar = i
                # print("top spar is boom " + str(i + 1))
                break

        # print("the length of boom_areas is: " + str(len(boom_areas)))

        booms = booms[:n_top_spar] + [[self.h_a / 2, self.h_a / 2]] + booms[n_top_spar:]
        distances = distances[:n_top_spar] + [np.pi * self.h_a / 4] + distances[n_top_spar:]

        area_before_top_spar = ((self.t_skin * distance) / 6) * (
                2 + (booms[n_top_spar - 2][1] / booms[n_top_spar - 1][1])) + \
                               ((self.t_skin * (distances[n_top_spar] - distances[n_top_spar - 1])) / 6) * \
                               (2 + (booms[n_top_spar][1] / booms[n_top_spar - 1][1]))

        area_after_top_spar = ((self.t_skin * distance) / 6) * (
                2 + (booms[n_top_spar + 2][1] / booms[n_top_spar + 1][1])) + \
                              ((self.t_skin * (distances[n_top_spar + 1] - distances[n_top_spar])) / 6) * \
                              (2 + (booms[n_top_spar][1] / booms[n_top_spar + 1][1]))

        area_top_spar = ((self.t_skin * (distances[n_top_spar + 1] - distances[n_top_spar])) / 6) * \
                        (2 + (booms[n_top_spar + 1][1] / booms[n_top_spar][1])) + \
                        ((self.t_skin * (distances[n_top_spar] - distances[n_top_spar - 1])) / 6) * \
                        (2 + (booms[n_top_spar - 1][1] / booms[n_top_spar][1])) + \
                        ((self.t_spar * self.h_a) / 6)

        boom_areas = boom_areas[:n_top_spar - 1] + [area_before_top_spar] + [area_top_spar] + [
            area_after_top_spar] + boom_areas[n_top_spar + 1:]

        # print("the length of boom_areas is: " + str(len(boom_areas)))

        for i in range(len(booms)):
            if booms[i][1] < 0 and booms[i][0] < self.h_a / 2:
                n_bottom_spar = i
                # print("bottom spar at " + str(i + 1)
                break

        booms = booms[:n_bottom_spar] + [[self.h_a / 2, -self.h_a / 2]] + booms[n_bottom_spar:]
        distances = distances[:n_bottom_spar] + [
            ((math.pi * self.h_a) / 4) + 2 * (
                math.sqrt((self.h_a / 2) ** 2 + (self.c_a - (self.h_a / 2)) ** 2))] + distances[
                                                                                      n_bottom_spar:]

        area_before_top_spar = ((self.t_skin * distance) / 6) * (
                2 + (booms[n_bottom_spar - 2][1] / booms[n_bottom_spar - 1][1])) + \
                               ((self.t_skin * (distances[n_bottom_spar] - distances[n_bottom_spar - 1])) / 6) * \
                               (2 + (booms[n_bottom_spar][1] / booms[n_bottom_spar - 1][1]))

        area_after_top_spar = ((self.t_skin * distance) / 6) * (
                2 + (booms[n_bottom_spar + 2][1] / booms[n_bottom_spar + 1][1])) + \
                              ((self.t_skin * (distances[n_bottom_spar + 1] - distances[n_bottom_spar])) / 6) * \
                              (2 + (booms[n_bottom_spar][1] / booms[n_bottom_spar + 1][1]))

        area_top_spar = ((self.t_skin * (distances[n_bottom_spar + 1] - distances[n_bottom_spar])) / 6) * \
                        (2 + (booms[n_bottom_spar + 1][1] / booms[n_bottom_spar][1])) + (
                                (self.t_skin * (distances[n_bottom_spar] - distances[n_bottom_spar - 1])) / 6) * \
                        (2 + (booms[n_bottom_spar - 1][1] / booms[n_bottom_spar][1])) + \
                        ((self.t_spar * self.h_a) / 6)

        boom_areas = boom_areas[:n_bottom_spar - 1] + [area_before_top_spar] + [area_top_spar] + [
            area_after_top_spar] + boom_areas[n_bottom_spar + 1:]

        return booms, distances, boom_areas, n_top_spar, n_bottom_spar

    def calculate_boom_area(self, booms):
        """
        Calculates the area of all the booms.

        :param booms: List, containing all the boom positinos in the y-z coordinate-frame.
        :return: List, containing floats that represent the area of all the booms in the cross-section.
        """
        boom_areas = []
        distance = self._calculate_boom_distance()

        for i in range(len(booms)):
            before = ((self.t_skin * distance) / 6) * (2 + (booms[-1 + i][1] / booms[i][1]))
            # print("the ratio with the boom before: " + str(booms[-1 + i][1] / booms[i][1]))

            if i == len(booms) - 1:
                after = ((self.t_skin * distance) / 6) * (2 + (booms[0][1] / booms[i][1]))
                # print("the final ratio equals: " + str(booms[0][1] / booms[i][1]))

            else:
                after = ((self.t_skin * distance) / 6) * (2 + (booms[i + 1][1] / booms[i][1]))

            boom_areas.append(before + after + (self.stiff_area / self.n_boom))

        return boom_areas

    # def calculate_boom_area_z(self, booms, n_top, n_bottom):
    #     """
    #     Calculates the area of all the booms.
    #
    #     :param booms: List, containing all the boom positinos in the y-z coordinate-frame.
    #     :return: List, containing floats that represent the area of all the booms in the cross-section.
    #     """
    #     boom_areas = []
    #     distance = self._calculate_boom_distance()
    #
    #     for i in range(len(booms)):
    #         before = ((self.t_skin * distance) / 6) * (2 + (booms[-1 + i][0] / booms[i][0]))
    #         # print("the ratio with the boom before: " + str(booms[-1 + i][1] / booms[i][1]))
    #
    #         if i == len(booms) - 1:
    #             after = ((self.t_skin * distance) / 6) * (2 + (booms[0][0] / booms[i][0]))
    #             # print("the final ratio equals: " + str(booms[0][1] / booms[i][1]))
    #         else:
    #             after = ((self.t_skin * distance) / 6) * (2 + (booms[i + 1][0] / booms[i][0]))
    #
    #         if i == n_top or i == n_bottom:
    #             spar_addition = (self.t_spar * self.h_a) / 6
    #             boom_areas.append(before + after + spar_addition)
    #         else:
    #             boom_areas.append(before + after + (self.stiff_area / self.n_boom))
    #
    #     return boom_areas

    @staticmethod
    def calculate_cg(pos, area):
        """
        Calculates the position of the z-coordinate of the center of gravity of the cross-section.

        :param pos: List, containing the positions off all the booms in the cross-section.
        :param area: List, containing the area off all the booms in the cross-section.
        :return: Float, represents the z-coordinate of the center of gravity.
        """
        nom = 0
        denom = 0
        print(len(pos))
        print(len(area))

        for i in range(len(pos)):
            nom += area[i] * pos[i][0]
            denom += area[i]

        z_pos = nom / denom

        # print(pos)
        for i in range(len(pos)):
            pos[i][0] = -1 * pos[i][0] + z_pos

        # print(pos)
        return z_pos, pos

    @staticmethod
    def calculate_moi(pos, area):
        """
        Calculates the moments of inertia arond both the z- and y-axis.

        :param pos: List, containing the positions of all the booms in the cross-section.
        :param area: List, containing the area of all the booms in the cross-section.
        :return: Floats, representing the moments of inertia around the z- and y-axis.
        """
        izz = 0
        iyy = 0

        for i in range(len(pos)):
            izz += area[i] * (pos[i][1] ** 2)
            iyy += area[i] * (pos[i][0] ** 2)

        return izz, iyy


# Simple get functions for specific parameters if they are required.
def get_booms(n_booms):
    ideal = Idealization(n_booms)
    booms, distances = ideal.set_boom_locations()
    areas = ideal.calculate_boom_area(booms)
    booms, distances, areas_y, top, bottom = ideal.add_spar_booms(booms, distances, areas)
    # areas_z = ideal.calculate_boom_area_z(booms, top, bottom)

    return booms, areas, ideal


def get_cg(n_booms):
    booms, areas, ideal = get_booms(n_booms)
    z_pos, booms = ideal.calculate_cg(booms, areas)

    return z_pos, booms, ideal, areas


def get_moi(n_booms):
    z_pos, booms, ideal, areas = get_cg(n_booms)
    izz, iyy = ideal.calculate_moi(booms, areas)

    return izz, iyy


def get_boom_information(n_booms):
    ideal = Idealization(n_booms)
    booms, distances = ideal.set_boom_locations()
    areas = ideal.calculate_boom_area(booms)
    booms, distances, areas, top, bottom = ideal.add_spar_booms(booms, distances, areas)

    return distances, areas, top, bottom


# def plotting(booms):
#     for i in range(len(booms)):
#         plt.scatter(booms[i][0], booms[i][1], color= 'blue', label="Cross-section")
#
#     plt.show()
#
# plotting(b)

