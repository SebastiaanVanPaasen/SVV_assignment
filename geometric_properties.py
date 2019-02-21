class IdealizedStructure():  # does NOT create instances, methods are static and values are class values
    # <editor-fold desc="CLASS METHODS">
    # the methods are all static class methods and should not be used on instances as none will be initialized by the class
    @staticmethod
    def stiffner_area(t, h, w):  # Ai
        return w * h - ((w - t) * (h - t))

    @staticmethod
    def wingbox_idealize(h_aileron, cord_aileron, n_stiffner):
        contour_length = math.pi * h_aileron / 2 + 2 * math.sqrt(
            (cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)
        distance_stiffner = contour_length / n_stiffner
        pos_stiffner = np.array([])
        pos_stiffner = np.append(pos_stiffner, 0)
        for i in range(1, n_stiffner):
            pos_stiffner = np.append(pos_stiffner, pos_stiffner[-1] + distance_stiffner)
        pos_spar = np.array([math.pi * h_aileron / 4, contour_length - math.pi * h_aileron / 4])
        return pos_stiffner, pos_spar, contour_length

    @staticmethod
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

    @staticmethod
    def area_units(stiffner, spar, a_stiffner, a_spar, t_skin, distance):
        a = t_skin * distance  # assume pure bending moment
        if stiffner:
            a = a + a_stiffner
        if spar:
            a = a + a_spar / 2
        return a

    @staticmethod
    def boom_area(aileron, stiffner, a_spar, n_discretize):  # pass aileron, stiffner as tupel, see below
        t_skin, h_aileron, cord_aileron = aileron
        t_stiffner, h_stiffner, w_stiffner = stiffner

        pos_stiffner, pos_spar, contour_length = IdealizedStructure.wingbox_idealize(h_aileron, cord_aileron, 13)
        distance_boom = contour_length / n_discretize

        pos_boom = np.array([])
        pos_boom = np.append(pos_boom, 0)
        for i in range(1, n_discretize):
            pos_boom = np.append(pos_boom, pos_boom[-1] + distance_boom)

        stiffner = False
        spar = False

        a_stiffner = IdealizedStructure.stiffner_area(t_stiffner, h_stiffner, w_stiffner)

        z_pos_booms = np.array([])
        y_pos_booms = np.array([])
        for i in range(n_discretize):
            y, z = IdealizedStructure.map_s_to_yz(i * distance_boom, h_aileron, cord_aileron)
            z_pos_booms = np.append(z_pos_booms, z)
            y_pos_booms = np.append(y_pos_booms, y)

        boom_areas_z = np.array([])
        areas = np.array([])
        for i in range(n_discretize):
            areas = np.append(areas,
                              IdealizedStructure.area_units(stiffner, spar, a_stiffner, a_spar, t_skin, distance_boom))
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
                                     IdealizedStructure.single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin,
                                                                         distance_boom,
                                                                         z_pos_local_booms))

        cg_y_coordinate = IdealizedStructure.cg(areas, y_pos_booms)  # gives cg distance from LE
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
                                     IdealizedStructure.single_boom_area(stiffner, spar, a_stiffner, a_spar, t_skin,
                                                                         distance_boom,
                                                                         y_pos_local_booms))

        return boom_areas_y, y_pos_booms, z_pos_booms, boom_areas_z

    @staticmethod
    def moi(B, pos_booms):
        I = 0
        for i in range(len(B)):
            I = I + B[i] * pos_booms[i] ** 2
        return I

    @staticmethod
    def cg(A, pos_booms):
        Q = 0
        for i in range(len(A)):
            Q = Q + A[i] * pos_booms[i]
        return Q / sum(A)

    @staticmethod
    def map_s_to_yz(s, h_aileron, cord_aileron):  # not generalized for different wing geometries
        lower_side = False
        if s > 2 * (
                math.pi * h_aileron / 4 + math.sqrt(
            (cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2)) or s < 0.:
            print("S is outside bounds. s = " + str(s))
            print(
                s / math.pi * h_aileron / 2 + 2 * math.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2))

        if s >= math.pi * h_aileron / 4 + math.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2):
            diff = s - (math.pi * h_aileron / 4 + math.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2))
            s = math.pi * h_aileron / 4 + math.sqrt((cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2) - diff
            lower_side = True

        if 0. <= s <= math.pi * h_aileron / 4:
            y = h_aileron / 2 * (1 - np.cos(s / (h_aileron / 2)))
            z = h_aileron / 2 * np.sin(s / (h_aileron / 2))
            if lower_side:
                z = -z

        elif math.pi * h_aileron / 4 < s <= math.pi * h_aileron / 4 + math.sqrt(
                (cord_aileron - h_aileron / 2) ** 2 + (h_aileron / 2) ** 2):
            y = h_aileron / 2 + (s - math.pi * h_aileron / 4) * np.cos(
                np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
            z = h_aileron / 2 - (s - math.pi * h_aileron / 4) * np.sin(
                np.arctan(h_aileron / 2 / (cord_aileron - h_aileron / 2)))
            if lower_side:
                z = -z
        return y, z

    # </editor-fold>

    # <editor-fold desc="CLASS VARIABLES">
    n_discretize = SimulationData.IdealizedStructure_n_discretize
    aileron = (ProblemData.t_skin, ProblemData.h_aileron, ProblemData.cord_aileron)
    stiffner = (ProblemData.t_stiffner, ProblemData.h_stiffner, ProblemData.w_stiffner)
    a_spar = ProblemData.t_spar * ProblemData.h_spar
    boom_areas_y, y_pos_booms, z_pos_booms, boom_areas_z = boom_area(aileron, stiffner, a_spar,
                                                                     n_discretize)
    izz = moi(boom_areas_z, z_pos_booms)
    iyy = moi(boom_areas_y, y_pos_booms)

    # </editor-fold>
