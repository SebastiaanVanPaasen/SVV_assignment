import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, SmoothBivariateSpline
import internal_forces
from heapq import nlargest
import bottleneck as bn
import math
from mpl_toolkits.mplot3d import Axes3D

# note: TODO: !!!=================================================================================!!!
# TODO: requires specific settings in internal_forces.py for interpolation, data restructure to function properly.
# TODO: if figures look differently from report you are using wrong input.
# TODO: set internal_forces.d_x = 0.01
# TODO: in internal_forces.calc_reaction_forces() use y,z which are commented out
# TODO: verify that you are using correct coordinate frame by checking dy_le_g2, dz_le_g2 (not global!)
# note: TODO: !!!=================================================================================!!!

# <editor-fold desc="GET METHODS">
def read_table_vali(filename,
                    replace_comma_with_dot):
    # tables must be in form: 1 row header (description), separation by tabulators
    file = open(filename, "r+")
    file = file.readlines()
    n_colums = len(file[0].split('\t'))
    del file[0]
    if replace_comma_with_dot:
        for i in range(len(file)):
            file[i] = file[i].replace(",", ".")
    n_rows = len(file)
    data = np.zeros(shape=(n_colums, n_rows))
    for i in range(n_rows):
        l = file[i].split('\t')
        for j in range(n_colums):
            try:
                data[j, i] = l[j]
            except ValueError:
                print("File: " + str(filename) + " | Line: " + str(i) + " | Value: " + str(l[j]))
    return data


def get_deformation_val(case):
    # output: returns an array with deformation_data[node label][magnitude def][x def][y def][z def] and
    # output: deformation[node label][x total][y total][z total]
    filename = str(case) + ".txt"
    nodes = read_table_vali("inp.txt", False)
    deformation_data = read_table_vali(filename, True)
    deformation_data = np.delete(deformation_data, [1, 2, 3, 4], 0)

    x = nodes[1] + deformation_data[2]
    y = nodes[2] + deformation_data[3]
    z = nodes[3] + deformation_data[4]

    # note: one of the tries how to make sense of the warped cord
    deformation = np.array([nodes[0], x, y, z])  # add node[0] subarray at [0] for convienience (plotting)
    return deformation, deformation_data


def get_stress_val(
        case):  # output: returns a sorted, unique array with stress[element label][x][y][z][avg_von_mieses_stress]
    filename = str(case) + ".txt"
    stress_data = read_table_vali(filename, True)
    elements = get_elements('elements')
    index_sorted = np.argsort(stress_data[0], axis=0, kind="quicksort")
    stress_data = np.array([stress_data[0][index_sorted], stress_data[1][index_sorted], stress_data[2][index_sorted],
                            stress_data[3][index_sorted]])
    unique, index_unique = np.unique(stress_data[0], axis=0, return_index=True)

    avg_stress_1 = np.array([])
    avg_stress_2 = np.array([])
    ctotal = 0
    for i in range(len(unique)):
        avg_stress_at_element_1 = 0
        avg_stress_at_element_2 = 0
        c = 0

        for j in range(ctotal, len(stress_data[0])):
            if stress_data[0, j] != unique[i]:
                break

            avg_stress_at_element_1 = avg_stress_at_element_1 + stress_data[2, j]
            avg_stress_at_element_2 = avg_stress_at_element_2 + stress_data[3, j]

            c += 1
            ctotal += 1

        avg_stress_1 = np.append(avg_stress_1, avg_stress_at_element_1 / c)
        avg_stress_2 = np.append(avg_stress_2, avg_stress_at_element_2 / c)

    avg_von_mieses = 0.5 * (avg_stress_1 + avg_stress_2)
    stress = np.array([elements[0], elements[1], elements[2], elements[3], avg_von_mieses])

    return stress


def get_elements(case):
    filename = str(case) + ".txt"
    nodes = read_table_vali("inp.txt", False)
    elements = read_table_vali(filename, True)

    element_pos = np.array([elements[0], [], [], []])
    for i in range(len(elements[0])):
        element_nodes = np.array([])
        for j in range(1, 5):
            current_node = elements[j, i] - 1
            current_node_pos = np.array(
                [nodes[[1], int(current_node)], nodes[[2], int(current_node)], nodes[[3], int(current_node)]])
            element_nodes = np.append(element_nodes, np.array([current_node_pos]))
        element_nodes = element_nodes.reshape((4, 3))
        x = (element_nodes[0, 0] + element_nodes[1, 0] + element_nodes[2, 0] + element_nodes[3, 0]) / 4
        y = (element_nodes[0, 1] + element_nodes[1, 1] + element_nodes[2, 1] + element_nodes[3, 1]) / 4
        z = (element_nodes[0, 2] + element_nodes[1, 2] + element_nodes[2, 2] + element_nodes[3, 2]) / 4
        element_pos[1] = np.append(element_pos[1], x)
        element_pos[2] = np.append(element_pos[2], y)
        element_pos[3] = np.append(element_pos[3], z)
    return element_pos


def get_nodes():
    nodes = read_table_vali("inp.txt", False)
    return nodes


def get_le_te_index(b_nodes, b_elments, case):  # takes: boolean
    nodes = get_nodes()
    index_sorted_nodes = np.argsort(nodes[3])

    index_LE_nodes = np.flip(index_sorted_nodes, axis=0)[..., 0:68].astype(int)
    index_TE_nodes = index_sorted_nodes[..., 0:68].astype(int)

    elements = get_elements(case)
    index_sorted_elements = np.argsort(nodes[3])

    index_LE_elements = np.flip(index_sorted_elements, axis=0)[..., 0:68].astype(int)
    index_TE_elements = index_sorted_elements[..., 0:68].astype(int)

    if b_nodes and b_elments:
        return index_LE_nodes, index_TE_nodes, index_LE_elements, index_TE_elements

    return index_LE_nodes, index_TE_nodes


def get_stress_sim():
    stress_sim, boom_locations, x_slice = internal_forces.v_m, internal_forces.boom_locations, internal_forces.x_slice
    stress_sim, boom_locations, x_slice = np.asarray(stress_sim), np.asarray(boom_locations), np.asarray(x_slice)
    stress_sim = np.ndarray.flatten(stress_sim)
    s = np.array([])
    for i in range(len(stress_sim)):
        for j in stress_sim[i]:
            s = np.append(s, j)
    stress_sim = np.reshape(s, (len(x_slice), 30))

    x = np.array([])
    y = np.array([])
    z = np.array([])
    s = np.array([])
    dx = internal_forces.d_x / 2
    get_spar_loc = True
    for i in range(len(x_slice)):
        for j in range(len(boom_locations)):
            x = np.append(x, x_slice[i])
            y = np.append(y, boom_locations[j, 0])
            z = np.append(z, boom_locations[j, 1])  #
            if 0.173 / 2 - dx < abs(y[-1]) < 0.173 / 2 + dx and z[-1] > 0.:
                x = np.append(x, x_slice[i])
                y = np.append(y, 0.)
                z = np.append(z, 0.173 / 2)
                if get_spar_loc:
                    get_spar_loc = False
                    y_spar, z_spar = y[-1], z[-1]
            elif 0.173 / 2 - dx < abs(y[-1]) < 0.173 / 2 + dx and z[-1] < 0.:
                x = np.append(x, x_slice[i])
                y = np.append(y, 0.)
                z = np.append(z, -0.173 / 2)
        for k in range(len(stress_sim[i])):
            s = np.append(s, stress_sim[i, k])
    stress = np.array([range(len(x)), x, z, y,
                       s])  # append pseudo node array for plotting, note see also here for coordinate transform
    delete_wired_beam = np.array([])
    for i in range(len(stress[0])):
        if 0.075 < abs(stress[2, i]) < 0.1 and -0.01 < stress[3, i] < 0.01:
            delete_wired_beam = np.append(delete_wired_beam, i)
    stress = np.delete(stress, delete_wired_beam, axis=1)
    stress[3] = stress[3] - 1.5 * z_spar
    return stress


def get_deformation_sim():
    deformation_le = np.array(
        [internal_forces.x_slice, np.zeros((1, len(internal_forces.x_slice))), internal_forces.dy_le_g2,
         internal_forces.dz_le_g2])
    deformation_te = np.array(
        [internal_forces.x_slice, np.zeros((1, len(internal_forces.x_slice))), internal_forces.dy_te_g2,
         internal_forces.dz_te_g2])
    return deformation_le, deformation_te


# </editor-fold>


# <editor-fold desc="PLOT METHODS">
def plot_node_pos(nodes, LE, TE):
    color = []
    for i in nodes[0]:
        color.append((0., i / max(nodes[0]), 0.))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(nodes[1], nodes[2], nodes[3], zdir="z", s=2, c="grey")
    ax.scatter(LE[1], LE[2], LE[3], zdir="z", s=10, c="darkcyan")
    ax.scatter(TE[1], TE[2], TE[3], zdir="z", s=10, c="blue")

    max_range = np.array(
        [nodes[1].max() - nodes[1].min(), nodes[2].max() - nodes[2].min(), nodes[3].max() - nodes[3].min()]).max() / 2.0

    mid_x = (nodes[1].max() + nodes[1].min()) * 0.5
    mid_y = (nodes[2].max() + nodes[2].min()) * 0.5
    mid_z = (nodes[3].max() + nodes[3].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    return


def plot_deformation_alternative(def_val, def_sim,
                                 title):
    # use case variable like "UR1", only will work with U...

    color1 = []
    for i in np.abs(def_val[-1]):
        color1.append((1 - 0.8 * (abs(i) / np.max(np.abs(def_val[-1]))), 0., 0.))  # high deformation -> dark color

    color2 = []
    for i in np.abs(def_sim[-1]):
        color2.append((0., 0., 1 - (abs(i) / np.max(np.abs(def_sim[-1])))))

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # note: def_sim[xtotal][x def][y def][z def][def mag]
    # note: def_val[xtotal][0][y def][z def][def mag]

    ax.scatter(def_val[0], def_val[2], def_val[3], zdir="z", s=10, c=color1, label="Validation data")
    ax.scatter(def_sim[0] + def_sim[1], def_sim[2], def_sim[3], zdir="z", s=10, c=color2, label="Model data")

    max_range = np.array(
        [def_val[0].max() - def_val[0].min(), def_val[2].max() - def_val[2].min(),
         def_val[3].max() - def_val[3].min()]).max() / 2.0

    mid_x = (def_val[0].max() + def_val[0].min()) * 0.5
    mid_y = (def_val[2].max() + def_val[2].min()) * 0.5
    mid_z = (def_val[3].max() + def_val[3].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X - span [m]')
    ax.set_ylabel('Z - cord [m]')
    ax.set_zlabel('Y [m]')
    ax.set_title(title)
    ax.legend(loc='center right')
    plt.show()
    return


def plot_deformation(deformation, deformation_data, title, c):
    # use case variable like "UR1", only will work with U
    if c[0]:
        color = []
        for i in np.abs(deformation_data[1]):
            color.append((1 - (abs(i) / np.max(np.abs(deformation_data[1]))), 0., 0.))
    else:
        choice = ["darkcyan", "blue", "green", "darkcyan", "orange"]
        color = choice[c[1]]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(deformation[1], deformation[2], deformation[3], zdir="z", s=10, c=color)

    max_range = np.array(
        [deformation[1].max() - deformation[1].min(), deformation[2].max() - deformation[2].min(),
         deformation[3].max() - deformation[3].min()]).max() / 2.0

    mid_x = (deformation[1].max() + deformation[1].min()) * 0.5
    mid_y = (deformation[2].max() + deformation[2].min()) * 0.5
    mid_z = (deformation[3].max() + deformation[3].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    plt.show()
    return


def plot_alt_deformation(deformation1, deformation2, title):
    # use case variable like "UR1", only will work with U...

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(deformation1[1], deformation1[2], deformation1[3], zdir="z", s=5, c="black", label="undeformed")
    ax.scatter(deformation2[1], deformation2[2], deformation2[3], zdir="z", s=5, c="red", label="deformed")

    max_range = np.array(
        [deformation1[1].max() - deformation1[1].min(), deformation1[2].max() - deformation1[2].min(),
         deformation1[3].max() - deformation1[3].min()]).max() / 2.0

    mid_x = (deformation1[1].max() + deformation1[1].min()) * 0.5
    mid_y = (deformation1[2].max() + deformation1[2].min()) * 0.5
    mid_z = (deformation1[3].max() + deformation1[3].min()) * 0.5

    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.legend(loc="lower right")
    plt.show()
    return


def plot_von_mieses_stress(stress, title):
    # stress = stress/1000. # mm -> m
    color = []
    maxplot= np.max(nlargest(100,np.abs(stress[4]),key=None))
    midplot = np.min(nlargest(1000, np.abs(stress[4]), key=None))
    minplot = np.min(nlargest(2000, np.abs(stress[4]), key=None))

    for i in np.abs(stress[4]):
        if i>maxplot:
            color.append((0., 0.,0.))# max black
        elif i>midplot:
            color.append((1 - 0.95 * (abs(i) / maxplot), 0.3, 0.3))# large red
        elif i>minplot:
            color.append((0.3, 1 - 0.95 * (abs(i) / midplot), 0.3))# mid green
        else:
            color.append((0.3, 0.3, 1 - 0.95 * (abs(i) / midplot)))# min blue

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(stress[1], stress[2], stress[3], zdir="z", s=2, c=color)

    max_range = np.array(
        [stress[1].max() - stress[1].min(), stress[2].max() - stress[2].min(),
         stress[3].max() - stress[3].min()]).max() / 2.0

    mid_x = (stress[1].max() + stress[1].min()) * 0.5
    mid_y = (stress[2].max() + stress[2].min()) * 0.5
    mid_z = (stress[3].max() + stress[3].min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    plt.show()
    return


# </editor-fold>


# <editor-fold desc="VALIDATION METHODS">
def reshape_data_grid(data):  # must be 4D array, discards 1st subarray
    flat_data = np.array([])
    for i in range(len(data[0])):
        flat_data = np.append(flat_data, np.array([data[1, i], data[2, i], data[3, i]]))
    data = np.reshape(flat_data, (len(data[0]), 3), order="A")
    return data


def alt_reshape_data_grid(data):  # must be 4D array, discards 1st subarray
    flat_data = np.array([])
    for i in range(len(data[0])):
        flat_data = np.append(flat_data, np.array([data[1][0][i], data[2][i], data[3][i]]))
    data = np.reshape(flat_data, (len(data[0]), 3), order="A")
    return data


def reverse_reshape_data_grid(data):  # must be ND array, adds 1st subarray in normal form
    flat_data = np.array([])
    for i in range(3):
        for j in range(len(data)):
            flat_data = np.append(flat_data, data[j, i])
    data = np.reshape(flat_data, (3, len(data)), order="A")
    data = np.array([range(len(data[0])), data[0], data[1], data[2]])
    return data


# TODO: wonder who wrote the FEA, like, honestly guys!

def coordiant_transform(x, y, z):
    x = x + 0
    y = y + 0
    z = z + 0
    return x, y, z  # return with switched coordinates as needed


def results_stress():
    stress_sim = get_stress_sim()
    # plot_von_mieses_stress(stress_sim)
    for i in ["SLC1"]:  # "SR1", "SR2", "SLC1", "SLC2"
        stress_vali = get_stress_val(
            i)  # put into Pa # output: returns a sorted, unique array with stress[element label][x][y][z][avg_von_mieses_stress]
        stress_vali[4] = stress_vali[4] * 1000 ** 2
        delete_wired_beam = np.array([])
        for i in range(len(stress_vali[0])):
            if -100 < stress_vali[2, i] < 100 and -0.9 * 0.173 / 2 < stress_vali[3, i] < 0.9 * 0.173 / 2:
                delete_wired_beam = np.append(delete_wired_beam, i)
        stress_vali = np.delete(stress_vali, delete_wired_beam, axis=1)

        sim_points = reshape_data_grid(stress_sim[0:4])
        vali_points = reshape_data_grid(stress_vali[0:4])

        stress_sim_at_vali = griddata(sim_points, stress_sim[4], vali_points, method="nearest")
        # note higher interpol methods than -narest- are not available for values over an unregular 3D-Spline

        MPE = 100. / len(stress_sim_at_vali) * (sum(stress_vali[4]) - sum(stress_sim_at_vali)) / sum(stress_vali[4])
        print("von Mises stress MPE over whole surface: " + str(
            int(10000. * MPE)/10.) + " %")  # note UNCOMMENT FOR RESULTS

        local_PE = (stress_vali[4] - stress_sim_at_vali) / stress_vali[4]
        local_pos = reverse_reshape_data_grid(vali_points)
        stress_PE = np.array([local_pos[0], local_pos[1], local_pos[2], local_pos[3], local_PE])

        select, s = 1, 3200
        plot_von_mieses_stress(stress_sim[..., range(0, len(stress_PE[0]), select)], 'Model von Mises stress')
        plot_von_mieses_stress(stress_PE[..., range(0, len(stress_PE[0]), 1)],
                               'von Mises stress difference function')  # note UNCOMMENT FOR RESULTS

        local_stress = stress_vali
        index_sorted = np.argsort(stress_vali[3])

        index_le = np.flip(index_sorted, axis=0)[..., 0:68].astype(int)
        local_stress_le = local_stress[..., index_le]

        index_te = index_sorted[..., 0:68].astype(int)
        local_stress_te = local_stress[..., index_te]

        index_x_sorted_le = np.argsort(local_stress_le[1])
        local_stress_le = local_stress_le[..., index_x_sorted_le]

        index_x_sorted_te = np.argsort(local_stress_te[1])
        local_stress_te = local_stress_te[..., index_x_sorted_te]

        fig, ax = plt.subplots(nrows=4, ncols=1, sharex=True, sharey=True)

        ax[0].plot(stress_sim[1, range(0, 28 * len(internal_forces.x_slice), 28)],
                   stress_sim[-1, range(0, 28 * len(internal_forces.x_slice), 28)] / (1000 ** 2) / s,
                   label='Model LE', marker='+', c='darkred', markevery=20,
                   linestyle='--', linewidth=1)
        ax[0].plot(local_stress_le[1] / 1000, local_stress_le[-1] / (1000 ** 2), label='Validation LE', marker='^',
                   c='darkblue',
                   markevery=20, linewidth=1)
        ax[0].axis(ymin=0, ymax=290000 / (1000 ** 2))

        fig.text(0.04, 0.5, 'von Miss stress [MPa]', va='center', rotation='vertical')

        ax[0].set_title('Spanwise von Mises stress')

        ax[1].plot(stress_sim[1, range(4, 28 * len(internal_forces.x_slice), 28)],
                   stress_sim[-1, range(7, 28 * len(internal_forces.x_slice), 28)] / (1000 ** 2) / s,
                   label='Model LE', marker='+', c='darkred', markevery=20,
                   linestyle='--', linewidth=1)

        ax[2].plot(stress_sim[1, range(9, 28 * len(internal_forces.x_slice), 28)],
                   stress_sim[-1, range(9, 28 * len(internal_forces.x_slice), 28)] / (1000 ** 2) / s,
                   label='Model LE', marker='+', c='darkred', markevery=20,
                   linestyle='--', linewidth=1)

        ax[3].plot(stress_sim[1, range(15, 28 * len(internal_forces.x_slice), 28)[0:-1]],
                   stress_sim[-1, range(15, 28 * len(internal_forces.x_slice), 28)[0:-1]] / (1000 ** 2) / s,
                   label='Model LE', marker='+', c='darkred', markevery=20,
                   linestyle='--', linewidth=1)
        ax[3].plot(local_stress_te[1] / 1000, local_stress_te[-1] / (1000 ** 2), label='Validation LE', marker='^',
                   c='darkblue',
                   markevery=20, linewidth=1, )
        ax[3].legend(loc="upper right")
        ax[3].set_xlabel('X - span [m]')

        fig.text(1. - 0.08, 0.20, 'TE', va='center', rotation='vertical')
        fig.text(1. - 0.08, 0.4, '3/4 cord', va='center', rotation='vertical')
        fig.text(1. - 0.08, 0.6, '1/4 cord', va='center', rotation='vertical')
        fig.text(1. - 0.08, 0.8, 'LE', va='center', rotation='vertical')

        plt.show()

    return


def results_deformation():
    nodes = get_nodes()
    def_le_sim, def_te_sim = get_deformation_sim()

    for case in ["ULC1"]: # "UR1", "UR2", "ULC1", "ULC2"
        deformation, deformation_data = get_deformation_val(case)
        def_mag_val = deformation_data[1]
        deformation_data = np.delete(deformation_data, 1, axis=0)  # delete [magnitude def] row
        plot_deformation(deformation,deformation_data,"Numerical model deformation",[True,0])

        # input: returns an array with deformation_data[node label][x def][y def][z def] and
        # input: deformation[node label][x total][y total][z total]

        delete_wired_beam = np.array([])
        for i in range(len(deformation[0])):
            if -65 < deformation[2, i] < 90 and - 30 < deformation[3, i] < 30:
                delete_wired_beam = np.append(delete_wired_beam, i)

        deformation = np.delete(deformation, delete_wired_beam, axis=1)
        deformation_data = np.delete(deformation_data, delete_wired_beam, axis=1)
        nodes = np.delete(nodes, delete_wired_beam, axis=1)
        def_mag_val = np.delete(def_mag_val, delete_wired_beam, axis=0)

        plot_alt_deformation(nodes, deformation, "FEA warping cord length")  # NOTE: rotated complete data 3D plot

        local_def_data = deformation_data
        local_def = deformation
        index_sorted = np.argsort(nodes[3])

        index_le = np.flip(index_sorted, axis=0)[..., 0:68].astype(int)
        local_def_data_le = local_def_data[..., index_le] / 1000.
        local_def_le = local_def[..., index_le] / 1000.

        index_te = index_sorted[..., 0:68].astype(int)
        local_def_data_te = local_def_data[..., index_te] / 1000.
        local_def_te = local_def[..., index_te] / 1000.

        index_x_sorted_le = np.argsort(local_def_le[1])
        local_def_data_le = local_def_data_le[..., index_x_sorted_le]
        local_def_le = local_def_le[..., index_x_sorted_le]

        index_x_sorted_te = np.argsort(local_def_te[1])
        local_def_data_te = local_def_data_te[..., index_x_sorted_te]
        local_def_te = local_def_te[..., index_x_sorted_te]

        # note: def_sim[xtotal][x def][y def][z def][def mag]
        # note: def_val[xtotal][0][y def][z def][def mag]
        # note: for plot LE, TE can be appended

        for i in range(len(local_def_le[2])):
            # print(local_def_le[1,i])
            if abs(local_def_le[1, i] - 0.554) < 0.001:
                z_offset_hinge_2 = local_def_le[2, i]
        set_model_above_vali = 0.1

        front = internal_forces.h_a / 2
        back = internal_forces.c_a - internal_forces.h_a / 2

        def_mag_sim_le = np.sqrt(def_le_sim[1] ** 2 + (def_le_sim[2]) ** 2 + def_le_sim[3] ** 2)  # -z_offset_hinge_2
        def_mag_sim_te = np.sqrt(def_le_sim[1] ** 2 + (def_le_sim[2]) ** 2 + def_le_sim[3] ** 2)

        def_mag_val_le = np.sqrt(local_def_data_le[1] ** 2 + (local_def_data_le[2] + 0.43837 * front) ** 2 + (
            local_def_data_le[3]) ** 2)
        def_mag_val_te = np.sqrt(local_def_data_te[1] ** 2 + (local_def_data_te[2] - 0.43837 * back - 0.01) ** 2 + (
            local_def_data_te[3]) ** 2)

        def_sim_le = np.array([def_le_sim[0], def_le_sim[1], def_le_sim[2], def_le_sim[3], def_mag_sim_le])
        def_sim_te = np.array([def_te_sim[0], def_te_sim[1], def_te_sim[2], def_te_sim[3], def_mag_sim_te])
        def_sim = np.array(
            [np.append(def_le_sim[0], def_te_sim[0]), np.append(def_le_sim[1], def_te_sim[1]),
             np.append(def_le_sim[2], def_te_sim[2] + 0.484),
             np.append(def_le_sim[3] + set_model_above_vali, def_te_sim[3] + set_model_above_vali),
             np.append(def_mag_sim_le, def_mag_sim_te)])

        def_val_le = np.array(
            [local_def_le[1], np.zeros((1, len(local_def_le[1]))), local_def_data_le[2]+ 0.43837 * front, #
             local_def_data_le[3], def_mag_val_le])
        def_val_te = np.array(
            [local_def_te[1], np.zeros((1, len(local_def_te[1]))), local_def_data_te[2]- 0.43837 * back - 0.01 , #
             local_def_data_te[3], def_mag_val_te])

        def_val = np.array([np.append(local_def_le[1], local_def_te[1]),
                            np.append(np.zeros((1, len(local_def_le[1]))), np.zeros((1, len(local_def_te[1])))),
                            np.append(local_def_data_le[2] + 0.43837 * front,
                                      local_def_data_te[2] + 0.484 - 0.43837 * back - 0.01),
                            np.append(local_def_data_le[3] + 0.102 * front, local_def_data_te[3] - 0.102 * front),
                            np.append(def_mag_val_le, def_mag_val_te)])

        # NOTE: MPE ========= MPE ========= MPE ========= MPE ========= MPE ========= MPE
        # NOTE: LE MPE
        sim_points = np.reshape(def_sim_le[0], (len(def_sim_le[0]), 1))
        vali_points = np.reshape(def_val_le[0], (len(def_val_le[0]), 1))
        def_sim_le_MPE = alt_reshape_data_grid(def_sim_le[0:5])
        def_val_le_MPE = alt_reshape_data_grid(def_val_le[0:5])

        MPE_le = np.array([])
        for i in [0, 1, 2]:
            def_sim_at_vali = griddata(sim_points, def_sim_le_MPE[..., i], vali_points, method="cubic", fill_value=0.)
            if np.sum(def_val_le_MPE[..., i]) == 0 and np.sum(def_sim_at_vali) == 0:
                MPE_le = np.append(MPE_le, 0.)
                continue
            elif np.sum(def_val_le_MPE[..., i]) == 0 or np.sum(def_sim_at_vali) == 0:
                print('Upps, exclude if this coms up')
                continue
            MPE_le = np.append(MPE_le, 100. / len(def_sim_at_vali) * (
                    np.sum(def_val_le_MPE[..., i]) - np.sum(def_sim_at_vali)) / np.sum(def_val_le_MPE[..., i]))

        # NOTE: TE MPE

        sim_points = np.reshape(def_sim_te[0], (len(def_sim_te[0]), 1))
        vali_points = np.reshape(def_val_te[0], (len(def_val_te[0]), 1))

        def_sim_te_MPE = alt_reshape_data_grid(def_sim_te[0:5])
        def_val_te_MPE = alt_reshape_data_grid(def_val_te[0:5])

        MPE_te = np.array([])
        for i in [0, 1, 2]:
            def_sim_at_vali = griddata(sim_points, def_sim_te_MPE[..., i], vali_points, method="cubic", fill_value=0.)

            if np.sum(def_val_te_MPE[..., i][0]) == 0. and np.sum(def_sim_at_vali[0]) == 0.:
                MPE_te = np.append(MPE_te, 0.)
                continue
            elif np.sum(def_val_te_MPE[..., i][0]) == 0.:
                print('Upps, exclude if this coms up')
                continue
            MPE_te = np.append(MPE_te, 100. / len(def_sim_at_vali) * (
                    np.sum(def_val_te_MPE[..., i]) - np.sum(def_sim_at_vali)) / np.sum(def_val_te_MPE[..., i]))

        MPE = np.sum(np.abs(MPE_le) + np.abs(MPE_le)) / (2 * 3)
        print("x,y,z averaged deformation MPE at LE & TE: " + str(
            int(10000 * MPE) / 100.) + " %")
        print("x,y,z averaged deformation MPE at LE: " + str(int(1000 * np.sum(np.abs(MPE_le))) / 100.) + " %")
        # note: UNCOMMENT FOR RESULTS

        # NOTE: PLOT ========= PLOT ========= PLOT ========= PLOT ========= PLOT ========= PLOT

        plot_deformation_alternative(def_val, def_sim,
                                     "LE & TE deformation, model data offset in Y for clarity.")  # NOTE: 3D plot here

        fig, ax = plt.subplots(nrows=4, ncols=2, sharex=True, sharey=True)

        ax[0, 0].plot(def_sim_le[0], def_sim_le[4][0] * 1000, label='Model LE', marker='+', c='darkred', markevery=20,
                      linestyle='--', linewidth=1)
        ax[0, 0].plot(def_val_le[0], def_val_le[4] * 1000, label='Validation LE', marker='^', c='darkblue',
                      markevery=20, linewidth=1)

        ax[0, 0].set_ylabel('def. mag. [mm]')
        ax[0, 0].set_title('Deformation LE')

        ax[1, 0].plot(def_sim_le[0], def_sim_le[1][0] * 1000, label='Model LE', marker='+', c='darkred', markevery=20,
                      linestyle='--', linewidth=1)
        ax[1, 0].plot(def_val_le[0], def_val_le[1][0] * 1000, label='Validation LE', marker='^', c='darkblue',
                      markevery=20, linewidth=1, )

        ax[1, 0].set_ylabel('def. x [mm]')
        ax[1, 0].legend(loc="upper right")

        ax[2, 0].plot(def_sim_le[0], def_sim_le[2] * 1000, label='Model LE', marker='+', c='darkred', markevery=20,
                      linestyle='--', linewidth=1)
        ax[2, 0].plot(def_val_le[0], def_val_le[2] * 1000, label='Validation LE', marker='^', c='darkblue',
                      markevery=20, linewidth=1, )

        ax[2, 0].set_ylabel('def. y [mm]')

        ax[3, 0].plot(def_sim_le[0], def_sim_le[3] * 1000, label='Model LE', marker='+', c='darkred', markevery=20,
                      linestyle='--', linewidth=1)
        ax[3, 0].plot(def_val_le[0], def_val_le[3] * 1000, label='Validation LE', marker='^', c='darkblue',
                      markevery=20, linewidth=1, )

        ax[3, 0].set_ylabel('def. z [mm]')
        ax[3, 0].set_xlabel('X - span [m]')
        ax[3, 0].axis(xmin=-0.05 * 1.691, xmax=1.05 * 1.691)

        # NOTE: right side of plot

        ax[0, 1].plot(def_sim_te[0], def_sim_te[4][0] * 1000, label='Model TE', marker='+', c='darkorange',
                      markevery=20, linestyle='--', linewidth=1)
        ax[0, 1].plot(def_val_te[0], def_val_te[4] * 1000, label='Validation TE', marker='^', c='darkcyan',
                      markevery=20, linewidth=1, )

        ax[0, 1].set_title('Deformation TE')

        ax[1, 1].plot(def_sim_te[0], def_sim_te[1][0] * 1000, label='Model TE', marker='+', c='darkorange',
                      markevery=20, linestyle='--', linewidth=1)
        ax[1, 1].plot(def_val_te[0], def_val_te[1][0] * 1000, label='Validation TE', marker='^', c='darkcyan',
                      markevery=20, linewidth=1, )

        ax[1, 1].legend(loc="upper right")

        ax[2, 1].plot(def_sim_te[0], def_sim_te[2] * 1000, label='Model TE', marker='+', c='darkorange', markevery=20,
                      linestyle='--', linewidth=1)
        ax[2, 1].plot(def_val_te[0], def_val_te[2] * 1000, label='Validation TE', marker='^', c='darkcyan',
                      markevery=20, linewidth=1, )

        ax[3, 1].plot(def_sim_te[0], def_sim_te[3] * 1000, label='Model TE', marker='+', c='darkorange', markevery=20,
                      linestyle='--', linewidth=1)
        ax[3, 1].plot(def_val_te[0], def_val_te[3] * 1000, label='Validation TE', marker='^', c='darkcyan',
                      markevery=20, linewidth=1, )

        ax[3, 1].set_xlabel('X - span [m]')
        ax[3, 1].axis(xmin=-0.05 * 1.691, xmax=1.05 * 1.691)
        plt.show()
        return


def validate(): #TODO run this for validation results
    # plot_node_pos()
    # for i in ["UR1", "UR2", "ULC1", "ULC2"]:
    #     plot_deformation(i)
    #
    # for i in ["SR1", "SR2", "SLC1", "SLC2"]:
    #     stress = get(i)
    #     plot_von_mieses_stress(stress)
    # get_elements('elements')
    # __, __ = get_deformation_sim()
    results_stress()
    results_deformation()
    return


# </editor-fold>

# validate()
