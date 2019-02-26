
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# <editor-fold desc="SCRAP">
# def clist(c):
#     lst = []
#     for i in range(100):
#         lst.append(((c[0] - i * 2.) / 255, (c[1] + i * 2.) / 255., (c[2] - i * 2.) / 255))
#     return lst
# </editor-fold>


# <editor-fold desc="GET METHODS">
def read_table_vali(filename,
                    replace_comma_with_dot):  # tables must be in form: 1 row header (description), separation by tabulators
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


def get_deformation(case):
    # output: returns an array with deformation_data[node label][magnitude def][x def][y def][z def] and deformation[node label][x total][y total][z total]
    filename = str(case) + ".txt"
    nodes = read_table_vali("inp.txt", False)
    deformation_data = read_table_vali(filename, True)
    deformation_data = np.delete(deformation_data, [1, 2, 3, 4], 0)
    print(deformation_data)

    x = nodes[1] + deformation_data[2]
    y = nodes[2] + deformation_data[3]
    z = nodes[3] + deformation_data[4]

    deformation = np.array([nodes[0], x, y, z])  # add empty subarray at [0] for convienience (plotting)
    return deformation, deformation_data


def get_von_mieses_stress(
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


# </editor-fold>


# <editor-fold desc="PLOT METHODS">
def plot_node_pos():
    nodes = read_table_vali("inp.txt", False)
    color = []
    for i in nodes[0]:
        color.append((0., i / max(nodes[0]), 0.))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(nodes[1], nodes[2], nodes[3], zdir="z", s=2, c=color)

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


def plot_deformation(case):  # use case variable like "UR1", only will work with U...
    deformation, deformation_data = get_deformation(case)

    color = []
    for i in deformation_data[1]:
        color.append((0., i / max(deformation_data[1]), 0.))

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    ax.scatter(deformation[1], deformation[2], deformation[3], zdir="z", s=2, c=color)

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
    plt.show()
    return


def plot_von_mieses_stress(case):
    stress = get_von_mieses_stress(case)

    color = []
    for i in stress[4]:
        color.append((0., i / max(stress[4]), 0.))

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
    plt.show()
    return
# </editor-fold>


# <editor-fold desc="VALIDATION METHODS">
#TODO: {switch coordinate frames}, interpolate, compute new value at comparison points, calculate MPE

def coordiant_transform(x,y,z):
    x = x + 0
    y = y + 0
    z = z + 0
    return x,y,z # return with switched coordinates as needed

def validate():
    # plot_node_pos()
    # for i in ["UR1", "UR2", "ULC1", "ULC2", ]:
    #     plot_deformation(i)
    for i in ["SR1", "SR2", "SLC1", "SLC2", ]:
        plot_von_mieses_stress(i)
    # get_elements('elements')

# </editor-fold>




validate()
