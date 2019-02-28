import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import bottleneck as bn


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
    # output: returns an array with deformation_data[node label][magnitude def][x def][y def][z def] and
    # output: deformation[node label][x total][y total][z total]
    filename = str(case) + ".txt"
    nodes = read_table_vali("inp.txt", False)
    deformation_data = read_table_vali(filename, True)
    deformation_data = np.delete(deformation_data, [1, 2, 3, 4], 0)

    x = nodes[1] + deformation_data[2]
    y = nodes[2] + deformation_data[3]
    z = nodes[3] + deformation_data[4]

    deformation = np.array([nodes[0], x, y, z])  # add node[0] subarray at [0] for convienience (plotting)
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


def get_nodes():
    nodes = read_table_vali("inp.txt", False)
    return nodes


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
    ax.scatter(LE[1], LE[2], LE[3], zdir="z", s=10, c="red")
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


def plot_deformation(deformation, deformation_data, title, c):  # use case variable like "UR1", only will work with U...
    if c[0]:
        colooor = []
        for i in deformation_data[1]:
            colooor.append((1.-i / max(deformation_data[1]), 0., 0.))
    else:
        choice = ["red", "blue", "green", "cyan", "orange"]
        colooor = choice[c[1]]

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    colooor=np.ravel(colooor).reshape((len(deformation[0]),3))
    ax.scatter(deformation[1], deformation[2], deformation[3], zdir="z", s=10, c=colooor)

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


def plot_von_mieses_stress(stress):
    print(stress)
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
def reshape_data_grid(data):  # must be 4D array, discards 1st subarray
    flat_data = np.array([])
    for i in range(len(data[0])):
        flat_data = np.append(flat_data, np.array([data[1, i], data[2, i], data[3, i]]))
    data = np.reshape(flat_data, (len(data[0]), 3), order="A")
    return data


# TODO: {switch coordinate frames}, interpolate, compute new value at comparison points, calculate MPE

def coordiant_transform(x, y, z):
    x = x + 0
    y = y + 0
    z = z + 0
    return x, y, z  # return with switched coordinates as needed


def results_deformation():
    nodes = get_nodes()
    for case in ["UR1", "UR2", "ULC1", "ULC2", ]:
        deformation, deformation_data = get_deformation(case)
        deformation_data = np.delete(deformation_data, 1, axis=0)
        combined_edges_def, combined_edges_def_data = np.array([]), np.array([])
        for i in [0, 1, 2]:  # note: for x,y,z axis view
            i += 1
            local_def = deformation
            local_def_data = deformation_data
            index_sorted = np.argsort(nodes[i])

            index_LE = np.flip(index_sorted, axis=0)[..., 0:68].astype(int)
            local_def_LE = local_def[..., index_LE]
            local_def_data_LE = local_def_data[..., index_LE]

            index_TE = index_sorted[..., 0:68].astype(int)
            local_def_TE = local_def[..., index_TE]
            local_def_data_TE = local_def_data[..., index_TE]

            combined_edges_def = np.append(combined_edges_def, local_def_LE)
            combined_edges_def_data = np.append(combined_edges_def_data, local_def_data_LE)
            combined_edges_def = np.append(combined_edges_def, local_def_TE)
            combined_edges_def_data = np.append(combined_edges_def_data, local_def_data_TE)

        print(deformation)
        print(deformation_data)
        plot_deformation(deformation, deformation_data, case + ": Full 3D", [True, 0])  # complete data 3D plot
        plot_deformation(combined_edges_def, combined_edges_def_data, case + ": LE & TE 3D plot",
                         [True, 0])  # LE,TE 3D plot

        f, ax = plt.subplots(2, 2)

        ax[0, 0].plot(local_def[0], local_def[1])
        ax[0, 0].set_title('Axis [0,0]')

        ax[0, 1].scatter(local_def[0], local_def[2])
        ax[0, 1].set_title('Axis [0,1]')

        ax[1, 0].plot(local_def[1], local_def[2])
        ax[1, 0].set_title('Axis [1,0]')

        # ax[1, 1].plot_deformation(deformation, deformation_data, case + ": Full 3D", [True, 0])  # complete data 3D plot
        # ax[1, 1].set_title('Axis [1,1]')

        plt.show()

    # for case in ["UR1", "UR2", "ULC1", "ULC2", ]:
    #     for j in [0, 1, 2]:
    #         for i in range(2, 4):  # note: prints deformation LE, TE viewed from y axis or z axis depending on j
    #             nodes = get_nodes()
    #             deformation, deformation_data = get_deformation(case)
    #             deformation_data = np.delete(deformation_data, 1, axis=1)
    #             index_sorted = np.argsort(nodes[i])
    #             # print(index_sorted)
    #             # print(nodes[0,index_sorted])
    #             index_LE = np.flip(index_sorted, axis=0)[..., 0:68].astype(int)
    #             index_TE = index_sorted[..., 0:68].astype(int)
    #             #
    #             # print(index_LE)
    #             # print(index_TE)
    #             #
    #             # print(deformation)
    #             deformation = deformation[..., index_LE]
    #             deformation_data = deformation_data[..., index_LE]
    #
    #             x_def = deformation[1]
    #             y_def = deformation[2]
    #             z_def = deformation[3]
    #
    #             plot_deformation(deformation, deformation_data)
    #
    #             if j == 1:
    #                 plt.scatter(x_def, y_def, label="x_def, y_def")
    #             elif j == 2:
    #                 plt.scatter(x_def, z_def, label="x_def, z_def")
    #             else:
    #                 plt.scatter(y_def, z_def, label="y_def, z_def")
    #         # plot_deformation(deformation, deformation_data)
    #         # plt.show()
    # plt.show()

    # <editor-fold desc="SCRAP">
    # LE = np.array([np.flip(nodes[0][index_sorted],axis=0),np.flip(nodes[1][index_sorted],axis=0),np.flip(nodes[2][index_sorted],axis=0),np.flip(nodes[3][index_sorted],axis=0)])
    # TE = np.array([nodes[0][index_sorted], nodes[1][index_sorted],nodes[2][index_sorted],nodes[3][index_sorted]])
    # # print(LE)
    # # print(TE)
    #
    #
    # LE = LE[...,0:68]
    # TE = TE[..., 0:68]
    # print(LE)
    # print(TE)
    # #print(nodes_sorted)
    # plot_node_pos(nodes,LE,TE)

    # nodes=reshape_data_grid(get_nodes())
    # deformation, deformation_data = get_deformation("UR1")
    # deformation = np.delete(deformation_data,1,axis=1)
    #
    # deformation = reshape_data_grid(deformation)

    # test_points= np.random.rand(len(get_nodes()[0]),3)
    # test_points= test_points + nodes
    #
    # nodes = np.delete(nodes,-1,axis=0)
    # test_points = np.delete(test_points,-1,axis=0)
    # # print(len(nodes))
    # # print(len(deformation))
    # # print(len(test_points))
    # def_interpoled=griddata(nodes,deformation,test_points)
    # print(def_interpoled)
    # </editor-fold>
    return


def validate():
    # plot_node_pos()
    # for i in ["UR1", "UR2", "ULC1", "ULC2"]:
    #     plot_deformation(i)
    #
    # for i in ["SR1", "SR2", "SLC1", "SLC2"]:
    #     stress = get_von_mieses_stress(i)
    #     plot_von_mieses_stress(stress)
    # get_elements('elements')

    results_deformation()

    return


# </editor-fold>


validate()