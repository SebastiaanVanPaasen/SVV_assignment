# <editor-fold desc="INIT">
import numpy as np
import math
import matplotlib.pyplot as plt
# </editor-fold>


class SimulationData:
    IdealizedStructure_n_discretize = 601  # only ODD integer, 13<n<600 [contour_length/w_stiffner]


class ProblemData:  # Storage for globally accessiblee data (input geometry,..), USE mKgs system!
    t_skin, h_aileron, cord_aileron = 0.0011, 0.173, 0.484
    t_stiffner, h_stiffner, w_stiffner = 0.0012, 0.014, 0.018
    t_spar, h_spar = 0.0025, h_aileron

#maybe create a class with general/utility functions like this one:
def read_table(filename):  # tables must be in form: 1 row header (description), separation by tabulators
    file = open(filename, "r+")
    file = file.readlines()
    n_colums = len(file[0].split('\t'))
    del file[0]
    n_rows = len(file)
    data=np.zeros(shape=(n_colums,n_rows))
    for i in range(n_rows):
        l = file[i].split('\t')
        for j in range(n_colums):
            data[j,i]=l[j]
    return data
