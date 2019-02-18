import numpy as np

# important x-locations along the beam
x_begin = 0
x_1 = 0.149
x_2 = 0.554
x_3 = 1.541
x_end = 1.691

# actuator distances
x_a = 0.272
x_a1 = x_2 - x_a/2
x_a2 = x_2 + x_a/2

# other parameters
l_a = 1.691
q = 2710
p = 37900
E = 1

x_coordinates = [x_begin, x_1, x_a1, x_2, x_a2, x_3, x_end]


y_forces = {'y_1': [1, x_1], 'y_2': [1, x_2], 'y_3': [1, x_3]}
z_forces = {'z_1': [1, x_1], 'a_1': [1, x_a1], 'z_2': [1, x_2], 'a_2': [p, x_a2], 'z_3': [1, x_3]}

sum_forces_y = []
sum_moments_x = []

for force in y_forces:
    sum_forces_y.append(y_forces[force][0])
    sum_moments_x.append(y_forces[force] * (y_forces[force][1] - x_begin))

















