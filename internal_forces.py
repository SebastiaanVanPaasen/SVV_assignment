import numpy as np
from SVV_assignment.SVV_assignment.geometric_properties import *

# aileron parameters
l_a = 1.691
h_a = 0.173

# force locations along the beam
x_begin = 0
x_1 = 0.149
x_2 = 0.554
x_3 = 1.541
x_end = 1.691

# actuator distances
x_a = 0.272
x_a1 = x_2 - x_a/2
y_a1 = np.cos(np.radians(h_a/2 - np.tan(np.radians(28))*h_a/2))
z_a1 = np.cos(np.radians(28))*h_a/2 + np.sin(np.radians(28))*h_a/2
x_a2 = x_2 + x_a/2

# load parameters
q = 2710
p = 37900
E = 73.1*10**9
G = 28*10**9
I_zz, I_yy = main()
delta_1 = 0.00681


class Force:

    def __init__(self, magnitude, direction, position):
        self.magnitude = magnitude
        self.direction = direction
        self.position = position

    def determine_force(self, direction):
        determinator = 0
        if direction == 'y':
            determinator = 1
        elif direction == 'z':
            determinator = 2

        if self.direction[determinator] != 0:
            return self.magnitude*self.direction[determinator]
        else:
            return 0

    def determine_moment(self, position):
        distance = self.position - position
        moments = np.zeros(3)

        moments[0] = self.magnitude*-1*(self.direction[1]*distance[2] + self.direction[2]*distance[1])
        moments[1] = self.magnitude*-1*(self.direction[0]*distance[2] + self.direction[2]*distance[0])
        moments[2] = self.magnitude*-1*(self.direction[0]*distance[1] + self.direction[1]*distance[0])

        return moments
    
    #rotation specifically in yz axis
    def rotate(self, angle):
        self.direction[1] = self.direction[1]*np.sin(angle)
        self.direction[2] = self.direction[2]*np.cos(angle)
        
        
y_1 = Force(1, np.array([0, 1, 0]), np.array([x_1, 0, 0]))
y_2 = Force(1, np.array([0, 1, 0]), np.array([x_2, 0, 0]))
y_3 = Force(1, np.array([0, 1, 0]), np.array([x_3, 0, 0]))
z_1 = Force(1, np.array([0, 0, 1]), np.array([x_1, 0, 0]))
z_2 = Force(1, np.array([0, 0, 1]), np.array([x_2, 0, 0]))
z_3 = Force(1, np.array([0, 0, 1]), np.array([x_3, 0, 0]))
a_1 = Force(1, np.array([0, 0, 1]), np.array([x_a1, y_a1, z_a1]))
q_ = Force(q*l_a, np.array([0, -1, 0]), np.array([l_a/2, 0, 0]))
p_ = Force(p, np.array([0, 0, -1]), np.array([x_a2, y_a1, z_a1]))

forces = [y_1, y_2, y_3, z_1, z_2, z_3, a_1, q_, p_]

sum_forces_y = []
sum_forces_z = []
#sum_moments = []
sum_moments_x = []
sum_moments_y = []
sum_moments_z = []

for force in forces:
    sum_forces_y.append(force.determine_force('y'))
    sum_forces_z.append(force.determine_force('z'))
    sum_moments_x.append(force.determine_moment([0, 0, 0])[0])
    sum_moments_y.append(force.determine_moment([0, 0, 0])[1])
    sum_moments_z.append(force.determine_moment([0, 0, 0])[2])
    
#print(sum_forces_y)
#print(sum_forces_z)
#print(sum_moments)

ce_eq_z = np.array([(l_a - x_1)**3/3, 
                    (l_a - x_2)**3/3 - (x_1+x_2)*(l_a - x_2)**2/2 + x_1*x_2*(l_a - x_2), 
                    (l_a - x_3)**3/3 - (x_1),
                    0.,0.,0.,0.,0.,1.])

ce_eq_y = np.array([0.,0.,0.,(l_a - x_1)**3/3,
                    (l_a - x_2)**3/3 - (x_1+x_2)*(l_a - x_2)**2/2 + x_1*x_2*(l_a - x_2),
                    (l_a - x_3)**3/3 - (x_1+x_3)*(l_a - x_3)**2/2 + x_1*x_3*(l_a - x_3),
                    (l_a - x_a1)**3/3 - (x_1+x_a1)*(l_a - x_a1)**2/2 + x_1*x_a1*(l_a - x_a1),
                    (l_a - x_a2)**3/3 - (x_1+x_a2)*(l_a - x_a2)**2/2 + x_1*x_a2*(l_a - x_a2),0.])

system = [sum_forces_y, sum_forces_z, sum_moments_x, sum_moments_y, sum_moments_z, ce_eq_y, ce_eq_z]

#solving reaction forces
sys_mat = np.zeros((7,7))
unk = np.zeros(7)
sys_vec = np.zeros(7)
for i in range(len(system)):
    for j in range(len(system[0])-2):
        sys_mat[i][j] = system[i][j]
    sys_vec[i] = system[i][-1] + system[i][-2]
    
unk = np.linalg.solve(sys_mat,sys_vec)
for i in range(len(forces)-2):
    forces[i].magnitude *= unk[i]
