import numpy as np

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
E = 1


class Force:

    def __init__(self, magnitude, direction, position):
        self.magnitude = magnitude
        self.direction = direction
        self.position = position
#        self.x = magnitude*direction[0]
#        self.y = magnitude*direction[1]
#        self.z = magnitude*direction[2]

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
#system = []
sum_forces_y = []
sum_forces_z = []
sum_moments = []

for force in forces:
    sum_forces_y.append(force.determine_force('y'))
    sum_forces_z.append(force.determine_force('z'))
    sum_moments.append(force.determine_moment([0, 0, 0]))

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

#system.append(sum_forces_y)
#system.append(sum_forces_z)
#for j in range(3):
#    for i in range(len(sum_moments)):
#        system.append(sum_moments[i][j])
#system.append(ce_eq_y)
#system.append(ce_eq_z)

#solving reaction forces
a = np.zeros((7,7))
unk = np.zeros(7)
b = np.zeros(7)

for i in range(7):
    a[0][i] = sum_forces_y[i]
    a[1][i] = sum_forces_z[i]
    a[2][i] = sum_moments[i][0]
    a[3][i] = sum_moments[i][1]
    a[4][i] = sum_moments[i][2]
    a[5][i] = ce_eq_y[i]
    a[6][i] = ce_eq_z[i]
    
b[0] = sum_forces_y[-1] + sum_forces_y[-2]
b[1] = sum_forces_z[-1] + sum_forces_z[-2]
b[2] = sum_moments[0][-1] + sum_moments[0][-2]
b[3] = sum_moments[1][-1] + sum_moments[1][-2]
b[4] = sum_moments[2][-1] + sum_moments[2][-2]
b[5] = ce_eq_y[-1] + ce_eq_y[-2]
b[6] = ce_eq_z[-1] + ce_eq_z[-2]

unk = np.linalg.solve(a,b)
