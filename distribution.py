import numpy as np


def torque_shear(torque, area):
    shear_flow = torque / (2 * area)

    return shear_flow


