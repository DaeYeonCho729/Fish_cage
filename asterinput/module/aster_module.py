import numpy as np

# Three function used by code_aster
def get_position_aster(table_aster):
    content = table_aster.EXTR_TABLE()
    original_x = content.values()['COOR_X']
    original_y = content.values()['COOR_Y']
    original_z = content.values()['COOR_Z']
    delta_x = content.values()['DX']
    delta_y = content.values()['DY'] 
    delta_z = content.values()['DZ']
    position = np.array([original_x, original_y, original_z]) + np.array([delta_x, delta_y, delta_z])
    return np.transpose(position)

def get_velocity_aster(table_aster):  # to get the velocity
    content = table_aster.EXTR_TABLE()
    velocity_x = content.values()['DX']
    velocity_y = content.values()['DY']
    velocity_z = content.values()['DZ']
    velocity = np.array([velocity_x, velocity_y, velocity_z])
    return np.transpose(velocity)

def get_displace_vector(table_aster):
    content = table_aster.EXTR_TABLE()
    delta_x = content.values()['DX']
    delta_y = content.values()['DY']
    delta_z = content.values()['DZ']
    return np.mean(np.transpose(np.array([delta_x, delta_y, delta_z])),axis=0)

def get_N_aster(table_aster):
    content = table_aster.EXTR_TABLE()
    # N = content.values()['N']
    N_x = content.values()['DX']
    N_y = content.values()['DY']
    N_z = content.values()['DZ']
    N = np.array([N_x, N_y, N_z])
    return np.transpose(N)