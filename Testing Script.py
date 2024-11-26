import numpy as np

def panel_normal(tilt):     # it calculates the normal vector
    u = np.array([0,0,1])
    v = np.array([-np.cos(tilt),np.sin(tilt),0])
    return np.cross(v,u)


def test_fun(tilt):
    return np.array([np.sin(tilt), np.cos(tilt), 0])


print(panel_normal(3))
print(test_fun(3))
