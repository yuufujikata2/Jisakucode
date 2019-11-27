import numpy as np

def polargrid(dr):
    ntheta = 100
    nphi = 200
    coordx = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    coordy = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    coordz = np.zeros(nphi * (ntheta + 1) - 2 * (nphi - 1) )
    number = 0
    for phi in range (nphi):
        if phi == 0 :
            for theta in range (ntheta + 1):
                coordx[number] = dr * np.sin(np.pi / ntheta * theta) * np.cos(2 * np.pi / nphi *  phi )
                coordy[number] = dr * np.sin(np.pi / ntheta * theta) * np.sin(2 * np.pi / nphi *  phi )
                coordz[number] = dr * np.cos(np.pi / ntheta * theta)
                number += 1

        else:
            for theta in range (1,ntheta):
                coordx[number] = dr * np.sin(np.pi / ntheta * theta) * np.cos(2 * np.pi / nphi *  phi )
                coordy[number] = dr * np.sin(np.pi / ntheta * theta) * np.sin(2 * np.pi / nphi *  phi )
                coordz[number] = dr * np.cos(np.pi / ntheta * theta)
                number += 1
    coord = np.array([coordx,coordy,coordz])
    coord = coord.T
    return coord,number-1

