import numpy as np

def make_environment():

    #default value

    #potential
    pot_type = "cubic"
    pot_bottom = -1
    pot_show_f = False
    pot_region = np.array((2/np.sqrt(3),2/np.sqrt(3),2/np.sqrt(3)))
    r2_ratio = 0.5
    bound_rad = pot_region * r2_ratio
    radius = np.sqrt(pot_region[0] **2 + pot_region[1] **2 + pot_region[2] **2 )
    region = (radius,radius,radius)
    gridpx,gridpy,gridpz = 200,200,200

    #logmesh
    nr = 201
    a = 0.01
    b = radius / (np.e**(a * ( nr - 1)) - 1)
    rofi = np.array([b * (np.e**(a * i) - 1) for i in range(nr)])

    #surface integral
    si_method = "lebedev_py"
    radial_pot_show_f = False
    new_radial_pot_show_f = False

    #basis
    node_open = 1
    node_close = 2
    LMAX = 4

    #read input
    try:
        fr = open("input",mode="r")
        lines = fr.readlines()
        for i in range(len(lines)):
            line = lines[i].split()
            if line[0] in "#":
                continue
            if line[0] == "pot_type":
                pot_type = line[1]
            if line[0] == "pot_region":
                pot_region = np.array((float(line[1]),float(line[2]),float(line[3])))
            if line[0] == "r2_ratio":
                r2_rario = float(line[1])
            if line[0] == "pot_mesh_point":
                gridpx,gridpy,gridpz = int(line[1]),int(line[2]),int(line[3])
            if line[0] == "pot_show_f":
                if line[1] == "True":
                    pot_show_f = True
            if line[0] == "logmesh_point":
                nr = float(line[1])
            if line[0] == "logmesh_a":
                a = float(line[1])
            if line[0] == "si_method":
                si_method = line[1]
            if line[0] == "radial_pot_show_f":
                if line[1] == "True":
                    radial_pot_show_f = True
            if line[0] == "new_radial_pot_show_f":
                if line[1] == "True":
                    new_radial_pot_show_f = True
            if line[0] == "node_open":
                node_open = int(line[1])
            if line[0] == "node_close":
                node_close = int(line[1])
            if line[0] == "LMAX":
                LMAX = int(line[1]) + 1
            if line[0] == "a":
                a = float(line[1])
    
    except FileNotFoundError:
        print("Waring! There is no input file. All values are default.")

    #make environment
    bound_rad = pot_region * r2_ratio
    x,y,z = grid(gridpx,gridpy,gridpz,region)
    xx,yy,zz = np.meshgrid(x,y,z)
    radius = np.sqrt(pot_region[0] **2 + pot_region[1] **2 + pot_region[2] **2 )
    region = (radius,radius,radius)
    b = radius / (np.e**(a * ( nr - 1)) - 1)
    rofi = np.array([b * (np.e**(a * i) - 1) for i in range(nr)])



    return pot_region,bound_rad,radius,region,nr,gridpx,gridpy,gridpz,x,y,z,xx,yy,zz,a,b,rofi,pot_type,pot_bottom,pot_show_f,si_method,radial_pot_show_f,new_radial_pot_show_f,node_open,node_close,LMAX


def grid (nx,ny,nz,region):
    x = np.linspace(-region[0],region[0],nx)
    y = np.linspace(-region[1],region[1],ny)
    z = np.linspace(-region[2],region[2],nz)
    return x,y,z


