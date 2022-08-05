from numpy import *

def normal(x_, th_,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis) :

    r = (r_med - r_min)*exp(-(x_**2/2/c**2)) + (r_max - r_med)*exp(-(x_**2/2/c**2))*cos(folds_num*th_) + r_min
    dr_dx = -x_*(r_max - r_med)*exp(-x_**2/(2*c**2))*cos(folds_num*th_)/c**2 - x_*(r_med - r_min)*exp(-x_**2/(2*c**2))/c**2
    dr_dth = -folds_num*(r_max - r_med)*exp(-x_**2/(2*c**2))*sin(folds_num*th_)

    X = [x_, r*cos(th_), r*sin(th_)]
    dX_dx = [1, dr_dx*cos(th_), dr_dx*sin(th_)]
    dX_dth = [0, dr_dth*cos(th_) - r*sin(th_), dr_dth*sin(th_) + r*cos(th_)]

    #n = cross(dX_dx,dX_dth)
    n = cross(dX_dth, dX_dx)
    n_un = n / linalg.norm(n)

    return n_un

def normal_vec(surf_tot, L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis) :
    if surf_tot == 1 :
        x = linspace(-L/2, L/2, dis)
        x_i = linspace(-L_inf/2, L_inf/2, dis)

        th = linspace(0, 2*pi, dis)
        th_i = linspace(0, 2*pi, dis)
    else :
        x = linspace(0, L/2, dis)
        x_i = linspace(0, L_inf/2, dis)

        th = linspace(0, pi, dis)
        th_i = linspace(0, pi, dis)

    n_x = []
    n_y = []
    n_z = []

    for i in x :
        n_x_i = []
        n_y_i = []
        n_z_i = []
        for j in th :
            n_x_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[0])
            n_y_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[1])
            n_z_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[2])
        n_x.append(n_x_i)
        n_y.append(n_y_i)
        n_z.append(n_z_i)

    # NORMAL VECTORS FOR INF. SURFACE

    n_x_inf = []
    n_y_inf = []
    n_z_inf = []

    for i in x_i :
        n_x_inf_i = []
        n_y_inf_i = []
        n_z_inf_i = []
        for j in th_i :
            n_x_inf_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[0])
            n_y_inf_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[1])
            n_z_inf_i.append(normal(i,j,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)[2])
        n_x_inf.append(n_x_inf_i)
        n_y_inf.append(n_y_inf_i)
        n_z_inf.append(n_z_inf_i)

    return n_x, n_y, n_z, n_x_inf, n_y_inf, n_z_inf
