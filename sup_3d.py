from axes import set_axes_equal
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import surf2stl
from normals import normal_vec

# PARAMETERS
###############################################

L = 20
L_inf = 16
r_min = 0.5
r_med = 2
r_max = 2.5
c = 2  # standard desviation
folds_num = 4
thickness = 0.2
dis = 100
name = 'Test_05_05_2022.stl'
surf_tot = 1 # 1 to set x range [-L/2,L/2], th range [0,2pi]
             # 2 to set x range [0,L/2], th range [0,pi]

# RADIUS FUNCTION r(x,θ)
# AND PARTIAL DERIVATES dr/dx, dr/dθ
###############################################

def r(x_, th_) :

    gauss1 = (r_med - r_min)*exp(-(x_**2/2/c**2))
    gauss2 = (r_max - r_med)*exp(-(x_**2/2/c**2))

    r = gauss1 + gauss2*cos(folds_num*th_) + r_min

    return r

# 3D SPACE
###############################################

fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

# SET X AND THETA
###############################################

if surf_tot == 1 :
    x = linspace(-L/2, L/2, dis)
    x_i = linspace(L_inf/2, -L_inf/2, dis)

    th = linspace(0, 2*pi, dis)
    th_i = linspace(0, 2*pi, dis)
else :
    x = linspace(0, L/2, dis)
    x_i = linspace(0, L_inf/2, dis)

    th = linspace(0, pi, dis)
    th_i = linspace(0, pi, dis)

# CALCULATE NORMAL VECTORS
###############################################

n_x, n_y, n_z, n_x_inf, n_y_inf, n_z_inf = normal_vec(surf_tot,L,L_inf,r_min,r_med,r_max,c,folds_num,thickness,dis)

# SET X, Y AND Z
###############################################

th, x = meshgrid(th, x)
th_i, x_i = meshgrid(th_i, x_i)

rad = r(x,th)
rad_i = r(x_i,th_i)

# MEDIUM SURFACE

x_med = x
y_med = rad * cos(th)
z_med = rad * sin(th)

# INFERIOR SURFACE

x_inf = x_i - (thickness/2 * array(n_x_inf))
y_inf = rad_i * cos(th) - (thickness/2 * array(n_y_inf))
z_inf = rad_i * sin(th) - (thickness/2 * array(n_z_inf))

# SUPERIOR SURFACE

x_sup = x + (thickness/2 * array(n_x))
y_sup = y_med + (thickness/2 * array(n_y))
z_sup = z_med + (thickness/2 * array(n_z))

# CLOSE SURFACE
###############################################

if surf_tot == 1 :
    u = linspace(0,  2*pi, 100)
else :
    u = linspace(0,  pi, 100)

R = linspace(0, r_min+thickness/2, 100)
h = L/2
x_c = h * outer(ones(size(u)), ones(size(u)))
y_c = outer(R, cos(u))
z_c = outer(R, sin(u) )

R_i = linspace(0, r_min-thickness/2, 100)
h_i = L_inf/2
x_c_i = h_i * outer(ones(size(u)), ones(size(u)))
y_c_i = outer(R_i, cos(u))
z_c_i = outer(R_i, sin(u) )

# WRITE .STL FILES
###############################################

#surf2stl.write('internal.stl', x_inf, y_inf, z_inf)
#surf2stl.write('external1.stl', x_sup, y_sup, z_sup)

surf2stl.multiwrite(name,[x_inf,x_sup,x_c,-x_c,x_c_i,-x_c_i],[y_inf,y_sup,y_c,y_c,y_c_i,y_c_i],[z_inf,z_sup,z_c,z_c,z_c_i,z_c_i])

# 3D PLOT
###############################################
ax.plot_surface(x_inf,y_inf,z_inf)
#ax.plot_surface(x_med,y_med,z_med)
ax.plot_surface(x_sup,y_sup,z_sup)

if surf_tot == 1 :
    ax.plot_surface(x_c,y_c,z_c)
    ax.plot_surface(-x_c,y_c,z_c)
    ax.plot_surface(x_c_i,y_c_i,z_c_i)
    ax.plot_surface(-x_c_i,y_c_i,z_c_i)

else :
    ax.plot_surface(x_c,y_c,z_c)
    ax.plot_surface(x_c_i,y_c_i,z_c_i)

set_axes_equal(ax)
plt.show()
