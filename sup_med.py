from axes import set_axes_equal
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# PARAMETERS
###############################################

L = 16
r_min = 0.25
r_med = 2.0
r_max = 2.67
c = 1  # standard desviation
folds_num = 4
thickness = 0.2
dis = 200
file_name = 'prueba.npz'

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

x = linspace(-L/2, L/2, dis)
th = linspace(0, 2*pi, dis)

# SET X, Y AND Z
###############################################

th, x = meshgrid(th, x)

rad = r(x,th)

# MEDIUM SURFACE

x_med = x
y_med = rad * cos(th)
z_med = rad * sin(th)

# 3D PLOT
###############################################

ax.plot_surface(x_med,y_med,z_med)

set_axes_equal(ax)
plt.show()

savez(file_name,x_med=x_med,y_med=y_med,z_med=z_med,L=L,r_min=r_min,r_med=r_med,r_max=r_max,c=c,folds_num=folds_num,thickness=thickness,dis=dis)
