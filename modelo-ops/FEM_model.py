from axes import set_axes_equal
from numpy import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import openseespy.opensees as ops
from calcular_volumen import calcular_volumen

#import time
#inicio = time.time()

# Código a medir


# PARAMETERS
###############################################

# SHAPE PARAMETERS

sacale_factor = 1

L = 16 * sacale_factor
r_min = 0.5 * sacale_factor
r_med = 2.0 * sacale_factor
r_max = 2.67 * sacale_factor

c = 1 * sacale_factor # standard desviation
folds_num = 3
thickness = 0.2 * sacale_factor
dis = 125

# MATERIAL PARAMETERS

E = 1980   #MPa  (N/mm)
nu = 0.15
rho = 0.

# CHECK r_min, r_med, r_max
###############################################

if r_min > r_med :
    r_min = r_med

if r_med > r_max :
    r_med = r_max

# 3D SURFACE
###############################################

# RADIUS FUNCTION r(x,θ)

def r(x_, th_) :

    if (r_med - r_min) > (r_max - r_med) :

        gauss1 = (r_med - r_min)*exp(-(x_**2/2/c**2))
        gauss2 = (r_max - r_med)*exp(-(x_**2/2/c**2))

    else :

        gauss1 = (r_max - r_med)*exp(-(x_**2/2/c**2))
        gauss2 = (r_med - r_min)*exp(-(x_**2/2/c**2))

    r = gauss1 + gauss2*cos(folds_num*th_) + r_min

    return r

# 3D SPACE

fig = plt.figure()
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

# SET X AND THETA

x = linspace(-L/2, L/2, dis)
th = linspace(0, 2*pi, dis)

# SET X, Y AND Z

th, x = meshgrid(th, x)

rad = r(x,th)

# MEDIUM SURFACE

x_med = x
y_med = rad * cos(th)
z_med = rad * sin(th)

# 3D PLOT

#ax.plot_surface(x_med,y_med,z_med)

#set_axes_equal(ax)
#plt.show()

# FEM MODEL
###############################################

#from numpy import load, zeros, cross

# DATA

x_med = x_med
y_med = y_med
z_med = z_med

cm_to_mm = 10.
thickness = thickness * cm_to_mm

Ni = x_med.shape[0]
Nj = x_med.shape[1]

# INITIALIZE 3D -> 6 GDL PER NODE

ops.model('basic', '-ndm', 3, '-ndf', 6)

# CREATE NYLON

material_tag = 1
ops.nDMaterial('ElasticIsotropic', material_tag, E, nu, rho)

# CREATE SECTION FOR ELEMENTS (THICKNESSES)

section_tag = 1
#ops.section('PlateFiber', section_tag, material_tag, thickness)
ops.section('ElasticMembranePlateSection', section_tag, E, nu, thickness, rho)

node_numbers = zeros((Ni, Nj), dtype=int)

#print(f"Ni = {Ni} Nj = {Nj}")

# ASSIGN NODE NUMBER

node_number = 1
for i in range(Ni):
    for j in range(Nj):
        x = x_med[i,j] * cm_to_mm
        y = y_med[i,j] * cm_to_mm
        z = z_med[i,j] * cm_to_mm
        
        if j == Nj-1:
            node_numbers[i,j] = node_numbers[i,0]
        else:
            node_numbers[i,j] = node_number
            ops.node(node_number, x, y, z)
            node_number += 1
        #print(f"i = {i} j = {j} assigned={node_numbers[i,j]} ({x}, {y}, {z})")

# NUMBER OF ELEMENTS

Nele_i = Ni - 1
Nele_j = Nj - 1

Nelem = Nele_i * Nele_j

element_number = 1
for ei in range(Nele_i):
    for ej in range(Nele_j):
        ni = int(node_numbers[ei,ej])
        nj = int(node_numbers[ei+1,ej])
        nk = int(node_numbers[ei+1,ej+1])
        nl = int(node_numbers[ei,ej+1])
        # print(f"Element {element_number} ({ni} {nj} {nk} {nl}) ")

        # print(f"   x = {x_med[ei,ej]} {x_med[ei+1,ej]} {x_med[ei+1,ej+1]} {x_med[ei,ej+1]}")
        # print(f"   y = {y_med[ei,ej]} {y_med[ei+1,ej]} {y_med[ei+1,ej+1]} {y_med[ei,ej+1]}")
        # print(f"   z = {z_med[ei,ej]} {z_med[ei+1,ej]} {z_med[ei+1,ej+1]} {z_med[ei,ej+1]}")

        ops.element('ShellMITC4', element_number, ni, nj, nk, nl, section_tag)
        # ops.element('ShellDKGQ', element_number, ni, nj, nk, nl, section_tag)

        element_number += 1
        
# RECORD

ops.recorder("PVD","disp","disp")
ops.record()

# FIX NODES i = 0

for j in range(Nj-1):
    n = int(node_numbers[0, j])
    print(f"Fixing node # {n}")
    ops.fix(n, 1, 1, 1, 1, 1, 1)

# APPLY LOAD AT NODES i = Ni-1

tsTag =  1
ops.timeSeries('Linear', tsTag)

patternTag = 1
ops.pattern("Plain", patternTag, tsTag)

Ncargas = Nj-1
for j in range(Ncargas):
    n = int(node_numbers[Ni-1, j])
    print(f"Loading node # {n}")
    ops.load(n, 1.0/Ncargas, 0.0, 0.0, 0.0, 0.0, 0.0)

# CREATE SOE

ops.system("UmfPack")
#ops.system("ProfileSPD")

# CREATE DOF NUMBER

ops.numberer("RCM")

# CREATE CONSTRAIN HANDLER

ops.constraints("Plain")

# CREATE INTEGRATOR

ops.integrator("LoadControl", 1.0)

# CREATE ALGORITHM

ops.algorithm("Linear")

# CREATE ANALYSIS OBJECT

ops.analysis("Static")

# PERFORM THE ANALYSIS

ops.analyze(1)

# DISPLACEMENTS

ux = ops.nodeDisp(n,1)
uy = ops.nodeDisp(n,2)
uz = ops.nodeDisp(n,3)
print(f"Displacement at nodes n = {n} ({ux} {uy} {uz})")

# STIFFNESS

k = 1.0 / ux
print(f"Rigidez es k = {k}")

# CALCULATE VOLUME VARIATION

x_new = zeros(x_med.shape)
y_new = zeros(y_med.shape)
z_new = zeros(z_med.shape)

for i in range(Ni):
    for j in range(Nj):
        x_new[i,j] = x_med[i,j] + ops.nodeDisp(int(node_numbers[i,j]),1)
        y_new[i,j] = y_med[i,j] + ops.nodeDisp(int(node_numbers[i,j]),2)
        z_new[i,j] = z_med[i,j] + ops.nodeDisp(int(node_numbers[i,j]),3)

V_i = calcular_volumen(x_med, y_med, z_med)
V_f = calcular_volumen(x_new, y_new, z_new)
dV = V_i - V_f

print(f"dV = {dV}")

#fin = time.time()
#print(fin-inicio) # 1.0005340576171875
