import numpy as np
import matplotlib.pyplot as plt
from matplotlib.text import Text
from mpl_toolkits.mplot3d import Axes3D
from numpy.random import rand

print_lables=True
#print_lables=False

F = False
T = True

nsd = 2
lplot_global_coords = T
lplot_low_order_only = T
plot_edge_nodes = T
plot_face_nodes = T
if nsd == 2:
    plot_vol_nodes  = F
else:
    plot_vol_nodes  = T


    
#
# USER: DO NO TOUCH from here on!
#
plt.close('all')
fig = plt.figure()
ax3d = fig.add_subplot(projection='3d')
ax3d.set_box_aspect(aspect = (1,1,1))


if (lplot_global_coords == True):
    coords = np.loadtxt('COORDS_GLOBAL.dat', usecols=range(4))

    x=coords[:,0];
    y=coords[:,1];
    z=coords[:,2];
    ip=coords[:,3];

    scatter = ax3d.scatter(x, y, z, marker='x',  picker=True)
    lplot_low_order_only = F
    plot_edge_nodes = F
    plot_face_nodes = F
    
    if (print_lables == True):
        for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
            ax3d.text(xcoords, ycoords, zcoords, int(label))

#
# Load coords
#
# 1) Low order:
if (lplot_low_order_only == True):
    coords_lo = np.loadtxt('COORDS_LO.dat', usecols=range(4))
    
    x=coords_lo[:,0];
    y=coords_lo[:,1];
    z=coords_lo[:,2];
    ip=coords_lo[:,3];
    
    scatter = ax3d.scatter(x, y, z, marker='o',  picker=True)
    
    if (print_lables == True):
        for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
            ax3d.text(xcoords, ycoords, zcoords, int(label))

            
if (plot_edge_nodes == True):
    # 2) high order edges
    del x, y, z, ip
    coords_ho = np.loadtxt('COORDS_HO_edges.dat', usecols=range(4))
    
    x=coords_ho[:,0];
    y=coords_ho[:,1];
    z=coords_ho[:,2];
    ip=coords_ho[:,3];
    
    #ax3d.scatter(x, y, z, marker='x')
    line = ax3d.scatter(x, y, z, marker='x', picker=True, pickradius=5)  # 5 points tolerance
    
    if (print_lables == True):
        for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
            ax3d.text(xcoords, ycoords, zcoords, int(label))
        

if (plot_face_nodes == True):
    # 3) high order faces
    del x, y, z, ip
    coords_ho = np.loadtxt('COORDS_HO_faces.dat', usecols=range(4))
    
    x=coords_ho[:,0];
    y=coords_ho[:,1];
    z=coords_ho[:,2];
    ip=coords_ho[:,3];
    
    ax3d.scatter(x, y, z, marker='s')
    
    if (print_lables == True):
        for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
            ax3d.text(xcoords, ycoords, zcoords, int(label))
        

if (plot_vol_nodes == True):
    # 4 high order internal
    del x, y, z, ip
    coords_ho = np.loadtxt('COORDS_HO_vol.dat', usecols=range(4))
    
    x=coords_ho[:,0];
    y=coords_ho[:,1];
    z=coords_ho[:,2];
    ip=coords_ho[:,3];
    
    ax3d.scatter(x, y, z, marker='s')
    
    if (print_lables == True):
        for xcoords, ycoords, zcoords, label in zip(x, y, z, ip):
            ax3d.text(xcoords, ycoords, zcoords, int(label))

# Make axes limits 
my_aspect_ratio = min(max(x)/max(z), max(x)/max(y))
ax3d.set_box_aspect((my_aspect_ratio, 1, 1))

xmax = max(x)
xmin = min(x)
eps = 0.1*(xmax - xmin)
ax3d.axes.set_xlim3d(left=xmin-eps, right=xmax+eps) 
#ax3d.axes.set_ylim3d(bottom=min(y)-eps, top=max(y)+eps)
#ax3d.axes.set_zlim3d(bottom=min(z)-eps, top=max(z)+eps)

plt.show()
