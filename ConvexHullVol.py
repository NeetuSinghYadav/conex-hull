from scipy.spatial import Voronoi,Delaunay, ConvexHull
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import MDAnalysis
import sys


u=MDAnalysis.Universe(sys.argv[1])

sel='name CA'

atoms=u.select_atoms(sel)

#print (atoms.positions)


CA_Coords=atoms.positions

#print (CA_Coords)

fig= plt.figure()
fig.set_size_inches(8,8)
ax=fig.add_subplot(111, projection= '3d')

hull = ConvexHull(CA_Coords)
ax.plot(CA_Coords[...,0], CA_Coords[...,1], CA_Coords [...,2],'o', label = 'CA_Coords')


for simplex in hull.simplices:
	ax.plot(CA_Coords[simplex, 0], CA_Coords[simplex, 1], CA_Coords[simplex, 2], 'k-')

ax.legend()

	
v = Voronoi(CA_Coords)	
print ("volume=", ConvexHull(CA_Coords).volume)	
		
def voronoi_volumes(CA_Coords):
    v = Voronoi(CA_Coords)
    vol = np.zeros(v.CA_Coords)
    for i, reg_num in enumerate(v.point_region):
        indices = v.regions[reg_num]
        if -1 in indices: # some regions can be opened
            vol[i] = np.inf
        else:
            vol[i] = ConvexHull(v.vertices[indices]).volume
    return vol

vol1=voronoi_volumes(CA_Coords)

for i in vol1:
	if i !='inf':
		print (i)




