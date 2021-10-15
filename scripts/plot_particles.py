# import ast # parse list from file

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from mpl_toolkits import mplot3d

# fig = plt.figure()
# ax = plt.axes(projection='3d')
# # ax.set_xlim([0,100])
# # ax.set_ylim([0,100])
# # ax.set_zlim([0,100])

# for i in range(6):

#     position_data = pd.read_csv("frames/output_"+str(i)+".csv")
    
#     x = position_data.iloc[:,0][::1000]    
#     y = position_data.iloc[:,1][::1000]
#     z = position_data.iloc[:,2][::1000]
#     ax.scatter3D(x, y, z, c = z)

# plt.show()

import ast # parse list from file

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tables
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# from numba import vectorize

GRID_SIZE=300

frames_rendered = 7
x,y,z = [],[],[]
hfile = tables.open_file("positions.h5")
positions = np.array(hfile.root['positions'])
shp = positions.shape
print(shp)
for i in range(shp[0]):
    position_data = positions[i]
    x.append(position_data[:,0])   
    y.append(position_data[:,1])
    z.append(position_data[:,2])


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

#ax.set_xlim([-150,150])
#ax.set_ylim([-50,150])
#ax.set_zlim([-150,150])
title = ax.set_title('3D Test')
def get_cube():   
    phi = np.arange(1,10,2)*np.pi/4
    Phi, Theta = np.meshgrid(phi, phi)

    x = np.cos(Phi)*np.sin(Theta)
    y = np.sin(Phi)*np.sin(Theta)
    z = np.cos(Theta)/np.sqrt(2)
    return x,y,z


a = GRID_SIZE/4
b = GRID_SIZE/10
c = GRID_SIZE/4
cube_x,cube_y,cube_z = get_cube()


position_plots = ax.scatter(x[0],y[0],z[0],s=3,c=y[0])
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
def animate_vector(i):
    title.set_text("Frame: " + str(i))
    position_plots._offsets3d = (x[i], y[i], z[i])
    #position_plots._offsets3d = (x[i], y[i], z[i])
    #osition_plots.set_array(y[i])
        
myAnimation = FuncAnimation(fig, animate_vector, frames=np.arange(0,shp[0]), blit = False, interval=50)
#ax.plot_surface(cube_x*a+1, cube_y*b+b/2, cube_z*c+1,alpha=0.1)
hfile.close()
print("Saving animation")
myAnimation.save('R3_'+str(GRID_SIZE)+'-'+str(frames_rendered)+'.gif', fps=10000)