import pylab as plt
import numpy as np
import csv
import os
import matplotlib.animation as animation
#plt.rcParams['animation.ffmpeg_path'] = r"C:\some_path\ffmpeg.exe"   # if necessary

# Generate data for plotting
Lx = Ly = 5e-5
Nx = Ny = 100
Nt = 475
x = np.linspace(-Lx, Lx, Nx)
y = np.linspace(-Lx, Lx, Nx)

array = np.ndarray((Nx, Ny))


def some_data(i):   # function returns a 2D data array
    print(i)
    csvfile = open("output/dataout%d.csv" % (i,))
    reader = csv.reader(csvfile, dialect='excel')
    for i, row in enumerate(reader):
        for j, value in enumerate(row):
            array[i][j] = float(value)
    csvfile.close()
    return array
    

fig = plt.figure()
extent = np.linspace(-1e-4, 1e-4, 100)
cont = plt.contourf(x, y, some_data(0), levels=extent)    # first image on screen
plt.colorbar()


# animation function
def animate(i):
    global cont
    z = some_data(i)
    for c in cont.collections:
        c.remove()
    cont = plt.contourf(x, y, z,levels=extent)
    return cont


writervideo = animation.PillowWriter(fps=60)
anim = animation.FuncAnimation(fig, animate, frames=Nt, repeat=False)
anim.save("animation.gif", writervideo)













"""


from np.core.fromnumeric import shape

minx, miny, maxx, maxy = -5e-5, -5e-5, 5e-5, 5e-5
numcells = 11
fig = plt.figure()
xdata, ydata = np.linspace(minx, maxx, numcells),\
               np.linspace(miny, maxy, numcells)
datalist = []
with open("output/testex%d.csv" % (0,)) as csvfile:
    reader = csv.reader(csvfile, dialect='excel')
    for row in reader:
        datalist.append([])
        for value in row:
            datalist[-1].append(float(value))
array0 = np.array(datalist)
print(shape(array0))
ln = plt.contourf(xdata, ydata, array0)


def init():
    ln.set_xlim(-5e-5, 5e-5)
    ln.set_ylim(-5e-5, 5e-5)
    return ln,


def update(frame):
    data = []
    with open("output/testex%d.csv" % (frame,)) as csvfile:
        reader = csv.reader(csvfile, dialect='excel')
        for row in reader:

            data.append([])
            for value in row:
                data[-1].append(float(value))
    array = np.array(data)
    ln.set_array(array)
    ln.show()
    return ln


ani = FuncAnimation(fig, update, frames=range(100),
                    init_func=init, blit=True)
plt.show()
"""