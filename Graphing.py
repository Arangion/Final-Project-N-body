import matplotlib.pyplot as plt
import matplotlib.animation as animation

scatterX = []
scatterY = []
scatterZ = []

def updatePoints(position):

    scatterX.append(position[:,0].flatten())
    scatterY.append(position[:,1].flatten())
    scatterZ.append(position[:,2].flatten())


def plotPoints3D(runTime):
    sim = plt.figure()
    ax = sim.add_subplot(111, projection="3d")
    sca, = ax.plot([],[],[],"o", markersize=2)
    ax.set_xlim3d(-20,20)
    ax.set_ylim3d(-20,20)
    ax.set_zlim3d(-20,20)
    def update_plot(i, x, y, z):
        sca.set_data(x[i], y[i])
        sca.set_3d_properties(z[i])
        return sca,

    plot = animation.FuncAnimation(sim, update_plot, runTime, fargs=(scatterX, scatterY, scatterZ), interval=1000/30, blit=True)

    plt.show()

def plotPoints2D(runTime):
    sim = plt.figure()
    ax = plt.subplots()
    sca, = ax.plot([],[], "o", markersize=2)

    def update_plot(i, x, y):
        sca.set_data(x[i], y[i])
        return sca,
    
    plot = animation.FuncAnimation(sim, update_plot, runTime, fargs=(scatterX, scatterY), interval = 50, blit=True)

    plt.show()