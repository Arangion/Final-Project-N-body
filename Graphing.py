import matplotlib.pyplot as plt
import matplotlib.animation as animation

scatterX = []
scatterY = []
scatterZ = []

def updatePoints(position):
    # List of position arrays for graphing
    scatterX.append(position[:,0].flatten())
    scatterY.append(position[:,1].flatten())
    scatterZ.append(position[:,2].flatten())


def plotPoints3D(runTime):
    # Plotting Animation for 3D representation
    sim = plt.figure()
    ax = sim.add_subplot(111, projection="3d")
    sca, = ax.plot([],[],[],"o", markersize=2)
    ax.set_xlim3d(-10,10)
    ax.set_ylim3d(-10,10)
    ax.set_zlim3d(-10,10)
    def update_plot(i, x, y, z):
        # Updates xyz positions for each timestep
        sca.set_data(x[i], y[i])
        sca.set_3d_properties(z[i])
        return sca,
    # Runs the animation of the simulation from tStart to tEnd
    plot = animation.FuncAnimation(sim, update_plot, runTime, fargs=(scatterX, scatterY, scatterZ), interval=1000/30, blit=True)

    plt.show()

def plotPoints2D(runTime):
    # Plotting animation for 2D representation
    sim = plt.figure()
    ax = sim.subplots()
    sca, = ax.plot([],[], "o", markersize=2)
    ax.set_xlim(-1000,1000)
    ax.set_ylim(-1000,1000)


    def update_plot(i, x, y):
        # Updates xy position for each timestep
        sca.set_data(x[i], y[i])
        return sca,
    
    # Runs the animation of the simulation from tStart to tEnd
    plot = animation.FuncAnimation(sim, update_plot, runTime, fargs=(scatterX, scatterY), interval = 50, blit=True)

    plt.show()