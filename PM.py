import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sp
import math
import PP

def particleMesh(N, position, velocity, acceleration, mass, G):
    # I decided to set the volume the simulation occurs in as a cube with length the longest
    # between x,y,z as particles leaving the system requires an extra case to handle, and I
    # don't want to deal with that yet so by having the particles more centered around the
    # middle it prevents that

    sideLength = max((position[:,0].max() - position[:,0].min()), position[:,1].max() - position[:,1].min(), position[:,2].max() - position[:,2].min())
    
    volume = sideLength**3

    negative = math.floor(-sideLength/2)
    positive = math.ceil(sideLength/2)

    xc, cellLength = np.linspace(negative, positive, num=int(sideLength), retstep = True) #These are positions of the cell-points x that are used to calculate density
    yc = np.linspace(negative, positive, num=int(sideLength)) # the cell length is also used for finding which cell a particle is closest to
    zc = np.linspace(negative, positive, num=int(sideLength))
    cellVolume = cellLength**3

    # creates an initial meshgrid of the positions and cell points,
    # creates two arrays of the same dimensions, N * len(xc)
    # the two are then subtracted to find which cell point is closest, and using
    # the index of that value the cell point is found

    cX, pX = np.meshgrid(xc, position[:,0])
    closestX = pX - cX
    chargesX = np.nonzero(np.piecewise(closestX, [(np.abs(closestX) <= cellLength/2), (np.abs(closestX) > cellLength/2)], [1, 0]))

    #distanceX = np.abs(closestX).min(axis = 1)
    #pointAssignments = np.argmin(np.abs(closestX), axis = 1)
    cY, pY = np.meshgrid(yc, position[:,1])
    closestY = pY - cY
    chargesY = np.nonzero(np.piecewise(closestY, [(np.abs(closestY) <= cellLength/2), (np.abs(closestY) > cellLength/2)], [1, 0]))

    #distanceY = np.abs(closestY).min(axis = 1)

    cZ, pZ = np.meshgrid(zc, position[:,2])
    closestZ = pZ - cZ
    chargesZ = np.nonzero(np.piecewise(closestZ, [(np.abs(closestZ) <= cellLength/2), (np.abs(closestZ) > cellLength/2)], [1, 0], [1, 0]))
    #distanceZ = np.abs(closestZ).min(axis = 1)
    print(closestZ)
    #print(closestX[20], cellLength/2)
    #print(chargesX[0])
    #print(chargesY)
    #print(chargesZ)


    xm, ym, zm = np.meshgrid(xc, yc, zc)

    
    #Ideally this is what I have to more or less finish
    #cellPotentials = Density*4*math.pi*G
    #print(cellPotentials)
    #gradient = sp.fft.fftn(cellPotentials)
    #print(gradient)



    #This is just a plotting for xm,ym,zm that I used for testing
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.scatter(xm,ym,zm)
    plt.show()

def testPM(N, tStart, tEnd, tInc, G, softening):
    #Position is an array of N particles with [:,0],[:,1],[:,2] being the x,y,z cordinates respectively
    #Velocity follows the same structure as position
    #mass is an 2-D array of N particles and their corresponding mass value
    #Acceleration is structured the same as position
    
    position = 1 + 2 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    mass = 10 * abs(np.random.randn(N,1))
    acceleration = PP.calcAcc(position, mass, G, softening)
    particleMesh(N, position, velocity, acceleration, mass, G)

# 100 particles is my preffered testing value, as more or less 
# makes the numbers either too small or too big for testing purposes
# Timestep has not been implemented yet for particleMesh
# G is gravitation constent don't change
# softening is a leftover from the Particle-Particle method, doesn't have any use yet 
# but I left it as it can be useful for fixing acceleration calculation
# issues from objects being too close to another 


#testPM(100, 0.0, 100, 0.1, 6.674*(10**-11), 0.1)
