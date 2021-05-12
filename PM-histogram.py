import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sp
import math


def particleMesh(N, position, velocity, acceleration, mass, G):
    # I decided to set the volume the simulation occurs in as a cube with length the longest
    # between x,y,z as particles leaving the system requires an extra case to handle, and I
    # don't want to deal with that yet so by having the particles more centered around the
    # middle it prevents that

    #sideLength = max((position[:,0].max() - position[:,0].min()), position[:,1].max() - position[:,1].min(), position[:,2].max() - position[:,2].min())
    
    # volume = sideLength**3

    # negative = math.floor(-sideLength/2)
    # positive = math.ceil(sideLength/2)

    # xc, cellLength = np.linspace(negative, positive, num=int(sideLength), retstep = True) #These are positions of the cell-points x that are used to calculate density
    # yc = np.linspace(negative, positive, num=int(sideLength)) # the cell length is also used for finding which cell a particle is closest to
    # zc = np.linspace(negative, positive, num=int(sideLength))
    # cellVolume = cellLength**3

    # # creates an initial meshgrid of the positions and cell points,
    # # creates two arrays of the same dimensions, N * len(xc)
    # # the two are then subtracted to find which cell point is closest, and using
    # # the index of that value the cell point is found

    # cX, pX = np.meshgrid(xc, position[:,0])
    
    # closestX = pX - cX
    # chargesX = np.nonzero(np.piecewise(closestX, [(np.abs(closestX) <= cellLength/2), (np.abs(closestX) > cellLength/2)], [1, 0]))
    # alt_chargesX = np.argmin(np.abs(closestX),1)

    # #distanceX = np.abs(closestX).min(axis = 1)
    # #pointAssignments = np.argmin(np.abs(closestX), axis = 1)
    # cY, pY = np.meshgrid(yc, position[:,1])
    # closestY = pY - cY
    # chargesY = np.nonzero(np.piecewise(closestY, [(np.abs(closestY) <= cellLength/2), (np.abs(closestY) > cellLength/2)], [1, 0]))
    # alt_chargesY = np.argmin(np.abs(closestY),1)

    # #distanceY = np.abs(closestY).min(axis = 1)

    # cZ, pZ = np.meshgrid(zc, position[:,2])
    # closestZ = pZ - cZ
    # chargesZ = np.nonzero(np.piecewise(closestZ, [(np.abs(closestZ) <= cellLength/2), (np.abs(closestZ) > cellLength/2)], [1, 0], [1, 0]))
    # alt_chargesZ = np.argmin(np.abs(closestZ),1)

    #distanceZ = np.abs(closestZ).min(axis = 1)
    
    # print(closestX[20], cellLength/2)
    # print(chargesX[0])
    # print(chargesY)
    # print(chargesZ)

    # This is where I am having an issue: if you look at the plot xm,ym,zm create it is a 3-d array of values
    # the issue is I don't conceptually understand how these work so that I could then add the mass of the particles
    # to the corresponding cell point which will allow for me to find its density so I can put it into the equation:
    # potential gravity = 4*G*pi*density
    # I would then get the potential gravity of each cell to create a gradient based on the values of the entire cube
    # Which would then allow me to interpolate those forces back onto the particles depending on their position

    #xm, ym, zm = np.meshgrid(xc, yc, zc)
    sideLength=12
    
    # Generate histograms with and without density function and weights (i.e, mass)
    #h, e = np.histogramdd(position,bins=int(sideLength),density=False)
    h, e = np.histogramdd(position,bins=int(sideLength),density=False,weights=mass.flatten())
    #h,e = np.histogramdd(position,bins=int(sideLength),density=True,weights=mass.flatten()*G) # Note: when done this way, the weights/densities are so small they're invisible?
    #h, e = np.histogramdd(position,bins=int(sideLength),density=True)

    # calculate the axes based on the bins of the histogram
    ex, ey, ez = e
    xm, ym, zm = np.meshgrid(np.linspace(ex[0], ex[-1], int(sideLength)),
                      np.linspace(ey[0], ey[-1], int(sideLength)),
                      np.linspace(ez[0], ez[-1], int(sideLength)))
    
    #Ideally this is what I have to more or less finish
    #cellPotentials = Density*4*math.pi*G
    #print(cellPotentials)
    #gradient = sp.fft.fftn(cellPotentials)
    #print(gradient)



    #This is just a plotting for xm,ym,zm that I used for testing

    # plot the 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')

    # use a scatterplot to display the particle densities at each bin location
    ax.scatter(xm.flatten(),ym.flatten(),zm.flatten(),s=h.flatten())
    plt.show()


def calcAcc(position, mass, G, softening):

    # dX,dY,dZ are the distances between each particle on their respective axes
    # dR is the equation accounting for the displacement between the particles with 
    # a softening value utilized as small values will return infinity without it
    # aX, aY, aZ is then the accelerations in each axis with gravity included
    # matmul uses the matrix of masses to do matrix multiplication on the displacements
    # to finish the acceleration equation

    dX = (position[:,0] - position[:,0:1])
    dY = (position[:,1] - position[:,1:2])
    dZ = (position[:,2] - position[:,2:3])

    dR = 1/((dX**2 + dY**2 + dZ**2 + softening**2)**1.5)

    aX = np.matmul(G * (dX * dR), mass)
    aY = np.matmul(G * (dY * dR), mass)
    aZ = np.matmul(G * (dZ * dR), mass)

    return(np.concatenate((aX, aY, aZ), axis=1))




def testPM(N, tStart, tEnd, tInc, G, softening):
    #Position is an array of N particles with [:,0],[:,1],[:,2] being the x,y,z cordinates respectively
    #Velocity follows the same structure as position
    #mass is an 2-D array of N particles and their corresponding mass value
    #Acceleration is structured the same as position
    
    position = 3 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    mass = abs(np.random.randn(N,1))
    acceleration = calcAcc(position, mass, G, softening)
    particleMesh(N, position, velocity, acceleration, mass, G)

# 100 particles is my preffered testing value, as more or less 
# makes the numbers either too small or too big for testing purposes
# Timestep has not been implemented yet for particleMesh
# G is gravitation constent don't change
# softening is a leftover from the Particle-Particle method, doesn't have any use yet 
# but I left it as it can be useful for fixing acceleration calculation
# issues from objects being too close to another 
testPM(10000, 0.0, 100, 0.1, 6.674*(10**-11), 0.1)
