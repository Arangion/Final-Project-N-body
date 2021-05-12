import numpy as np
import Graphing
import pandas as pd

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

def particleParticle(position, velocity, acceleration, mass, tInc, G, softening, tStart, tEnd):
    
    Graphing.updatePoints(position)
    t = tStart
    num = 0
    p1 = position  
    weird = 0
    while(t < tEnd):    

        velocity += (acceleration * (tInc/2.0))

        position += (velocity * tInc)
        
        acceleration = calcAcc(position, mass, G, softening)
        
        velocity += (acceleration * (tInc/2.0))

        t += tInc

        Graphing.updatePoints(position)

   # Graphing.plotPoints3D(int(tEnd/tInc))
    Graphing.plotPoints2D(int(tEnd/tInc))



