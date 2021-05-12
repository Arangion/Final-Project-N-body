from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sp
import math
import PP
import PM

def main(N, tStart, tEnd, tInc, G, softening):
    # Creates matricies of position[x,y,z], velocity[x,y,z], and mass[m] for each particle
    position = np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    #mass =  abs(np.random.randn(N,1))
    mass = 10 * np.ones((N,1))/N #-> Causes mass of particles to be 10, and each particle has the same mass
    acceleration = PP.calcAcc(position, mass, G, softening)
    
    PP.particleParticle(position,velocity, acceleration, mass, tInc, G, softening, tStart, tEnd)


main(100, 0.0, 10.0, 0.01, 6.674*(10**-11), 0.01)
