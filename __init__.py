import numpy as np
from PP import particleParticle as PP
from PM import ParticleMesh as PM

def main(N, m, f, tStart, tEnd, tInc, G, softening):
    # Creates matricies of position[x,y,z], velocity[x,y,z], and mass[m] for each particle
    # Position and velocity are arrays of N particles with [:,0],[:,1],[:,2] being the x,y,z cordinates respectively
    # mass is an 2-D array of N particles and their corresponding mass value, 
    position = np.random.randn(N,3)
    
    velocity = np.random.randn(N,3)
    
    mass = m * np.ones((N,1))/N #-> Causes total mass of particles to be 10, and each particle has the same mass
    
    f (N, position, velocity, mass, tInc, G, softening, tStart, tEnd)


#Inputs of Main: 
#   N           =   Number of Particles
#   m           =   Total Mass of System
#   f           =   Particle Function (PP or PM)
#   tStart      =   Start time of Simulation
#   tEnd        =   End time of Simulation
#   tInc        =   Timestep of simulation
#   G           =   Gravitational Constant
#   softening   =   Softening Factor(0.1 or 0.01)

main(100, 10, PM, 0.0, 10.0, 0.01, 1.0, 0.01)
