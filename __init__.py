import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Particle:

    def __init__(self, id, mass, velocity, acceleration, position):
        self.id = str(id)
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.acceleration = acceleration
    
    def getPosition(self):
        return self.position

    def getVelocity(self):
        return self.velocity

    def getMass(self):
        return self.mass



class Simulation:
    def __init__(self, N, tStart, tEnd, tInterval, G):
        self.N = N
        self.time = tStart
        self.timeEnd = tEnd
        self.timeInc = tInterval
        self.Grav = G
        self.particles = {}
    
    def addParticle(self, particle):
            self.particles[particle.id] = particle
    
    def getParticle(self, id):
        return self.particles[id]

    def incrementTime(self):
        self.time += self.timeInc

    def updateDrift(self, acceleration, position):
        j = 0
        for i in self.particles.values():
            i.acceleration = acceleration[j]
            i.position = position[j]

    def updateKick(self, velocity):
        j = 0
        for i in self.particles.values():
            i.velocity = velocity[j]
            j += 1



        


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

    
#def getForce():
    
      


def main(N, tStart, tEnd, tInc, G, softening):
    # Creates the simulation object that will be run to record results of the simulation
    S = Simulation(N, tStart, tEnd, tInc, G)

    # Creates matricies of position[x,y,z], velocity[x,y,z], and mass[m] for each particle
    position = 2 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    print(velocity[0])
    mass = 10 * abs(np.random.randn(N,1))
    #mass = 10 * np.ones((N,1))/N -> Causes mass of particles to be 10, and each particle has the same mass


    acceleration = calcAcc(position, mass, G, softening)
    for i in range(N):
        S.addParticle(Particle(str(i), mass[i], velocity[i], acceleration[i], position[i]))
    

    sim = plt.figure()

    while(S.timeEnd > S.time):

        velocity += acceleration * tInc/2
        position += velocity * tInc
        acceleration = calcAcc(position, mass, G, softening)
        
        S.updateDrift(acceleration, position)

        velocity += acceleration * tInc/2

        S.updateKick(velocity)
        S.incrementTime()














main(1000, 0, 1000, 1, 6.674*(10**-11), 0.1)
