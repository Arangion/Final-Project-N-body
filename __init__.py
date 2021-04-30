import numpy as np
import matplotlib.pyplot as plt


class Particle:

    def __init__(self, id, mass, velocity, position):
        self.id = str(id)
        self.mass = mass
        self.velocity = velocity
        self.position = position
        self.force = []
    
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


        


def calcAcc(position, mass, G):
    #print(position)
    print(mass)
    #print(position[:,0])
    #print(position[:,0:1])
    dX = (position[:,0] - position[:,0:1])
    dY = (position[:,1] - position[:,1:2])
    dZ = (position[:,2] - position[:,2:3])
    #m = (mass.T * mass)
    
    #print('mass', m)
    dR = (np.sqrt(dX**2 + dY**2 + dZ**2))**3
    #print('distance', dR)
    print((1/dR) @ mass)
    #print(G * (((1/dR) @ mass) * dX))

    #print(dX)
    #print(y)
    #print(z)

    
def getForce():
    pass
      


def main(N, tStart, tEnd, tInc, G):
    S = Simulation(N, tStart, tEnd, tInc, G)
    #position = np.random.randn(N,3)
    position = 2 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    mass = 10 * abs(np.random.randn(N))
    acceleration = calcAcc(position, mass, G)
    

    #for i in range(N):
        #S.addParticle(Particle(i, 5, velocity[i], position[i]))


    #print('Particle 4 positions', S.getParticle('4').getPosition())
    #print('Particle 4 position X', S.getParticle('4').getPosition()[0])
    #print('Particle 4 position X', S.getParticle('4').getPosition()[1])
    #print('Particle 4 position X', S.getParticle('4').getPosition()[2])
    #print(S.getParticle('4').mass)

main(5, 0, 1000, 1, 6.674*(10**-11))
