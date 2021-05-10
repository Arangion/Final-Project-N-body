from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy as sp
import math



class Particles:

    def __init__(self, id, position, velocity, acceleration):
        self.id = str(id)
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
    
    def addParticles(self, particle):
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

class Cell:
    def _init__(self, id, particles):
        self.id = id
        self.particles = particles
    
    def addParticle(self, particle):
        self.particles.append(particle)

    def removeParticle(self, particle):
        self.particles.remove(particle)
        

class Mesh:
    def __init__(self, meshDimensions, pNum):
        self.volume = self.findVolume(meshDimensions)
        self.cellCount = pNum//10
        cellDimensions = [self.x/cellCount, self.y/cellCount, self.z/cellCount]
        self.cells = self.genCells(cellDimensions, cellCount)

    def findVolume(self, meshDimensions):
        self.x = meshDimensions[0]
        self.y = meshDimensions[1]
        self.z = meshDimensions[2]
        return(x * y * z)
    
    def genCells(self, cellDimensions, n):
        temp = []
        for i in range(n):
            temp.append([cellDimensions[0], cellDimensions[1], cellDimensions[2] ])

        


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

    
def particleParticle(position, velocity, acceleration, mass, tInc, G, softening):
        p = position
        v = velocity
        a = acceleration
        
        v += acceleration * tInc/2
        p += v * tInc
        a = calcAcc(p, mass, G, softening)
        v += a * tInc/2

        return(p, v, a)

def particleMesh(N, position, velocity, acceleration, mass, G):

    sideLength = max((position[:,0].max() - position[:,0].min()), position[:,1].max() - position[:,1].min(), position[:,2].max() - position[:,2].min())
    
    volume = sideLength**3
    
    
    #print(volume, '=', sideLength)
    
    cellSide = math.ceil((N/5)**(1/3))
    
    cellTot = cellSide**3

    print(cellTot)
    print(cellSide)
    print(cellTot//cellSide)
    xTmp = 0
    yTmp = 0
    xCounter = 0
    cellPoints = np.zeros((cellTot, cellSide))
    #print(cellPoints)
    print(cellPoints[0])
    cellLength = (sideLength/3)
    cellVolume = (cellLength**3)
    midPoint = (cellLength/2)
    #print(midPoint, cellLength, sideLength, cellSide, cellTot)
    for i in range(0, cellTot, cellTot//cellSide):
        while(yTmp < cellSide):
            cellPoints[i+(cellSide*yTmp)] = midPoint + (cellLength*(xCounter)), midPoint + (cellLength*(yTmp)), midPoint
            cellPoints[i+1+(cellSide*yTmp)] = midPoint + (cellLength*(xCounter)), midPoint + (cellLength*(yTmp)), midPoint + cellLength
            cellPoints[i+2+(cellSide*yTmp)] = midPoint + (cellLength*(xCounter)), midPoint + (cellLength*(yTmp)), midPoint + (cellLength*2)
            if(xTmp == (cellSide-1)):
                xTmp = 0
                yTmp += 1
            else:
                xTmp += 1
        print(i)
        xCounter += 1
    
    print(cellPoints)
    print(N, cellTot)
    tempMatrix = np.zeros((N, cellTot))
    closestCell = np.zeros((N, 2))
    print(tempMatrix)

    for i in range(N):
        tempMatrix[i] = np.sqrt((((position[i,0] - cellPoints[:,0])**2) + ((position[i,1] - cellPoints[:,1])**2) + ((position[i,1] - cellPoints[:,1])**2)))
    print(tempMatrix)
    
    closestCell[:,0] = np.amin(tempMatrix, keepdims = True) 
    closestCell[:,1] = np.argmin(tempMatrix, axis = 1)
    
    Density = np.zeros(cellTot)

    for i in range(cellTot):
        Density[i] +=  mass[np.argwhere(i == closestCell[1])[1]]
        Density[i] = Density[i]/cellVolume
    
    print(closestCell)

    cellPotentials = Density*4*math.pi*G
    gradient = sp.fft.fftn(cellPotentials)
    return gradient
        









    
def testPM(N, tStart, tEnd, tInc, G, softening):
    position = 2 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    mass = 10 * abs(np.random.randn(N,1))
    acceleration = calcAcc(position, mass, G, softening)
    particleMesh(N, position, velocity, acceleration, mass, G)


    


def particleParticleMesh():
    pass
    
def animatePlot(n, p, scatters, N):
    for i, j in zip(scatters, p[str(n)].position):
        for k in range(N):
            i._offsets3d(j[k,0], j[k,1],j[k,2])
    return scatters



def main(N, tStart, tEnd, tInc, G, softening):
    # Creates the simulation object that will be run to record results of the simulation
    S = Simulation(N, tStart, tEnd, tInc, G)

    # Creates matricies of position[x,y,z], velocity[x,y,z], and mass[m] for each particle
    position = 2 * np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    #print(velocity[0])
    mass = 10 * abs(np.random.randn(N,1))
    #mass = 10 * np.ones((N,1))/N -> Causes mass of particles to be 10, and each particle has the same mass

    acceleration = calcAcc(position, mass, G, softening)

    S.addParticles(Particles(0, position, velocity, acceleration))

    num = 1
    #print(S.timeEnd, S.time)
    while(S.timeEnd > S.time):
        position, velocity, acceleration = particleParticle(position, velocity, acceleration, mass, tInc, G, softening)        
        #print(position, '\n', velocity, '\n', acceleration)
        S.addParticles(Particles(num, position, velocity, acceleration))
        num += 1
        S.incrementTime()

    #print(S.particles.keys())
    sim = plt.figure()
    ax = sim.add_subplot(projection="3d")
    #for i in range(int(tEnd/tInc)):
        #ax.scatter(S.particles[str(i)].position[:,0], S.particles[str(i)].position[:,1], S.particles[str(i)].position[:,1])
    scatters = [ ax.scatter(S.particles[str(i)].position[:,0], S.particles[str(i)].position[:,1], S.particles[str(i)].position[:,2]) for i in range(int(tEnd/tInc))]
    print(S.particles['0'].position)
    #print(scatters)
    #ax.set_xlim3d([-100,100])
    #ax.set_ylim3d([-100,100])
    #ax.set_zlim3d([-100,100])
    plot = animation.FuncAnimation(sim, animatePlot, int(tEnd/tInc), fargs=(S.particles, scatters, N), interval=50)
    plt.show()
    


    #print(position)
    #print(velocity)





testPM(50, 0.0, 100, 0.1, 6.674*(10**-11), 0.1)
#main(10, 0.0, 100, 0.1, 6.674*(10**-11), 0.1)
