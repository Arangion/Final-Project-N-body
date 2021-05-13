import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.interpolate as ip
import math
import Graphing

def _power_shift(k, n):
    save = np.seterr(divide = 'ignore')
    a = np.where(k == 0, 0, k**n)
    np.seterr(**save)

    return a

def wavenumber(v):
    N = v[0]
    i = np.indices(v)

    return(np.where(i > N/2, i-N, i))

def potential_integration(x, y, z, F, g, bL):
    deltaX = F[int(g[x[1],0]), int(g[y[0],1]), int(g[z[0],2])] - F[int(g[x[2],0]), int(g[y[0],1]), int(g[z[0],2])]
    deltaY = F[int(g[x[0],0]), int(g[y[1],1]), int(g[z[0],2])] - F[int(g[x[0],0]), int(g[y[2],1]), int(g[z[0],2])]
    deltaZ = F[int(g[x[0],0]), int(g[y[0],1]), int(g[z[1],2])] - F[int(g[x[0],0]), int(g[y[0],1]), int(g[z[2],2])]

    return (np.array([deltaX/(2*bL), deltaY/(2*bL), deltaZ/(2*bL)]))


def potential_calc_(F, bL, g, i, kv):
    
    gxyz = potential_integration([i, i+1, i-1], [i, i+1, i-1], [i, i+1, i-1], F, g, bL)
    
    gxyz_ = potential_integration([i, i+1, i-1], [i, i+1, i-1], [i+1, i+2, i], F, g, bL)

    gxy_z = potential_integration([i, i+1, i-1], [i+1, i+2, i], [i, i+1, i-1], F, g, bL)
    
    gxy_z_ = potential_integration([i, i+1, i-1], [i+1, i+2, i], [i+1, i+2, i], F, g, bL)
    
    gx_yz = potential_integration([i+1, i+2, i], [i, i+1, i-1], [i, i+1, i-1], F, g, bL)
    
    gx_yz_ = potential_integration([i+1, i+2, i], [i, i+1, i-1], [i+1, i+2, i], F, g, bL)

    gx_y_z = potential_integration([i+1, i+2, i], [i+1, i+2, i], [i, i+1, i-1], F, g, bL)

    gx_y_z_ = potential_integration([i+1, i+2, i], [i+1, i+2, i], [i+1, i+2, i], F, g, bL)

    gp = kv[0]*gxyz + kv[1]*gxyz_ + kv[2]*gxy_z + kv[3]*gxy_z_ + kv[4]*gx_yz + kv[5]*gx_yz_ + kv[6]*gx_y_z + kv[7]*gx_y_z_

    return(gp)

def potential_calc_rollu(F, bL, g, i, kv, N):
    if(i == N-1):

        gxyz = potential_integration([i, 0, i-1], [i, 0, i-1], [i, 0, i-1], F, g, bL)
    
        gxyz_ = potential_integration([i, 0, i-1], [i, 0, i-1], [0, 1, i], F, g, bL)
    
        gxy_z = potential_integration([i, 0, i-1], [0, 1, i], [i, 0, i-1], F, g, bL)
    
        gxy_z_ = potential_integration([i, 0, i-1], [0, 1, i], [0, 1, i], F, g, bL)

        gx_yz = potential_integration([0, 1, i], [i, 0, i-1], [i, 0, i-1], F, g, bL)

        gx_yz_ = potential_integration([0, 1, i], [i, 0, i-1], [0, 1, i], F, g, bL)

        gx_y_z = potential_integration([0, 1, i], [0, 1, i], [i, 0, i-1], F, g, bL)

        gx_y_z_ = potential_integration([0, 1, i], [0, 1, i], [0, 1, i], F, g, bL)

    if(i == N-2):

        gxyz = potential_integration([i, i+1, i-1], [i, i+1, i-1], [i, i+1, i-1], F, g, bL)
    
        gxyz_ = potential_integration([i, i+1, i-1], [i, i+1, i-1], [i+1, 0, i], F, g, bL)

        gxy_z = potential_integration([i, i+1, i-1], [i+1, 0, i], [i, i+1, i-1], F, g, bL)
    
        gxy_z_ = potential_integration([i, i+1, i-1], [i+1, 0, i], [i+1, 0, i], F, g, bL)
    
        gx_yz = potential_integration([i+1, 0, i], [i, i+1, i-1], [i, i+1, i-1], F, g, bL)
    
        gx_yz_ = potential_integration([i+1, 0, i], [i, i+1, i-1], [i+1, 0, i], F, g, bL)

        gx_y_z = potential_integration([i+1, 0, i], [i+1, 0, i], [i, i+1, i-1], F, g, bL)

        gx_y_z_ = potential_integration([i+1, 0, i], [i+1, 0, i], [i+1, 0, i], F, g, bL)

    gp = (kv[0]*gxyz + kv[1]*gxyz_ + kv[2]*gxy_z + kv[3]*gxy_z_ + kv[4]*gx_yz + kv[5]*gx_yz_ + kv[6]*gx_y_z + kv[7]*gx_y_z_)*-1

    return(gp)


def interpolate(F, N, bL, xBins, yBins, zBins, g, position):
    del_x = lambda i: xBins[int(g[i,0])] - position[i,0]
    del_y = lambda i: yBins[int(g[i,1])] - position[i,1]
    del_z = lambda i: zBins[int(g[i,2])] - position[i,2]
    
    pA = np.zeros((N,3))
    for i in range(N):
        kv = []
        kv.append((bL-del_x(i))*(bL-del_y(i))*(bL-del_z(i)))
        kv.append((bL-del_x(i))*(bL-del_y(i))*(del_z(i)))
        kv.append((bL-del_x(i))*(del_y(i))*(bL-del_z(i)))
        kv.append((bL-del_x(i))*(del_y(i))*(del_z(i)))
        kv.append((del_x(i))*(bL-del_y(i))*(bL-del_z(i)))
        kv.append((del_x(i))*(bL-del_y(i))*(del_z(i)))
        kv.append((del_x(i))*(del_y(i))*(bL-del_z(i)))
        kv.append((del_x(i))*(del_y(i))*(del_z(i)))

        if(i >= N-2):
            pA[i] = potential_calc_rollu(F, bL, g, i, kv, N)
        else:
            pA[i] = potential_calc_(F, bL, g, i, kv)
    return(pA)


def allocateDensity(position, sideLength, N):
    dr = position - np.floor(position)
    # Generate histograms with and without density function and weights (i.e, mass)
    
    CIC = np.zeros((sideLength,sideLength,sideLength))
    h, e = np.histogramdd(position, bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(1-dr[:,0]) * (1-dr[:,1]) * (1-dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(dr[:,0]) * (1-dr[:,1]) * (1-dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(dr[:,0]) * (dr[:,1]) * (1-dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(dr[:,0]) * (1-dr[:,1]) * (dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(1-dr[:,0]) * (dr[:,1]) * (1-dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(1-dr[:,0]) * (dr[:,1]) * (dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(1-dr[:,0]) * (1-dr[:,1]) * (dr[:,2]))
    CIC += h
    h, e = np.histogramdd(position,bins=int(sideLength), range = [[0, N], [0,N], [0,N]], weights=(dr[:,0]) * (dr[:,1]) * (dr[:,2]))
    CIC += h

    return CIC

def particleMesh(N, position, velocity, mass, G):
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

    # This is where I am having an issue: if you look at the plot xm,ym,zm create it is a 3-d array of values
    # the issue is I don't conceptually understand how these work so that I could then add the mass of the particles
    # to the corresponding cell point which will allow for me to find its density so I can put it into the equation:
    # potential gravity = 4*G*pi*density
    # I would then get the potential gravity of each cell to create a gradient based on the values of the entire cube
    # Which would then allow me to interpolate those forces back onto the particles depending on their position

    #xm, ym, zm = np.meshgrid(xc, yc, zc)
    sideLength=10
    
    
    #h, e = np.histogramdd(position,bins=int(sideLength),density=False,weights=mass.flatten())
    #h,e = np.histogramdd(position,bins=int(sideLength),density=True,weights=mass.flatten()) # Note: when done this way, the weights/densities are so small they're invisible?
    #h1, e1 = np.histogramdd(position,bins=int(sideLength),density=True)

    # calculate the axes based on the bins of the histogram
    density = allocateDensity(position, sideLength, N)
    density *= mass[0]
    density -= 1.0
    h, e = np.histogramdd(position, bins=int(sideLength))
    ex, ey, ez = e
    xm, ym, zm = np.meshgrid(np.linspace(ex[0], ex[-1], int(sideLength)),
                     np.linspace(ey[0], ey[-1], int(sideLength)),
                     np.linspace(ez[0], ez[-1], int(sideLength)))

    #print(xm.shape)
    #print(ym.shape)
    #print(zm.shape)
    PE = np.fft.fftn(density)
    waveshape = (sideLength,) *3
    wave = (wavenumber(waveshape)) * 2*math.pi/sideLength
    
    kernel = -(_power_shift((wave**2).sum(axis=0), -1.))
    #print(PE.shape)
    #print(kernel.shape)
    #int(sideLength)), np.linspace(PE, kernel[1], int(sideLength)), np.linspace(PE, kernel[2], int(sideLength)))
    #print(PEx)
    #print(Kx)
    
    F = np.fft.ifftn(PE*kernel).real * G/sideLength
    
    Grav = np.zeros((N,3))
    xBins = np.asarray(ex)
    yBins = np.asarray(ey)
    zBins = np.asarray(ez)
    bL = ex[1] - ex[0]
    for i in range(N):
        Grav[i] = (np.argmax(xBins >= position[i,0]) - 1), (np.argmax(yBins >= position[i,1]) - 1), (np.argmax(zBins >= position[i,2]) - 1)
    pG = interpolate(F, N, bL, xBins, yBins, zBins, Grav, position)
    #print(pG.shape)
    return(pG)

def testPM(N, tStart, tEnd, tInc, G):
    #Position is an array of N particles with [:,0],[:,1],[:,2] being the x,y,z cordinates respectively
    #Velocity follows the same structure as position
    #mass is an 2-D array of N particles and their corresponding mass value
    #Acceleration is structured the same as position
    
    position = np.random.randn(N,3)
    velocity = np.random.randn(N,3)
    #mass = abs(np.random.randn(N,1))
    mass = 10 * np.ones((N,1))/N
    #acceleration = calcAcc(position, mass, G, softening)
    acceleration = particleMesh(N, position, velocity, mass, G)
    t = tStart
    while(tEnd > t):
        velocity += (acceleration * (tInc/2.0))
        
        position += (velocity * tInc)
    
        acceleration = particleMesh(N, position, velocity, mass, G)

        velocity += (acceleration * (tInc/2.0))

        t += tInc

        Graphing.updatePoints(position)

    Graphing.plotPoints3D(int(tEnd/tInc))

# 100 particles is my preffered testing value, as more or less 
# makes the numbers either too small or too big for testing purposes
# Timestep has not been implemented yet for particleMesh
# G is gravitation constent don't change
# softening is a leftover from the Particle-Particle method, doesn't have any use yet 
# but I left it as it can be useful for fixing acceleration calculation
# issues from objects being too close to another 
testPM(1000, 0.0, 10, 0.01, 1.0)
