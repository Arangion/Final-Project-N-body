import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import scipy.interpolate as ip
import math
import Graphing

def _power_shift(k, n):
    # power_shift calcuates the wavenumber value to be used in kernel
    save = np.seterr(divide = 'ignore')
    a = np.where(k == 0, 0, k**n)
    np.seterr(**save)

    return a

def wavenumber(v):
    # Caluclates wavenumber based on number of particles
    N = v[0]
    i = np.indices(v)

    return(np.where(i > N/2, i-N, i))

def potential_integration(x, y, z, F, g, bL):
    # calculates the x,y,z difference between potential for gradient
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
    # del_x, del_y, and del_z find the distance between the gridpoint and the particle
    del_x = lambda i: xBins[int(g[i,0])] - position[i,0]
    del_y = lambda i: yBins[int(g[i,1])] - position[i,1]
    del_z = lambda i: zBins[int(g[i,2])] - position[i,2]
    
    pA = np.zeros((N,3))
    # Each k value is the positional distances between the grid points and particle
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
            # potential_calc_ calculates the derivative of the grid distance 
            # to find the gradient of gravitation potential on the particle
            pA[i] = potential_calc_rollu(F, bL, g, i, kv, N)
        else:
            pA[i] = potential_calc_(F, bL, g, i, kv)
    return(pA)


def allocateDensity(position, sideLength, N):
    dr = position - np.floor(position)
    # Generate histograms of particels and gridpoints to determine density based on distance to particle
    # Follows Cloud in Cloud method, which allocates to the nearest 8 gridpoints in a 3d space
    
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

def calculateMesh(N, position, velocity, mass, G):

    #Set side length of grid to 10 units
    sideLength=10

    #Allocate Denisty based on particle distances to each gridpoint in mesh
    density = allocateDensity(position, sideLength, N)
    #add the mass to the gridpoint position based on desnity
    density *= mass[0]
    density -= 1.0

    #Mesh position histogram, e has the distances between each gridpoint xyz
    h, e = np.histogramdd(position, bins=int(sideLength))
    ex, ey, ez = e

    #A Fast Fournier Transform is used to move the density field to the frequency domain
    PE = np.fft.fftn(density)
    # Calculates the wave of the simulation, waveshape determines wavelength xyz
    waveshape = (sideLength,) *3
    wave = (wavenumber(waveshape)) * 2*math.pi/sideLength
    
    # kernel is the reciprocal wavenumer squared
    kernel = -(_power_shift((wave**2).sum(axis=0), -1.))
    
    # Force is then caluclated from the real Inverse Fast Fournier Transform allowing for Gravitation Constant
    F = np.fft.ifftn(PE*kernel).real * G/sideLength
    
    #Grav is an array for indexing the position in space of each particle, xbins,ybins,zbins are arrays of the bin lengths xyz
    Grav = np.zeros((N,3))
    xBins = np.asarray(ex)
    yBins = np.asarray(ey)
    zBins = np.asarray(ez)
    # bL is the length between grid points
    bL = ex[1] - ex[0]
    # Allocates potentials of Grid to positions in grid
    for i in range(N):
        Grav[i] = (np.argmax(xBins >= position[i,0]) - 1), (np.argmax(yBins >= position[i,1]) - 1), (np.argmax(zBins >= position[i,2]) - 1)
    # Interpolate calucaltes the force on each particle based on the distance from mesh grid points
    pG = interpolate(F, N, bL, xBins, yBins, zBins, Grav, position)

    return(pG)

def ParticleMesh(N, position, velocity, mass, tInc, G, softening, tStart, tEnd):

    acceleration = calculateMesh(N, position, velocity, mass, G)

    Graphing.updatePoints(position)
    
    t = tStart
    # Time step integration of particles to simulate movement
    while(tEnd > t):
        velocity += (acceleration * (tInc/2.0))
        
        position += (velocity * tInc)
    
        acceleration = calculateMesh(N, position, velocity, mass, G)

        velocity += (acceleration * (tInc/2.0))

        t += tInc

        Graphing.updatePoints(position)

    Graphing.plotPoints3D(int(tEnd/tInc)) # Run a 3D representation of the simulation
    #Graphing.plotPoints2D(int(tEnd/tInc)) # Run a 2D representation of the simulation