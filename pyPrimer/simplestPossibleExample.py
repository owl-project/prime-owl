#!/usr/bin/python3

# This example is supposed to serve as "simplest possible example" for how to
# use pypr. It creates a super-simple scene consisting of a single 'box' (made
# up of 12 triangles), then traces a set of NxN parallel rays against this box,
# and prints the results. This mostly serves to illustrate how to create a
# model of triangle mesh[es], how to set up arrays of rays and hits to be used
# in trace(), and how to interpret the results.

#import pypr
import math
import numpy as np
import pypr

float_t=np.float32
int_t=np.int32

# a hard-coded box shape: vertex array
vertices = np.array([
    [-1.,-1.,-1.],
    [+1.,-1.,-1.],
    [-1.,+1.,-1.],
    [+1.,+1.,-1.],
    [-1.,-1.,+1.],
    [+1.,-1.,+1.],
    [-1.,+1.,+1.],
    [+1.,+1.,+1.]
    ], dtype='float32')

# a hard-coded box shape: index array
indices = np.array([
     [0,1,3],  [2,3,0],
     [5,7,6],  [5,6,4],
     [0,4,5],  [0,5,1],
     [2,3,7],  [2,7,6],
     [1,5,7],  [1,7,3],
     [4,0,2],  [4,2,6]
     ], dtype='int32')

def main():
    print("#pyPrimer: creating context...")
    pypr.context_create()
    
    # -----------------------------------------------------------------------------
    # create pyPrimer context, and define a geometry (using our hard-coded box)
    # -----------------------------------------------------------------------------
    vtxCount = vertices.shape[0]
    idxCount = indices.shape[0]
    print("#pyPrimer: creating triangle mesh: %d vertices, %d indices" % (vtxCount, idxCount))
    mesh  = pypr.create_mesh(0, vertices, indices)
    model = pypr.create_model([mesh])
    # -----------------------------------------------------------------------------
    # creating a grid of sampelRes x sampleRes parallel rays starting on a XY
    # plane and going in Z direction
    # -----------------------------------------------------------------------------
    sampleRes=4
    plane00=np.array([-2,-2,-3], dtype=float_t)
    planeDU=np.array([4,0,0], dtype=float_t)
    planeDV=np.array([0,4,0], dtype=float_t)
    rayDir=[0,0,1]
    # compute rays and hits as a python list
    rayType = [
        ('org_x', np.float32),
        ('org_y', np.float32),
        ('org_z', np.float32),
        ('t_min', np.float32),
        ('dir_x', np.float32),
        ('dir_y', np.float32),
        ('dir_z', np.float32),
        ('t_max', np.float32) ]
    rays=[]
    hitType = [
        ('primID', np.int32),
        ('geomID', np.int32),
        ('instID', np.int32),
        ('t', np.float32),
        ('u', np.float32),
        ('v', np.float32) ]
    hits=[]
    for i in range(0,sampleRes):
        for j in range(0,sampleRes):
            org_x = plane00[0] + (i+.5)/sampleRes * planeDU[0] + (j+.5)/sampleRes * planeDV[0]
            org_y = plane00[1] + (i+.5)/sampleRes * planeDU[1] + (j+.5)/sampleRes * planeDV[1]
            org_z = plane00[2] + (i+.5)/sampleRes * planeDU[2] + (j+.5)/sampleRes * planeDV[2]
            dir_x = rayDir[0]
            dir_y = rayDir[1]
            dir_z = rayDir[2]
            tmin = 0
            tmax = 1e20
            rays.append((org_x,org_y,org_z,tmin,dir_x,dir_y,dir_z,tmax))
            hits.append((-1,-1,-1,0,0,0))
    # ... and upload to numpy arrays
    rays=np.array(rays,dtype=rayType)
    hits=np.array(hits,dtype=hitType)
    print('rays to be traced:')
    print(rays)
    # now trace ...
    model.trace(rays,hits)
    # and done:
    print('hits received as a result of this trace (in full AoS form):')
    print(hits)
    print('and just the per-ray hit.primIDs of the prims hit by these rays:')
    print(hits['primID'])

main()
