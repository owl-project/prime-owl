#!/usr/bin/python3

import pypr
from numba import njit
import numpy as np
import sys

def reshape_dim(arr, vec_dim):
    assert arr.shape[0] % vec_dim == 0
    num_elements = arr.shape[0] // vec_dim
    return arr.reshape((num_elements,vec_dim))

def read_arr_of_vecs(filename,vec_dim,dtype):
    data_arr=np.fromfile(filename,dtype=dtype)
    data_arr=reshape_dim(data_arr,vec_dim)
    return data_arr

@njit
def frame_fromnormal(N, e0, e1, T):  # modifies T
  dx0 = np.cross(e0,N)
  dx1 = np.cross(e1,N)
  if np.dot(dx0,dx0) > np.dot(dx1,dx1):
    dx = dx0
  else:
    dx = dx1
  dx = dx/np.linalg.norm(dx) # normalize
  dy = np.cross(N,dx)
  dy = dy/np.linalg.norm(dy) # normalize
  T[0,:] = dx
  T[1,:] = dy
  T[2,:] = N
  # return np.array([dx,dy,N])

@njit
def xfm_vector(N, v, e0, e1, T): # assumes v is a row vector
  frame_fromnormal(N, e0, e1, T) # modifies T
  return v @ T

# ASSUMPTION: number of vertices per face is always the same
# ASSUMPTION: we don't care about "repeated" vertices... 
# i.e. we aren't really going by the original OBJ's indexing scheme.
def read_obj(filename, verts_per_face, int_type=np.int32, flt_type=np.float32):
    vtx_list = []
    idx_list = []
    curr_idx = 0
    with open(filename) as f:
        for line_num, line in enumerate(f.readlines()):
            trimmed_ln = line.strip()
            if(trimmed_ln): # check for string that is all whitespace 
                if trimmed_ln[0] == '#': # comment line
                    continue 
                tokens = trimmed_ln.split()
                if tokens[0] == "f":
                    assert len(tokens) - 1 == verts_per_face
                    idx_list.extend([curr_idx, curr_idx+1, curr_idx+2, curr_idx+3])
                    curr_idx += 4
                elif tokens[0] == "v":
                    for t in tokens[1:]:
                        vtx_list.append(t)
                else:
                    print("Unsupported OBJ feature in line: " + str(line_num))
                    print(line)
    
    vtx_arr = np.array(vtx_list,dtype=flt_type).reshape((-1,3))
    idx_arr = np.array(idx_list,dtype=int_type).reshape((-1,4))
    return vtx_arr, idx_arr

# For now, assuming we read the vertex array four at a time to obtain a quad.
# I.e. we aren't keeping track of face indices.
@njit
def create_rays(vtx_arr, samples_per_patch):
    num_verts = vtx_arr.shape[0] 
    assert num_verts % 4 == 0

    flt_type = vtx_arr.dtype 

    # pre-allocate arrays for later use
    e0 = np.array([1,0,0],dtype=flt_type)
    e1 = np.array([0,1,0],dtype=flt_type)
    T = np.zeros((3,3),dtype=flt_type)
    dir_sample = np.zeros(3,dtype=flt_type)

    origins = []
    directions = []
    tmins = []
    tmaxs = []
    for i in range(0,num_verts,4):
        for j in range(samples_per_patch):
            v0 = vtx_arr[i,:]
            v1 = vtx_arr[i+1,:]
            v3 = vtx_arr[i+3,:]
            du = v1 - v0
            dv = v3 - v0

            ksi1 = np.random.random() 
            ksi2 = np.random.random() 
            N = np.cross(du,dv)
            N = N/np.linalg.norm(N) # unit normal
            origin = v0 + ksi1*du + ksi2*dv 

            r1 = np.random.random() 
            r2 = np.random.random() 
            dirSampleZ = r1;
            dirSampleX = np.cos(2*np.pi*r2)*np.sqrt(1-r1*r1);
            dirSampleY = np.sin(2*np.pi*r2)*np.sqrt(1-r1*r1);
            # dirSample = np.array([dirSampleX, dirSampleY, dirSampleZ],flt_type)
            dir_sample[0] = dirSampleX
            dir_sample[1] = dirSampleY
            dir_sample[2] = dirSampleZ
            dir_transformed = xfm_vector(N, dir_sample, e0, e1, T) # modifies T

            tmin = 1.0e-3;
            tmax = 1.0e20;
            tmins.append(tmin)
            tmaxs.append(tmax)

            origins.append(origin[0])
            origins.append(origin[1])
            origins.append(origin[2])
            directions.append(dir_transformed[0])
            directions.append(dir_transformed[1])
            directions.append(dir_transformed[2])

    org_arr = np.array(origins,dtype=flt_type)
    dir_arr = np.array(directions,dtype=flt_type)
    tmin_arr = np.array(tmins,dtype=flt_type)
    tmax_arr = np.array(tmaxs,dtype=flt_type)

    return org_arr.reshape((-1,3)), dir_arr.reshape((-1,3)), tmin_arr, tmax_arr

def get_patch_num_type(quadmesh_vtxs):
    num_verts = quadmesh_vtxs.shape[0] 
    assert num_verts % 4 == 0
    num_patches = num_verts // 4  # integer division
    return num_patches, quadmesh_vtxs.dtype

# Assumes that each quad has exactly four unique vertices
@njit
def integrate_patches(hits, num_patches, flt_t):
    patch_values = np.zeros(num_patches,dtype=flt_t)

    num_hits = hits.shape[0]
    assert num_hits % num_patches == 0
    num_rays_per_patch = num_hits // num_patches # integer division
    for patchID in range(num_patches): # iterate over patches
        unoccluded = 0
        for i in range(num_rays_per_patch):
            if hits[patchID*num_rays_per_patch+i,2] < 0:
                unoccluded += 1
        patch_values[patchID] = unoccluded / num_rays_per_patch # float divison

    return patch_values

if __name__ == "__main__":

    flt_t = np.float32
    int_t = np.int32

    # load quadmesh
    geom_filename = "example_data/boxes_quadmesh.obj"
    quadmesh_vtxs, quadmesh_idxs = read_obj(geom_filename,4)

    samps_per_patch = 64 # TODO: read this from a config file

    np.random.seed()

    # JIT warmup (for Numba)
    num_warmup_verts = 4
    num_warmup_rands = 4*num_warmup_verts*samps_per_patch 
    create_rays(quadmesh_vtxs[0:num_warmup_verts], samps_per_patch) 

    # Run our actual computation after the warmup
    orgs, dirs, tmins, tmaxs = create_rays(quadmesh_vtxs, samps_per_patch)

    # load triangle mesh
    trimesh_verts = read_arr_of_vecs("example_data/boxes_trimesh.vertices",3,flt_t)
    trimesh_idxs = read_arr_of_vecs("example_data/boxes_trimesh.indices",3,int_t)
    print("Shape trimesh inds: " + str(trimesh_idxs.shape))

    pypr.context_create()
    mesh=pypr.create_mesh(0,trimesh_verts,trimesh_idxs)
    model=pypr.create_model([mesh])

    assert orgs.shape[0] == dirs.shape[0]
    num_rays = orgs.shape[0]
    hit_ints=np.ndarray((num_rays,3),dtype=int_t)
    hit_floats=np.ndarray((num_rays,3),dtype=flt_t)
    
    model.trace(orgs,dirs,tmins,tmaxs,hit_ints,hit_floats)

    num_patches, vtx_dtype = get_patch_num_type(quadmesh_vtxs)
    
    # JIT warmup
    num_hits = hit_ints.shape[0]
    assert num_hits % num_patches == 0
    num_rays_per_patch = num_hits // num_patches # integer division
    num_warmup_patches = 1
    num_warmup_hits = num_warmup_patches*num_rays_per_patch
    patch_values = integrate_patches(hit_ints[0:num_warmup_hits], num_warmup_patches, vtx_dtype)

    patch_values = integrate_patches(hit_ints, num_patches, vtx_dtype)

    np.save("quadmesh_vtxs.npy", quadmesh_vtxs)
    np.save("quadmesh_idxs.npy", quadmesh_idxs)
    np.save("patch_values.npy", patch_values)
