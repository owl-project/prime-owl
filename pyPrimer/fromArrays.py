#!/usr/bin/python3

import pypr
import math
import sys
import numpy as np

def reshape_dim(arr, vec_dim):
    assert arr.shape[0] % vec_dim == 0
    num_elements = arr.shape[0] // vec_dim
    return arr.reshape((num_elements,vec_dim))

def read_arr_of_vecs(filename,vec_dim,dtype):
    data_arr=np.fromfile(filename,dtype=dtype)
    data_arr=reshape_dim(data_arr,vec_dim)
    print("data in " + filename + " has shape " + str(data_arr.shape) + " and looks like this:")
    print(data_arr)
    return data_arr

def main():

    flt_type = np.float32
    int_type = np.int32

    if len(sys.argv) != 2:
        print("Incorrect number of arguments!")
        print("Usage: python3 fromArrays.py <path_to_inputfiles>")
        print("Terminating!")
        sys.exit(1)

    fileNameBase=sys.argv[1]+"fromArrays-testData"

    print("#fromArrays: creating context...")
    pypr.context_create()

    print("#fromArrays: reading triangle mesh...")
    mesh_vertices = read_arr_of_vecs(fileNameBase+".vertices",3,flt_type)
    mesh_indices = read_arr_of_vecs(fileNameBase+".indices",3,int_type)
    
    print("#fromArrays: creating mesh...")
    mesh=pypr.create_mesh(0,mesh_vertices,mesh_indices)
    model=pypr.create_model([mesh])
    
    print("#fromArrays: reading rays...")
    ray_orgs = read_arr_of_vecs(fileNameBase+".ray_orgs",3,flt_type)
    ray_dirs = read_arr_of_vecs(fileNameBase+".ray_dirs",3,flt_type)
    ray_tmins=np.fromfile(fileNameBase+".ray_tmins",dtype=flt_type)
    ray_tmaxs=np.fromfile(fileNameBase+".ray_tmaxs",dtype=flt_type)
    print("ray_tmins: " + str(ray_tmins))
    print("ray_tmaxs: " + str(ray_tmaxs))
    num_rays=ray_orgs.shape[0]
    
    print("#fromArrays: reading reference hits...")
    ref_hit_ints=np.fromfile(fileNameBase+".hit_ids",dtype=int_type)
    print("ref_hit_ints: " + str(ref_hit_ints))
    ref_hit_coords = read_arr_of_vecs(fileNameBase+".hit_coords",3,flt_type)


    print("#fromArrays: allocing hits...")
    hit_ints=np.ndarray((num_rays,3),dtype=int_type)
    hit_floats=np.ndarray((num_rays,3),dtype=flt_type)
    
    print("#fromArrays: tracing...")
    model.trace(ray_orgs,ray_dirs,ray_tmins,ray_tmaxs,
                    hit_ints,hit_floats)
    
main()
