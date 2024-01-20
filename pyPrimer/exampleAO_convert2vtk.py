import numpy as np
import pyvista as pv

vertices = np.load("quadmesh_vtxs.npy")
indices = np.load("quadmesh_idxs.npy")
nodecount_col = 4*np.ones((indices.shape[0],1),dtype=indices.dtype)
# Concatenation is necessary because PyVista assumes 
# columns start with the number of nodes per face
faces = np.concatenate([nodecount_col,indices],axis=1) 
faces = faces.flatten()
cell_vals = np.load("patch_values.npy")
# define a surface or as PyVista calls it a 'mesh' object
surf = pv.PolyData(vertices, faces)
surf.cell_data['occluded_fraction']   = cell_vals

surf.save("exampleAO.vtk")
