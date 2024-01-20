# pyPrimer README
## Setup
**IMPORTANT**: In order to run most of the Python scripts in this directory, the `PYTHONPATH` environment variable must be set to the directory to the directory that contains the `pypr.*.so` file. 

* This file might be called `pypr.cpython-39-x86_64-linux-gnu.so`, for example, on a 64-bit x86 Linux 64 bit system using Python 3.9.
* I am assuming that the path is `../build`, but you may need to use something different if you created a different build directory with a different name/path when you ran CMAKE.

⚠️ The name of the `.so` file will depend on your operating system, Python version, and possibly other factors, e.g. your compiler.

## Dependencies
* [Numba](https://numba.pydata.org/) is required for `exampleAO.py` and  `exampleAO.ipynb`
  * Install with `conda install numba` or `pip install numba`
* [PyVista](https://numba.pydata.org/) is required for `exampleAO_convert2vtk.py`, and the "Visualization" section of `exampleAO.ipynb`
  * Install with `conda install -c conda-forge pyvista` or `pip install pyvista`

## Examples
* `exampleAO.py`: Runs the example application. Writes output as `patch_values.npy`, `quadmesh_idxs.npy`, and `quadmesh_vtxs.npy`
* `exampleAO_convert2vtk.py`: Converts the output of `exampleA0.py` to `.vtk` format, for use with tools like [ParaView](https://www.paraview.org/) and [VTK](https://vtk.org/).
* `exampleAO.ipynb`: Example application and visualization in [Jupyter Notebook](https://jupyter.org/) format.
