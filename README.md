# Phonopy-Spectroscopy

Phonopy-Spectroscopy adds the capability to simulate vibrational spectra to the Phonopy code.[[1](#Ref1)]

The software package consists of a Python package, `phonopy_spectroscopy`, and command-line scripts for performing typical calculations.


## Upgrade

The current status of the development code is as follows:

* The "backend" Python API has been completely rewritten with a new object model.
* The infrared (IR) capability has been rewritten around computing the IR dielectric function, with plans to implement broad support for single-crystal and powder spectra measured in both on- and off-axis geometries.
* The Raman capability has been expanded to include complete functionality for simulating single-crystal and powder Raman experiments, including with energy-dependent Raman tensors and non-analytical corrections for non-centrosymmetric materials.
* Unit tests have been implemented with good code coverage.
* The command-line programs and associated API are currently being rewritten.


## Installation

The code requires the following:
* Python >= 3.10
* The typical scientific computing Python stack: NumPy[[2](#Ref3)], SciPy[[3](#Ref3)], Pandas[[4](#Ref4)], PyYaml[[5](#Ref5)] and H5py[[6](#Ref6)].
* The Phonopy and spglib Python libraries[[1](#Ref1), [7](#Ref7)].
* Numba[[8](#Ref8)] is optional but recommended.

We recommend using Anaconda (or Miniconda)[[9](#Ref9)] and installing the required packages in a virtual environment:

```bash
% conda create -n phonopy-spectroscopy
% conda activate phonopy-spectroscopy
% conda install -c conda-forge phonopy pandas numba
```

The current code can be cloned from the `develop` branch of the GitHub:

```bash
% git clone -b develop https://github.com/skelton-group/Phonopy-Spectroscopy
```

This code does not currently ship with a `setup.py` script or as a PIP/Conda package.
After cloning or downloading and unpacking the repository, add the `lib` subfolder to your `PYTHONPATH`, e.g.:

```bash
% export PYTHONPATH=${PYTHONPATH}:$(pwd)/Phonopy-Spectroscopy/lib
```

It is a good idea to run the unit tests to check your installation is working:

```bash
% cd Phonopy-Spectroscopy/tests
% python -m unittest discover
```


## Citation

Paper(s) on the new implementation and code are planned in the very near future.
For now, the citation for the previous version of Phonopy-Spectroscopy is:

J. M. Skelton, L. A. Burton, A. J. Jackson, F. Oba, S. C. Parker and A. Walsh, "Lattice dynamics of the tin sulphides SnS<sub>2</sub>, SnS and Sn<sub>2</sub>S<sub>3</sub>: vibrational spectra and thermal transport", *Physical Chemistry Chemical Physics* **19**, 12452 (**2017**), DOI: [10.1039/C7CP01680H](https://doi.org/10.1039/C7CP01680H) (open access)

If you use Phonopy-Spectroscopy in your work, please consider citing this paper and/or including a link to this GitHub repository when you publish your results.


## References

1. <a name="Ref1"></a>[https://github.com/phonopy/phonopy](https://github.com/phonopy/phonopy)
2. <a name="Ref2"></a>[https://www.numpy.org](https://www.numpy.org)
3. <a name="Ref3"></a>[https://scipy.org](https://scipy.org)
4. <a name="Ref4"></a>[https://pandas.pydata.org](https://pandas.pydata.org)
5. <a name="Ref5"></a>[https://pyyaml.org/wiki/PyYAML](https://pyyaml.org/wiki/PyYAML)
6. <a name="Ref6"></a>[https://www.h5py.org](https://www.h5py.org)
7. <a name="Ref7"></a>[https://spglib.readthedocs.io](https://spglib.readthedocs.io)
8. <a name="Ref8"></a>[https://numba.pydata.org](https://numba.pydata.org)
9. <a name="Ref7"></a>[https://www.anaconda.com](https://www.anaconda.com)
