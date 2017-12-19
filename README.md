# Phonopy-Spectroscopy

Phonopy-Spectroscopy is a project to add the capability to simulate vibrational spectra to the Phonopy code.[[1](#Ref1)]

The software package consists of a Python module, `SpectroscoPy`, along with a set of command-line scripts for working with output from Phonopy and VASP.


## Features

* Calculate infrared (IR) intensities from Phonopy or VASP calculations.

* Calculate Raman-activity tensors and scalar-averaged intensities within the far-from-resonance approximation.

* Prepare peak tables including assigning modes to irreducible representations (Phonopy interface).

* Output customisable simulated spectra with support for multiple unit systems and simulated instrumental broadening.

* Include first-principles mode linewidths from Phono3py[[2](#Ref2)] calculations (Phonopy interface).


## Installation

The code depends on the `NumPy`[[3](#Ref3)] and `PyYAML`[[4](#Ref4)] packages, and the Phonopy interface additionally requires the `Phonopy` Python library[[1](#Ref1)] and the H5py package.[[5](#Ref5)]
All four packages are available from PyPI (via `pip`) and on the Anaconda platform (`conda`).
Please see the documentation of the codes for instructions on how to install them on your system.

This code does not currently ship with a `setup.py` script.
After cloning or downloading and unpacking the repository, add it to your `PYTHONPATH` so that the command-line scripts can locate `SpectroscoPy`, e.g.:

`export PYTHONPATH=${PYTHONPATH}:/Volumes/Data/Repositories/Phonopy-Spectroscopy`

The command-line scripts are in the [Scripts](./Scripts) directory.
For convenience, you may wish to add this folder to your `PATH` variable, e.g.:

`export PATH=${PATH}:/Volumes/Data/Repositories/Phonopy-Spectroscopy/Scripts`

For a description of the command-line arguments the scripts accept, call them with the `-h` option, e.g.:

`phonopy-ir -h`


## Examples

1. [Benzene derivatives](./Examples/Benzene-Derivatives): Simulated IR spectra of isolated molecules using the Phonopy interface, compared to reference gas-phase spectra from the NIST database

2. [&alpha;-SiO<sub>2</sub>](./Examples/a-SiO2): Detailed simulation of the IR and Raman spectra of &alpha;-SiO<sub>2</sub> (quartz) using the Phonopy interface, including first-principles linewidths calculated using Phono3py, compared to spectra from the RRUFF database


## Citation

We haven't written a standalone paper on this code yet, but it is based on the formulae collected in:

J. M. Skelton, L. A. Burton, A. J. Jackson, F. Oba, S. C. Parker and A. Walsh, "Lattice dynamics of the tin sulphides SnS<sub>2</sub>, SnS and Sn<sub>2</sub>S<sub>3</sub>: vibrational spectra and thermal transport", *Physical Chemistry Chemical Physics* **19**, 12452 (**2017**), DOI: [10.1039/C7CP01680H](https://doi.org/10.1039/C7CP01680H) (this paper is open access)

If you use Phonopy-Spectroscopy in your work, please consider citing this paper and/or including a link to this GitHub repository when you publish your results.


## References

1. <a name="Ref1"></a>[https://atztogo.github.io/phonopy/](https://atztogo.github.io/phonopy/)
2. <a name="Ref2"></a>[https://atztogo.github.io/phono3py/](https://atztogo.github.io/phono3py/)
3. <a name="Ref3"></a>[http://www.numpy.org/](http://www.numpy.org/)
4. <a name="Ref4"></a>[http://pyyaml.org/wiki/PyYAML](http://pyyaml.org/wiki/PyYAML)
5. <a name="Ref5"></a>[http://www.h5py.org/](http://www.h5py.org/)
