# solpolpy
[![codecov](https://codecov.io/gh/punch-mission/solpolpy/branch/main/graph/badge.svg?token=835TUH7CKI)](https://codecov.io/gh/punch-mission/solpolpy)
[![CI](https://github.com/punch-mission/solpolpy/actions/workflows/CI.yml/badge.svg)](https://github.com/punch-mission/solpolpy/actions/workflows/CI.yml)
[![PyPI version](https://badge.fury.io/py/solpolpy.svg)](https://badge.fury.io/py/solpolpy)
[![DOI](https://zenodo.org/badge/498801276.svg)](https://zenodo.org/doi/10.5281/zenodo.10076326)

🚧🚧🚧 **UNDER DEVELOPMENT** 🚧🚧🚧

`solpolpy` is a solar polarization resolver based on [Deforest et al. 2022](https://doi.org/10.3847/1538-4357/ac43b6).
It converts between various polarization formats, e.g. from the native three triple version from observations
(also known as the MZP convention) to polarization brightness (pB) and total polarization (B), Stokes I, Q and U, etc.
An example of transforming the polarization basis using the LASCO/C2 images is 
shown in the image below.
The images at polarizing angles of -60°, 0° and +60° is shown in the top panel as 
Bm, Bz and Bp respectively.
The bottom panel shows the output of the `solpolpy` to convert the initial basis 
to the Stokes I, Q and U.
![Example result image](eg_image.png)

## Quickstart
`pip install solpolpy`

Note: We will be updating [the documentation](https://punch-mission.github.io/solpolpy/quickstart.html) soon! 

## Getting Help
Please open a discussion or issue for help. 

## Contributing
We encourage all contributions. If you have a problem with the code or would like to see a new feature, please open an issue. Or you can submit a pull request. 

## Code of Conduct
[Access here](CODE_OF_CONDUCT.md)

## Origin of the Name
`solpolpy` is just a combination of `sol` for solar, `pol` for polarization, and `py` for Python. 

