# solpolpy
[![codecov](https://codecov.io/gh/punch-mission/solpolpy/branch/main/graph/badge.svg?token=835TUH7CKI)](https://codecov.io/gh/punch-mission/solpolpy)

**UNDER DEVELOPMENT**

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
As this package is not currently released, you must clone the repo and install with `pip install .`. Then follow [the documentation](https://punch-mission.github.io/solpolpy/quickstart.html).

## Getting Help
Please contact [Ritesh Patel](mailto:ritesh.patel@swri.org) or [Marcus Hughes](mailto:marcus.hughes@swri.org). 

## Contributing
We encourage all contributions. If you have a problem with the code or would like to see a new feature, please open an issue. Or you can submit a pull request. 

## Code of Conduct
[Access here](CODE_OF_CONDUCT.md)

## Citation
Coming soon with the publication of a paper. 

## Origin of the Name
`solpolpy` is just a combination of `sol` for solar, `pol` for polarization, and `py` for Python. 

