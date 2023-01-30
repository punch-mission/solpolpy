# Quick Start

## Overview

`solpolpy` supports the conversion among the following polarization systems:

Npol: Set of images taken at more than three polarizing angles 
(e.g. 0°, 45°, 90° and 135°). 

MZP: Triplet of images taken at -60°, 0°, and +60° polarizing angles.

B<sub>T</sub>, B<sub>R</sub> : Pair of images with polarization along the 
tangential and radial direction with respect to the Sun respectively.

Stokes: Total (unpolarized) brightness (I), polarized brightness along vertical and horizontal axes (Q) and 
polarized brightness along &pm;45° (U) .

B, pB: Total (unpolarized) brightness and 'excess polarized' brightness images pair respectively.

B, pB, pB': Analogous to Stokes I, Q and U, but rotates with &alpha; around the 
Sun unlike the later with fixed frame of reference of the instrument.

B, θ, p: System with total (unpolarized) brightness, angle and degree of polarization.

The equations required for the transformations from one system to another can be found in [Deforest et al. 2022](https://doi.org/10.3847/1538-4357/ac43b6).

## Calling Sequence
```
import solpolpy as sp
output = sp.resolve(input_data, out_polarize_state)
```

input_data: Input NDCollection of the type listed above 

out_polarize_state: String of required output polarization state 

output: Output NDCollection of polarization state in 'out_polarize_state'
