Overview
========

``solpolpy`` supports the conversion among the polarization systems listed below:

(Use the keyword inside the parenthesis in string format for out_polarize_state in the calling sequence).

MZP (MZP): Triplet of images taken at -60°, 0°, and +60° polarizing angles.

B\ :sub:`T`, B\ :sub:`R` (BtBr): Pair of images with polarization along the tangential and radial direction with respect to the Sun respectively.

Stokes (Stokes): Total brightness (I), polarized brightness along vertical and horizontal axes (Q) and polarized brightness along ±45° (U) .

B, pB (BpB): Total brightness and 'excess polarized' brightness images pair respectively.

B, pB, pB' (Bp3): Analogous to Stokes I, Q and U, but rotates with α around the Sun unlike the later with fixed frame of reference of the instrument.

B, θ, p (Bthp): System with total brightness, angle and degree of polarization.

4 polarizer system (fourpol): For observations taken at sequence of four polarizer angles, i.e. 0°, 45°, 90° and 135°.

Npol (npol): Set of images taken at than three polarizing angles other than MZP

The equations required for the transformations from one system to another can be found in `Deforest et al. 2022`_.

Calling Sequence
================

.. code-block:: python

    import solpolpy as sp
    output = sp.resolve(input_data, out_polarize_state)

input_data: Input NDCollection of the type listed above or string of path to FITS files

out_polarize_state: String of required output polarization state; keywords inside the parenthesis listed above.

output: Output NDCollection of polarization state in 'out_polarize_state'

.. _Deforest et al. 2022: https://iopscience.iop.org/article/10.3847/1538-4357/ac43b6