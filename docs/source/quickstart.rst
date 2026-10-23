Quickstart
===========

.. note::

    We recommend first following the :doc:`example` if you're unfamiliar with the package.
    It is more comprehensive than this page.

To get started, install the package with ``pip install solpolpy``.

Resolve between systems
-------------------------

This package allows you to convert between different polarization systems relevant to solar physics.
You can find a detailed description of these systems with their equations in `Deforest et al. 2022`_.

The easiest way to interact with solpolpy is through the ``resolve`` method:

.. code-block:: python

    from solpolpy import resolve, load_data
    paths = ['path_to_image0.fits', 'path_to_image1.fits', 'path_to_image2.fits']
    out_system = "tBpB"
    input_collection = load_data(paths)
    output_collection = resolve(input_collection, out_system)

    # to access the data, just access the appropriate cube
    print(output_collection['tB'].data)

``resolve`` takes two (or more depending on the polarization systems) parameters:
1. the input data and 2. the desired output polarization system.
Valid polarization systems are described in the documentation.

.. autofunction:: solpolpy.resolve

The equations required for the transformations from one system to another can be found in `Deforest et al. 2022`_.

Loading PUNCH data
------------------

Use ``load_punch`` with a PTM, PIM, or PAM FITS path, or an already-loaded
PUNCH NDCube. The returned collection can be passed directly to ``resolve``:

.. code-block:: python

    from solpolpy import load_punch, resolve

    collection = load_punch("/data/PUNCH_L3_PIM_observation.fits")
    bp3 = resolve(collection, "bp3")

    # Already-loaded Punchbowl cubes work too.
    collection = load_punch(punch_cube)

For individual M/Z/P frames, supply three paths or cubes in any order:

.. code-block:: python

    collection = load_punch(["/data/PP1.fits", "/data/PM1.fits", "/data/PZ1.fits"])
    bp3 = resolve(collection, "bp3")

The standard supported products are:

=====  ========  =======================  ===============
Level  Product   Input layout             Collection keys
=====  ========  =======================  ===============
1      PM/PZ/PP  Separate M/Z/P frames    M, Z, P
2      PTM       M/Z/P OBSLAYRs           M, Z, P
3      PIM       M/Z/P OBSLAYRs           M, Z, P
3      PTM       tB/pB/pB-prime OBSLAYRs  tB, pB, pBp
3      PAM       tB/pB/pB-prime OBSLAYRs  tB, pB, pBp
=====  ========  =======================  ===============

``pBp`` is solpolpy's key for pB-prime, stored as ``Polar_pBp`` in PUNCH
headers. Plane labels in ``OBSLAYR1-3`` determine whether a mosaic contains
M/Z/P or tB/pB/pBp, so Level 2 and Level 3 PTM products are distinguished
correctly. Experimental PAM files saved with M/Z/P layers are also supported.
Individual frames are identified by ``POLAR`` and must specify
``POLARREF``. Choose frames from the intended observation on a common spatial
grid; the loader does not match times or reproject them.

The loader preserves metadata, spatial WCS (including distortion tables),
and masks. It decodes the FITS signal-to-noise extension into standard
deviations on the input collection; errors on loaded cubes are used directly.
Missing errors remain absent. Nonfinite signals and invalid errors
are masked, while finite negative signals remain valid. A boolean spatial
``mask=`` can add a shared mask. Observation-level quality flags such as
``OUTLIER`` and ``BADPKTS`` remain available in metadata for your own selection.
No Punchbowl installation is required when loading FITS paths.

.. autofunction:: solpolpy.load_punch

IMAX Effect
-------------

Wide field polarizing imagers such as the Wide Field Imager (WFI) of PUNCH suffer from
foreshortening of the polarizer angle across their field of view (FOV). The observed polarization angles
deviate from the ideal based on their spatial location.
The foreshortening effect draws parallels with IMAX 3D
presentation which use linear polarizer systems and a wide screen.
In ``solpolpy``, the IMAX effect is corrected in input data using equation 44 from `Deforest et al. 2022`_.
Details about IMAX effect on WFI data will be soon published as a research article.

solpolpy supports this correction for the ``imax_effect`` within the `~solpolpy.resolve` function.

Plotting results
-----------------

You can plot any input or output collection using the `plot_collection` function provided.

.. code-block:: python

    from solpolpy import plot_collection
    plot_collection(input_collection)
    plot_collection(output_collection)

It has many options:

.. autofunction:: solpolpy.plot_collection


.. _Deforest et al. 2022: https://iopscience.iop.org/article/10.3847/1538-4357/ac43b6
