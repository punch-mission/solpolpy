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
    out_system = "BpB"
    input_collection = load_data(paths)
    output_collection = resolve(input_collection, out_system)

    # to access the data, just access the appropriate cube
    print(output_collection['B'].data)

``resolve`` takes two (or more depending on the polarization systems) parameters:
1. the input data and 2. the desired output polarization system.
Valid polarization systems are described in the documentation.

.. autofunction:: solpolpy.resolve

The equations required for the transformations from one system to another can be found in `Deforest et al. 2022`_.

IMAX Effect
-------------

Wide field polarizing imagers such as the Wide Field Imager (WFI) of PUNCH suffer from
foreshortening of the polarizer angle across their field of view (FOV). The observed polarization angles
deviate from the ideal based on their spatial location.
The foreshortening effect draws parallels with IMAX 3D
presentation which use linear polarizer systems and a wide screen.
In ``solpolpy``, the IMAX effect is corrected in input data using equation 44 from `Deforest et al. 2022`_.
Details about IMAX effect on WFI data will be soon published as a research article.

solpolpy supports this correction with the ``imax_effect`` keyword on the `~solpolpy.resolve` function.
If set to true, it corrects for the IMAX effect and converts the apparent non-ideal angle to ideal MZP configuration.

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
