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

``resolve`` takes two parameters: 1. the input data and 2. the desired output polarization system.
Valid polarization systems are described in the documentation.

.. autofunction:: solpolpy.resolve

The equations required for the transformations from one system to another can be found in `Deforest et al. 2022`_.

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
