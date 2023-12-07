Development
============
We encourage all contributions. Please see our `contribution guide first <https://github.com/punch-mission/punch-mission/blob/main/contributing.md>`_.


We recommend working in a virtual environment.
This can be created by running ``python -m venv venv``. Then, activate the environment with ``source venv/bin/activate``.
You can then install the required packages with ``pip install -r requirements_dev.txt``.

If at any time you run into issues, please contact us by :doc:`following the guidelines here <help>`.

Building the docs
------------------
The docs are built using ``sphinx``. First, you must install it and the other documentation requirements with ::

    pip install -r ./docs/requirements.txt
    pip install -r requirements.txt

Then, navigate to the ``docs`` directory and run ``make html`` to build the docs.

Note, that this does not rerun the example Jupyter notebook. You must manually execute that if you want those figures
to update.

Running tests
-------------
To run the tests for this package, run ``pytest`` in the repository base directory.
