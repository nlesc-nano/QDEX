Installation Guide
==================

Because `miniBSE` relies on C++ extensions compiled against `Libint2` and `Eigen3`, **Conda / Micromamba** is the recommended installation method.

Method 1: Conda / Micromamba (Recommended)
------------------------------------------

Our `environment.yml` handles the installation of C++ compilers, `CMake`, `Eigen3`, and `Libint2`, saving you the hassle of system-level configurations.

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/nlesc-nano/miniBSE.git
      cd miniBSE

2. Create and activate the conda environment:

   .. code-block:: bash

      micromamba env create -f environment.yml
      micromamba activate minibse_env

   *(Note: The environment file automatically installs `miniBSE` in editable mode via pip).*

Method 2: Manual Pip Installation
---------------------------------

If you already have `CMake` (>= 3.16), a C++17 compiler, `Eigen3`, and `Libint2` installed natively on your operating system:

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/nlesc-nano/miniBSE.git
      cd miniBSE

2. Install Python dependencies and build the C++ extension:

   .. code-block:: bash

      pip install -r requirements.txt
      pip install -e .

Verifying the Installation
--------------------------

To verify that `miniBSE` is correctly installed, run the CLI help command:

.. code-block:: bash

   minibse --help

You can also run the built-in test suite:

.. code-block:: bash

   pytest tests/
