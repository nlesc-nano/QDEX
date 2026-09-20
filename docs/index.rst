.. miniBSE documentation master file

Welcome to miniBSE's Documentation!
===================================

**miniBSE** is a high-performance, lightweight post-DFT exciton solver and non-adiabatic molecular dynamics (NAMD) engine designed for calculating, analyzing, and propagating the excited states of molecules and semiconductor nanoclusters.

By combining the flexibility of a Python interface with a fast C++ backend powered by `Libint2` and `Eigen3`, `miniBSE` offers researchers a scalable tool to investigate:

* **Excitonic Structure**: Solves the Bethe-Salpeter Equation (BSE) under the Tamm-Dancoff Approximation (TDA) with screened electron-hole interactions.
* **Excitation Frameworks**: Choose from coupled resonant BSE, diagonal BSE, quasiparticle independent transitions, or raw DFT transitions.
* **Spin-Orbit Coupling (SOC)**: Fully relativistic 2-component spinor basis using separable HGH pseudopotential projectors with sparse block-diagonal assembly.
* **Carrier Cooling Dynamics (NAMD)**: High-throughput Non-Adiabatic Molecular Dynamics via the deterministic Pauli Master Equation (PME) or Classical Path Approximation Fewest Switches Surface Hopping (CPA-FSSH).
* **Photoluminescence & Recombination**: *Ab initio* Einstein spontaneous emission ($k_{\mathrm{rad}}$), multi-phonon Energy Gap Law ($k_{\mathrm{nr}}$), and Photoluminescence Quantum Yield (PLQY).
* **Wavefunction Analysis**: Spatial exciton descriptors ($d_{eh}$, Pearson $R$, charge-transfer character), Projected DOS (PDOS), Crystal Orbital Overlap Population (COOP), and Fuzzy Band unfolding.

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/quickstart
   getting_started/configuration

.. toctree::
   :maxdepth: 2
   :caption: Physics & Theory

   theory/bse_tda
   theory/kernels
   theory/gw_scissor
   theory/soc

.. toctree::
   :maxdepth: 2
   :caption: Non-Adiabatic Dynamics (NAMD)

   namd/index

.. toctree::
   :maxdepth: 2
   :caption: Electronic Structure Analysis

   analysis/index

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
