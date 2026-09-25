QDEX documentation
==================


.. figure:: /_static/figures/pipeline.svg
   :width: 100%
   :alt: pipeline

   QDEX workflow. Each stage links to a topic group below; optional stages (SOC, dynamics, recombination) reuse the orbitals, QP energies and kernels of the earlier stages.

QDEX connects post-DFT electronic structure, quasiparticle models, BSE/TDA excitations and carrier dynamics for molecules and semiconductor nanoclusters.

The topic pages place the physical equations beside the QDEX function, flags and YAML configuration that implement them. Model-status notes distinguish established results from QDEX approximations.

Quick entry points: :doc:`getting_started/installation`, :doc:`getting_started/quickstart`, :doc:`getting_started/configuration`, and :doc:`api/index`.

.. toctree::
   :maxdepth: 2
   :caption: Start here

   getting_started/installation
   getting_started/quickstart

.. toctree::
   :maxdepth: 2
   :caption: Electronic structure and relativity

   electronic_structure/index
   relativity/index

.. toctree::
   :maxdepth: 2
   :caption: Interactions, quasiparticles and excitons

   interactions/index
   quasiparticles/index
   excitons/index
   validation/index

.. toctree::
   :maxdepth: 2
   :caption: Analysis, dynamics and spectroscopy

   exciton_analysis/index
   dynamics/index
   spectroscopy/index
   recombination/index

.. toctree::
   :maxdepth: 2
   :caption: Workflows and reference

   workflows/index
   reference/index
   getting_started/configuration
   api/index

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
