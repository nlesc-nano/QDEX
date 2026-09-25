Cubes
=====

Part of :doc:`/exciton_analysis/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.nto``
* Callable: ``qdex.nto.analyze_nto_state``
* CLI: ``--nto, --nto-states``
* YAML: ``analysis.nto, analysis.nto_states``

.. code-block:: python

   analyze_nto_state(solver, vec, energy_ev, f_osc, state_index, coords, symbols, mu_ia=None, soc_U=None, top_n=3, context=None)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part5_exciton_analysis/index.rst:196-208``

4. 3D Volumetric Visualization (.cube Files)
--------------------------------------------

``QDEX`` generates 3D volumetric Gaussian ``.cube`` files for direct visualization in VMD, PyMOL, or ChimeraX:

* **Hole Density**: ``exciton_S1_hole.cube``: :math:`\rho_h(\mathbf{r}) = \sum_{ia} |X_{ia}^1|^2 |\phi_i(\mathbf{r})|^2`
* **Electron Density**: ``exciton_S1_elec.cube``: :math:`\rho_e(\mathbf{r}) = \sum_{ia} |X_{ia}^1|^2 |\phi_a(\mathbf{r})|^2`
* **Difference Density**: ``exciton_S1_diff.cube``: :math:`\Delta \rho(\mathbf{r}) = \rho_e(\mathbf{r}) - \rho_h(\mathbf{r})`

Positive values in the difference cube indicate regions of net electron accumulation (photo-induced negative charge), while negative values indicate regions of hole accumulation (photo-induced positive charge).

---
