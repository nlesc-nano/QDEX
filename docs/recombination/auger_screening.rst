Auger screening
===============

Part of :doc:`/recombination/index`.

.. rubric:: Theory and QDEX implementation

The detailed theory and worked equations follow below. The corresponding entry point is:

* Module: ``qdex.auger``
* Callable: ``qdex.auger.calculate_auger_rates``
* CLI: ``--auger, --auger-channel``
* YAML: ``auger.run, auger.channel``

.. code-block:: python

   calculate_auger_rates(C: np.ndarray, eps: np.ndarray, S: np.ndarray, atom_ao_ranges: List[Tuple[int, int]], coords: np.ndarray, atom_symbols: List[str], homo_idx: int, W_resta: Optional[np.ndarray]=None, material_name: Optional[str]='DEFAULT', eps_out: float=2.0, sigma_ev: float=0.05, broadening_mode: str='gaussian', channel: str='all', n_initial_elec: int=1, n_initial_hole: int=1, e_search_sigma_factor: float=4.0, spinor: bool=False, U_spinor_alpha: Optional[np.ndarray]=None, U_spinor_beta: Optional[np.ndarray]=None, lambda_reorg_ev: Optional[float]=None, temperature_k: float=300.0, eps_eff: Optional[float]=None, verbose: bool=True)

.. rubric:: Detailed derivations and reference material

.. rubric:: From ``docs/part8_auger/index.rst:263-273``

3. The Atom-Centered Resta Screening Contraction
------------------------------------------------

In standard quantum chemistry, evaluating :math:`V_{e' e_1, h e_2}` requires transforming four-center two-electron atomic orbital integrals :math:`(\mu \nu | \hat{W} | \lambda \sigma)`:

.. math::

   V_{e' e_1, h e_2} = \sum_{\mu \nu \lambda \sigma} C_{\mu e'}^* C_{\nu e_1} C_{\lambda h}^* C_{\sigma e_2} (\mu \nu | \hat{W} | \lambda \sigma)

For a nanocrystal with :math:`N_{\mathrm{ao}} = 10,000` basis functions, storing and transforming :math:`(\mu \nu | \lambda \sigma)` requires :math:`O(N_{\mathrm{ao}}^5)` operations and petabytes of memory, rendering *ab initio* Auger calculations completely intractable.


.. rubric:: From ``docs/part8_auger/index.rst:274-296``

Monopole Transition Charge Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To solve this scaling bottleneck, ``QDEX`` contracts the transition densities into atom-centered point charges using Mulliken-type population analysis:

.. math::

   q_A^{e_2 h} = \sum_{\lambda \in A} \sum_{\sigma=1}^{N_{\mathrm{ao}}} C_{\lambda e_2} S_{\lambda \sigma} C_{\sigma h}

.. math::

   q_B^{e' e_1} = \sum_{\mu \in B} \sum_{\nu=1}^{N_{\mathrm{ao}}} C_{\mu e'} S_{\mu \nu} C_{\nu e_1}

The continuous six-dimensional integral collapses into a fast, pairwise contraction over atomic sites :math:`A` and :math:`B`:

.. math::

   V_{\mathrm{dir}} = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{e_2 h} \; W_{AB}^{\mathrm{Resta}} \; q_B^{e' e_1} = (\mathbf{q}^{e_2 h})^T \, \mathbf{W}^{\mathrm{Resta}} \, \mathbf{q}^{e' e_1}

.. math::

   V_{\mathrm{exch}} = \sum_{A=1}^{N_{\mathrm{atoms}}} \sum_{B=1}^{N_{\mathrm{atoms}}} q_A^{e_1 h} \; W_{AB}^{\mathrm{Resta}} \; q_B^{e' e_2} = (\mathbf{q}^{e_1 h})^T \, \mathbf{W}^{\mathrm{Resta}} \, \mathbf{q}^{e' e_2}


.. rubric:: From ``docs/part8_auger/index.rst:297-314``

The Microscopic Resta Screening Kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The interaction matrix :math:`W_{AB}^{\mathrm{Resta}}` is evaluated using the **Resta model of electronic screening**:

.. math::

   W_{AB}^{\mathrm{Resta}} = \frac{1}{\epsilon_\infty \sqrt{R_{AB}^2 + a_{AB}^2}} + \frac{1 - \epsilon_\infty^{-1}}{\sqrt{R_{AB}^2 + a_{AB}^2}} \exp\left( -\frac{\sqrt{R_{AB}^2 + a_{AB}^2}}{\lambda_s} \right)

where:
* :math:`\epsilon_\infty` is the bulk high-frequency dielectric constant.
* :math:`a_{AB} = \frac{1}{2}(\eta_A^{-1} + \eta_B^{-1})` is the Ohno-Klopman damping parameter derived from atomic chemical hardness values :math:`\eta_A`.
* :math:`\lambda_s = 1 / k_s` is the valence Thomas-Fermi screening length:

.. math::

   k_s = \frac{\sqrt{\epsilon_\infty - 1}}{d_{\mathrm{NN}}}


.. rubric:: From ``docs/part8_auger/index.rst:315-335``

Local vs. Macroscopic Dielectric Screening in Auger Scattering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In optical absorption (BSE), the electron and hole interact across long distances (:math:`r \sim 15 - 35\text{ \AA} \gg \lambda_s \approx 1.5\text{ \AA}`), so the interaction is fully screened by the macroscopic dielectric constant :math:`\epsilon_\infty` (e.g. 4.8 for :math:`\mathrm{CsPbBr}_3`).

In Auger recombination, however, the process transfers a large energy :math:`\hbar \omega = E_g \approx 3.3\text{ eV}` and large momentum:

.. math::

   q = \frac{\sqrt{2 m_e^* E_g}}{\hbar} \approx 0.55\text{ \AA}^{-1}

This corresponds to a short collision impact parameter :math:`r_{\mathrm{Auger}} \sim 1/q \approx 1.8\text{ \AA}` (interatomic distance). Evaluating Resta's screening function at :math:`r = 1.8\text{ \AA}`:

.. math::

   \frac{1}{\epsilon(1.8\text{ \AA})} = \frac{1}{4.8} + \left(1 - \frac{1}{4.8}\right) e^{-0.65 \times 1.8} \approx 0.454 \implies \epsilon(r_{\mathrm{Auger}}) \approx 2.20

That local dielectric, about :math:`2.2` rather than :math:`\epsilon_\infty = 4.8`, is already the Resta kernel. It is not applied a second time. A constant :math:`\epsilon_{\mathrm{eff}} \approx 0.4\,\epsilon_{\mathrm{bulk}}` is the alternative used when the whole interaction is scaled by one number (Efros, as cited by Hou et al., Nat. Commun. 2019). Passing ``eps_eff`` does not multiply the Resta matrix.

---
