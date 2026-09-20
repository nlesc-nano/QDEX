Bethe-Salpeter Equation & TDA
==============================

The **Bethe-Salpeter Equation (BSE)** is the rigorous many-body Green's function framework for describing optical excitations, neutral two-particle (electron-hole) states, and excitonic binding in molecules and nanostructures.

Tamm-Dancoff Approximation (TDA)
--------------------------------

In standard linear-response BSE, the full Hamiltonian contains both resonant (creation of an electron-hole pair) and anti-resonant (simultaneous de-excitation) blocks:

.. math::

   \begin{pmatrix}
     \mathbf{A} & \mathbf{B} \\
     -\mathbf{B}^* & -\mathbf{A}^*
   \end{pmatrix}
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}
   = \Omega
   \begin{pmatrix}
     \mathbf{X} \\
     \mathbf{Y}
   \end{pmatrix}

Under the **Tamm-Dancoff Approximation (TDA)**, coupling to de-excitations is neglected (:math:`\mathbf{B} \approx 0`). This reduces the problem to a standard Hermitian eigenvalue problem:

.. math::

   \mathbf{A} \mathbf{X}_S = \Omega_S \mathbf{X}_S

where :math:`\Omega_S` is the excitation energy of exciton state :math:`S`, and :math:`\mathbf{X}_S = (X_{ia}^S)` represents the expansion coefficients in the single-particle electron-hole basis :math:`|ia\rangle = a_a^\dagger a_i |0\rangle`:

.. math::

   |\Psi_S\rangle = \sum_{i \in \mathrm{occ}} \sum_{a \in \mathrm{virt}} X_{ia}^S \, |ia\rangle

Matrix Elements of the Resonant Matrix A
----------------------------------------

The resonant matrix elements $A_{ia, jb}$ consist of three distinct contributions:

.. math::

   A_{ia, jb} = \left( \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} \right) \delta_{ij} \delta_{ab} + 2 K_{ia, jb}^x - K_{ia, jb}^d

1. **Quasiparticle Diagonal ($D_{\mathrm{QP}}$)**:
   The single-particle energy difference between virtual orbital $a$ and occupied orbital $i$, corrected by the GW quasiparticle scissor:
   
   .. math::

      \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} = (\varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}) + \Delta_{\mathrm{GW}}

2. **Bare Exchange ($K_{ia, jb}^x$)**:
   The unscreened repulsive electron-hole exchange interaction responsible for the singlet-triplet splitting and local-field effects:

   .. math::

      K_{ia, jb}^x = \iint \frac{\phi_i(\mathbf{r}) \phi_a(\mathbf{r}) \, \phi_j(\mathbf{r}') \phi_b(\mathbf{r}')}{|\mathbf{r} - \mathbf{r}'|} \, d\mathbf{r} \, d\mathbf{r}'

3. **Screened Direct Attraction ($K_{ia, jb}^d$)**:
   The attractive Coulomb interaction between the electron and the hole, screened by the dielectric response of the material:

   .. math::

      K_{ia, jb}^d = \iint \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) \, W(\mathbf{r}, \mathbf{r}') \, \phi_a(\mathbf{r}') \phi_b(\mathbf{r}') \, d\mathbf{r} \, d\mathbf{r}'

Excitation Frameworks
---------------------

`miniBSE` supports four excitation levels:

* **`bse` (Full TDA-BSE)**:
  Full diagonalization (or Davidson iterative solving) of the complete matrix $A_{ia, jb}$. Accurately captures configuration interaction, spatial exciton delocalization, and oscillator strength redistribution.
* **`diagonal_bse`**:
  Omits off-diagonal configuration mixing ($ia \neq jb$), evaluating:

  .. math::

     \Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}} + 2 K_{ia, ia}^x - K_{ia, ia}^d

  Extremely fast and memory-efficient; ideal for non-adiabatic dynamics in dense manifolds.
* **`independent_qp`**:
  Non-interacting quasiparticle transitions: $\Omega_{ia} = \varepsilon_a^{\mathrm{QP}} - \varepsilon_i^{\mathrm{QP}}$.
* **`independent_dft`**:
  Bare DFT Kohn-Sham orbital energy differences: $\Omega_{ia} = \varepsilon_a^{\mathrm{DFT}} - \varepsilon_i^{\mathrm{DFT}}$.

Davidson Iterative Eigensolver
------------------------------

For large systems where full dense diagonalization ($O(N_{\mathrm{pairs}}^3)$) is prohibitive, `miniBSE` provides an optimized **Davidson iterative solver**. It computes the lowest $k$ excited states using matrix-vector products $A \mathbf{v}$ without ever storing the full $N_{\mathrm{pairs}} \times N_{\mathrm{pairs}}$ matrix in memory.
