Spin-Orbit Coupling (SOC)
===========================

In heavy-element semiconductors (such as lead halide perovskites $\text{CsPb}X_3$, bismuth halides $\text{Cs}_3\text{Bi}_2\text{Br}_9$, or III-V/II-VI quantum dots $\text{InAs}$, $\text{HgTe}$), relativistic **Spin-Orbit Coupling (SOC)** fundamentally alters the electronic structure. 

In lead halide perovskites, SOC splits the Pb $6p$ conduction band into a lower $j = 1/2$ doublet and an upper $j = 3/2$ quartet, reducing the band gap by $\approx 0.65\text{ eV}$ and dramatically accelerating electron cooling.

Relativistic Pseudopotential Formalism
--------------------------------------

`miniBSE` evaluates SOC within the separable **Goedecker-Teter-Hutter (GTH)** pseudopotential framework:

.. math::

   V_{\mathrm{SO}} = \sum_I \sum_{l, m} \sum_{i, j} |p_i^{lm}\rangle \, k_{ij}^{l} \, \mathbf{L} \cdot \mathbf{S} \, \langle p_j^{lm}|

where:
* $|p_i^{lm}\rangle$ are Gaussian-type projector functions centered on atom $I$.
* $k_{ij}^l$ are material-specific SOC coupling coefficients parsed from `GTH_SOC_POTENTIALS.txt`.
* $\mathbf{L} = (L_x, L_y, L_z)$ are the orbital angular momentum operators in the complex spherical harmonic basis.
* $\mathbf{S} = \frac{1}{2} (\sigma_x, \sigma_y, \sigma_z)$ are the Pauli spin matrices.

Spinor Hamiltonian Assembly
---------------------------

The single-particle Hamiltonian in the 2-component spinor basis $(\alpha, \beta)$ is:

.. math::

   \mathbf{H}_{\mathrm{total}} = \begin{pmatrix}
     \mathbf{H}_0 - \frac{1}{2} \mathbf{H}_z & -\frac{1}{2} (\mathbf{H}_x - i \mathbf{H}_y) \\
     -\frac{1}{2} (\mathbf{H}_x + i \mathbf{H}_y) & \mathbf{H}_0 + \frac{1}{2} \mathbf{H}_z
   \end{pmatrix}

where $\mathbf{H}_0 = \operatorname{diag}(\varepsilon_m^{\mathrm{DFT}})$ contains the spin-free molecular orbital eigenvalues, and the Cartesian SOC blocks are:

.. math::

   \mathbf{H}_x = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_x \, \mathbf{B}_{\mathrm{mo}}^T, \quad
   \mathbf{H}_y = \frac{i}{2} \mathbf{B}_{\mathrm{mo}} \, \tilde{\mathbf{K}}_y \, \mathbf{B}_{\mathrm{mo}}^T, \quad
   \mathbf{H}_z = \frac{1}{2} \mathbf{B}_{\mathrm{mo}} \, \mathbf{K}_z \, \mathbf{B}_{\mathrm{mo}}^T

Here, $\mathbf{B}_{\mathrm{mo}} = \mathbf{C}_{\mathrm{act}}^T \mathbf{B}_{\mathrm{raw}}$ is the projection of the AO-projector overlap matrix into the active space.

High-Performance Sparse Assembly
--------------------------------

For large nanocrystals ($> 10,000$ AOs, $2,500$ MOs, $5,000$ spinors), naive loop-based assembly of $\mathbf{H}_{\mathrm{SO}}$ takes minutes. `miniBSE` utilizes three key optimizations:

1. **Sparse CSR Projector Representation**:
   All 1,750 angular projector blocks are pre-assembled into global sparse CSR matrices ($\mathbf{K}_x, \tilde{\mathbf{K}}_y, \mathbf{K}_z$) and cached across frames.
2. **Real-Valued DGEMM**:
   Exploiting the fact that $\mathbf{B}_{\mathrm{mo}}$ and the $K$ matrices are purely real (`float64`), complex matrix multiplications are replaced with real BLAS DGEMM.
3. **In-Place Block Construction**:
   Constructs $\mathbf{H}_{\mathrm{total}}$ directly into a contiguous $(2N_{\mathrm{mo}} \times 2N_{\mathrm{mo}})$ memory buffer, completely avoiding large temporary array allocations.

Diagonalizing $\mathbf{H}_{\mathrm{total}}$ yields the relativistic spinor energies $\varepsilon_k^{\mathrm{spinor}}$ and spinor expansion coefficients $\mathbf{U}$:

.. math::

   |\psi_k^{\mathrm{spinor}}\rangle = \sum_m \left[ U_{mk}^\alpha |\phi_m\rangle \otimes |\alpha\rangle + U_{mk}^\beta |\phi_m\rangle \otimes |\beta\rangle \right]
