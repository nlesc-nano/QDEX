Dielectric Screening & Kernels
===============================

The accuracy of Bethe-Salpeter excited states hinges on the treatment of the screened direct electron-hole attraction $W(\mathbf{r}, \mathbf{r}')$ and the bare exchange $v(\mathbf{r}, \mathbf{r}')$.

Electronic Resta Model
----------------------

In semiconductor nanoclusters, dielectric screening varies from the bulk high-frequency dielectric constant $\epsilon_\infty$ at large distances to unscreened vacuum behavior ($\epsilon = 1$) at short distances. 

`miniBSE` implements the **Resta model of electronic screening**:

.. math::

   W(r) = \frac{1}{\epsilon_\infty r} + \frac{1 - \epsilon_\infty^{-1}}{r} \exp\left( -\frac{r}{\lambda_s} \right)

where:
* $\epsilon_\infty$ is the optical dielectric constant of the bulk material.
* $\lambda_s$ is the **Thomas-Fermi screening length** of the valence electron gas:

.. math::

   \lambda_s = \sqrt{ \frac{\pi}{4 k_F} } = \left( \frac{\pi}{4 (3\pi^2 n_v)^{1/3}} \right)^{1/2}

where $n_v$ is the valence electron density of the crystal.

Asymptotic Limits:
* **Short-range ($r \ll \lambda_s$)**: $W(r) \to \frac{1}{r}$ (unscreened Coulomb attraction, preventing unphysical over-screening at atomic distances).
* **Long-range ($r \gg \lambda_s$)**: $W(r) \to \frac{1}{\epsilon_\infty r}$ (macroscopic bulk dielectric screening).

MNOK Atom-Centered Discretization
---------------------------------

To evaluate four-center electron-hole integrals efficiently without computing millions of atomic orbital two-electron integrals, `miniBSE` contracts transition densities into atom-centered charges using Mulliken population analysis:

.. math::

   q_A^{ia} = \sum_{\mu \in A} \sum_\nu C_{\mu i} S_{\mu \nu} C_{\nu a}

The electron-hole interaction matrix elements are then computed via the **Mataga-Nishimoto-Ohno-Klopman (MNOK)** damped interaction:

.. math::

   K_{ia, jb}^d = \sum_{A} \sum_{B} q_A^{ij} \, W_{AB}^{\mathrm{Resta}} \, q_B^{ab}

.. math::

   W_{AB}^{\mathrm{Resta}} = \frac{1}{\epsilon_\infty \sqrt{R_{AB}^2 + a_{AB}^2}} + \frac{1 - \epsilon_\infty^{-1}}{\sqrt{R_{AB}^2 + a_{AB}^2}} \exp\left( -\frac{\sqrt{R_{AB}^2 + a_{AB}^2}}{\lambda_s} \right)

where $a_{AB} = 2 / (\eta_A + \eta_B)$ is the Ohno-Klopman damping parameter determined by the atomic chemical hardness $\eta$ of atoms $A$ and $B$.

Bare Exchange Kernel ($K^x$)
----------------------------

The repulsive electron-hole exchange is unscreened:

.. math::

   K_{ia, jb}^x = \sum_{A} \sum_{B} q_A^{ia} \, \frac{1}{\sqrt{R_{AB}^2 + a_{AB}^2}} \, q_B^{jb}

This term is essential for singlet-triplet splitting and oscillator strength intensity borrowing.
