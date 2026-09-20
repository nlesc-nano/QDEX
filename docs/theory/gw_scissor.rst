Quasiparticle (GW) Scissor Model
==================================

Standard semi-local DFT functionals (such as PBE) severely underestimate the fundamental band gap ($E_g = \mathrm{IP} - \mathrm{EA}$) of semiconductor nanostructures due to self-interaction error and the absence of the derivative discontinuity.

In `miniBSE`, the single-particle gap is corrected using an analytical **Scaled GW Quasiparticle Model** that accurately accounts for both quantum confinement and dielectric mismatch.

The Scaled GW Model
-------------------

The total quasiparticle scissor shift $\Delta_{\mathrm{GW}}$ applied to the DFT eigenvalues is formulated as:

.. math::

   \Delta_{\mathrm{GW}}(R, \epsilon_{\mathrm{out}}) = \Delta_{\mathrm{bulk}} + \Delta_{\mathrm{conf}}(R) + \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}})

1. **Bulk GW Shift ($\Delta_{\mathrm{bulk}}$)**:
   The difference between the bulk GW quasiparticle gap and the bulk PBE gap:

   .. math::

      \Delta_{\mathrm{bulk}} = E_g^{\mathrm{GW, bulk}} - E_g^{\mathrm{PBE, bulk}}

2. **Finite-Size Confinement Shift ($\Delta_{\mathrm{conf}}$)**:
   Interpolates between a high-accuracy finite vacuum anchor (computed with PBE0+GW on a monomer or small cluster of radius $R_0$) and the bulk limit:

   .. math::

      \Delta_{\mathrm{conf}}(R) = \frac{A}{(R + \ell)^p}

   where $A$, $\ell$, and $p$ are material-specific scaling parameters calibrated against full GW benchmark calculations.

3. **External Dielectric Polarization Shift ($\Delta_{\mathrm{pol}}$)**:
   When the nanocrystal is embedded in a solvent or dielectric medium of permittivity $\epsilon_{\mathrm{out}}$, image charge polarization screens the quasiparticle self-energy:

   .. math::

      \Delta_{\mathrm{pol}}(R, \epsilon_{\mathrm{out}}) = -\frac{e^2}{R} \left( \frac{\epsilon_\infty - \epsilon_{\mathrm{out}}}{\epsilon_\infty + \epsilon_{\mathrm{out}}} \right)

HOMO / LUMO Scissor Split Fractions
-----------------------------------

The total scissor $\Delta_{\mathrm{GW}}$ is partitioned between the valence and conduction bands according to the orbital self-energy shifts:

.. math::

   \varepsilon_i^{\mathrm{QP}} = \varepsilon_i^{\mathrm{DFT}} - f_{\mathrm{homo}} \, \Delta_{\mathrm{GW}} \quad (i \in \mathrm{occ})

.. math::

   \varepsilon_a^{\mathrm{QP}} = \varepsilon_a^{\mathrm{DFT}} + f_{\mathrm{lumo}} \, \Delta_{\mathrm{GW}} \quad (a \in \mathrm{virt})

where $f_{\mathrm{homo}} + f_{\mathrm{lumo}} = 1$. In materials like $\text{CsPbBr}_3$, $f_{\mathrm{homo}} \approx 0.43$ and $f_{\mathrm{lumo}} \approx 0.57$.
