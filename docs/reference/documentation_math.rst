Writing mathematical notation in the documentation
===================================================

Sphinx does not reliably interpret an inline ``:math:`` role nested inside
reStructuredText strong or emphasis markup. Keep formatting roles adjacent:

.. code-block:: rst

   **Screened exchange** (:math:`\Delta\Sigma^{\mathrm{SEX}}`)
   **Anchor 1: vacuum cluster** (:math:`R_0`)

Do not wrap either complete expression in ``**...**``. The same rule applies
to table cells and numbered-list labels. Use ``.. math::`` for displayed
equations and leave a blank line before and after the directive.

After changing the docs, build with Sphinx warnings treated as errors and
run ``python audit/check_rendered_math.py <html-build-directory>``. Sphinx can
otherwise build successfully while silently showing a literal ``:math:``
fragment to readers.
