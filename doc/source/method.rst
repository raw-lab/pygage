How GAGE works
==============

GAGE's defining feature is a **two-level** test: a per-sample comparison of each
gene set against the whole-array background, combined across samples by a
meta-test. PyGAGE ports this 1:1 from the gage R sources.

1. Preparation
--------------

Per-gene fold changes are formed from reference vs sample columns
(``paired`` / ``unpaired`` / ``as.group`` / ``1ongroup``). See
:func:`pygage.core.GAGEPreparation.prepare_expression`.

2. Per-sample gene-set statistic
--------------------------------

For column *j* and set *S* of *n* present genes, with background mean ``mu_j`` and
variance ``s_j`` over all genes:

.. code-block:: text

   a  = var(S_j) / n          # set variance / set size
   b  = s_j       / n          # background variance / SET size  (GAGE's definition)
   df = (a + b)^2 / (a^2/(n-1) + b^2/(n-1))
   stat_j   = (mean(S_j) - mu_j) * (a + b)^(-1/2)
   p_up_j   = P(T_df > stat_j) ;  p_down_j = P(T_df < stat_j)

.. important::

   GAGE divides **both** variance terms by the *set* size ``n`` (``b = s/n``) and
   uses Welch degrees of freedom on ``n-1``. This is GAGE's statistic by
   definition — not a textbook two-sample Welch test — and PyGAGE matches it
   exactly. The z-test (PAGE-style) and a rank-based KS test are also available.

3. Cross-sample meta-summarisation
----------------------------------

The per-sample p-values are combined into one p-value per set. **Stouffer's Z** is
the default:

.. code-block:: text

   p.val     = Phi( sum_j qnorm(p_j) / sqrt(nc) )       # meta_method="stouffer"
   p.geomean = exp( -sum_j -log(p_j) / nc )
   stat.mean = mean_j stat_j
   q.val     = BH(p.val)                                 # per direction

A Fisher/gamma alternative is available with ``meta_method="fisher"``. Separate
**greater** and **less** tables report up- and down-regulated sets.
