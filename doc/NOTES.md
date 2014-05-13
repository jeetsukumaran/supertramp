
The Simulation Life-Cycle
=========================

#. DISPERSAL
    -   For each island:
        - If dispersal event occurs, select lineage and destination island; record incoming migrant.
    -   For each island:
        - Incorporate incoming migrants into population.
#. DIVERSIFICATION
    -   (Lineage Birth)
        For each tip lineage, $S_i$:
        - If speciation event occurs:
            - Split lineage $S_i$ into two daughter lineages, $S_j$ and $S_k$;
              $S_i$ is now no longer a tip lineage, being replaced by $S_j$ and $S_k$.
            - Set the habitat type or niche of $S_j$ to the habitat type or niche of $S_i$.
            - If niche evolution occurs, set the habitat type of $S_k$ to a new
              randomly-sampled habitat type; otherwise, set the habitat type of
              $S_k$ to $S_i$.
            - In all islands in which $S_i$ occurs, add $S_j$ as a lineage (in
              the original habitat of $S_i$).
            - In all islands in which $S_i$ occurs, select one at random, and
              add $S_k$ (in the habitat type of $S_k$).
            - Remove $S_i$ from all islands/habitats.
    -   (Lineage Death)
        For each habitat in each island:
        -   If the number of lineages > carrying capacity:
            -   Remove lineages at random until the number of lineages == carrying capacity
            Else:
            -   Remove each lineage with probability = $\mu * (1-(K-n)/K)$, where:
                - \mu : global death rate
                - K : carrying capacity of the habitat
                - n : number of lineages in the habitat

Notes
=====
-   The current speciation process results in both daughter lineages occurring
    in the same island, and possibly the same habitat, at the moment of
    speciation.  It is easy enough to implement an alternative, e.g. complete
    allopatry, where each island receives either one daughter or another. My
    intuition is that there should be not much difference with a reasonable
    dispersal rate. But I could be wrong.
-   The carrying capacity of each habitat is reached relatively quickly, so for
    most of the simulation, the "remove lineages at random" regime is in effect
    rather than the weighted death probability.


Speciation Model
================

    Speciation rate = $s_0 R^a$
    Extinction rate = $e_0 R^b$

    where:

        * $R$ = local richness in habitat
        * $a < 0$
        * $b > 0$

$s_0$ and $e_0$ are positive constants that can be tuned to change overall relative
rates of speciation and extinction
I suggest to start that $a=-0.5$ and $b=0.5$  which would make overall
diversification rate go down linearly with increasing richness.

To change rate to probability over a timestep, use logistic function.

Equilibrium Richness
--------------------

Equilibrium richness can be found by setting equal and solving for $R$:

```latex

    \frac{s0}{e0}
```

Effect of $s0$
--------------

Prediction:


1e-1   1e+6
1e-2   1e+5
1e-3   1e+4
1e-4   1e+3
1e-5   1e+2
1e-6   1e+1

(run for a 1e+6 gens)

So, order of magnitude of number of speciations: s0 * ngens * x
x = 10, or s0/e0

