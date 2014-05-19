#! /usr/bin/env python
"""
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

Equilibrium richness can be found by setting equal and solving for $R$.

\begin{align*}
    s_0 R^{a} = e_0 R^{b} \\
    \log(s_0) + a \log(R) = e_0 + b \log(R) \\
    (a-b) \log(R) = \log(e_0) - \log(s_0) \\
    \log(R) = \frac{ \log(e_0) - \log(s_0) } { a - b} \\
\end{align*}

"""

import math
import sys

def f(s0, e0, a=-0.5, b=0.5):
    logR = (math.log(e0) - math.log(s0)) / (a-b)
    R = math.exp(logR)
    print("s0 = {:20.5f}, e0 = {:20.5f}, logR = {:20.5f}, R = {:20.5f} [s0/e0 = {:20.5f} ]".format(s0, e0, logR, R, s0/e0))

for s0 in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 10000]:
    for e0 in [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0, 10, 100, 1000, 10000]:
        f(s0=s0, e0=e0)

