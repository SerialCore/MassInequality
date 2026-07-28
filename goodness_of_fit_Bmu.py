#!/usr/bin/env python3
"""
Goodness of fit for B(mu) in Eq. (10) of the manuscript.

Uses the six effective (non-vanishing) two-body bindings from Table I
and the three-parameter form
  B(mu) = a0 + a1 (mu/Lambda)^2 + a2 ln((mu/Lambda)^2),  Lambda = 200 MeV,
with the rounded coefficients quoted in the paper:
  a0 = -45.5,  a1 = -1.3,  a2 = -75.6  (MeV).

Reports: per-point residuals, SSR, rms, R^2, and unweighted chi^2/dof
(unit weights of 1 MeV^2).

Usage:
  python3 goodness_of_fit_Bmu.py
"""
from __future__ import annotations

import math
from pathlib import Path

# Table I: effective bindings only (vanishing nn, sn, cn, bn excluded from fit)
POINTS = [
    # (label, mu [MeV], B [MeV])
    ("ss", 164.0, -40.6),
    ("cs", 267.2, -76.0),
    ("bs", 305.6, -91.0),
    ("cc", 720.0, -258.9),
    ("bc", 1089.7, -352.0),
    ("bb", 2240.0, -566.8),
]

LAMBDA = 200.0  # MeV
A0, A1, A2 = -45.5, -1.3, -75.6  # paper Eq. (10)
N_PAR = 3


def B_fit(mu: float) -> float:
    x2 = (mu / LAMBDA) ** 2
    return A0 + A1 * x2 + A2 * math.log(x2)


def main() -> None:
    labels, mus, Bs = zip(*POINTS)
    fits = [B_fit(mu) for mu in mus]
    res = [b - f for b, f in zip(Bs, fits)]

    n = len(POINTS)
    dof = n - N_PAR
    ssr = sum(r * r for r in res)
    rms = math.sqrt(ssr / n)
    chi2_dof = ssr / dof  # unit weights (1 MeV^2)

    mean_B = sum(Bs) / n
    sst = sum((b - mean_B) ** 2 for b in Bs)
    r2 = 1.0 - ssr / sst

    lines = [
        "Goodness of fit for B(mu) [Eq. (10)]",
        "=" * 64,
        f"Form: B = {A0} + {A1} (mu/{LAMBDA})^2 + {A2} ln((mu/{LAMBDA})^2)",
        f"Data: {n} effective bindings, {N_PAR} free parameters, dof = {dof}",
        "",
        f"{'ij':>4} {'mu':>8} {'B_data':>10} {'B_fit':>10} {'residual':>10}",
    ]
    for lab, mu, b, f, r in zip(labels, mus, Bs, fits, res):
        lines.append(f"{lab:>4} {mu:8.1f} {b:10.1f} {f:10.2f} {r:10.2f}")
    lines += [
        "",
        f"SSR (sum of squared residuals) = {ssr:.2f} MeV^2",
        f"rms residual                   = {rms:.2f} MeV",
        f"R^2                            = {r2:.4f}",
        f"chi^2/dof (unit 1 MeV weights) = {chi2_dof:.2f}",
        "",
        "Note: the large unweighted chi^2/dof reflects O(10-25) MeV residuals",
        "in a simple three-parameter phenomenological form (especially in the",
        "light ss/cs/bs sector). R^2 ~ 0.99 shows that the overall concave",
        "trend vs reduced mass is well captured; the fit is used to motivate",
        "Jensen-like inequalities, not as a precision mass formula.",
    ]
    text = "\n".join(lines) + "\n"
    print(text)
    out = Path(__file__).resolve().parent / "goodness_of_fit_Bmu_results.txt"
    out.write_text(text)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
