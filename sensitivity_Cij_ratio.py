#!/usr/bin/env python3
"""
Sensitivity of Table-II sigma and Eq.(12) baryon masses to
r = C_ij / C_i\bar{j}  (notebook: gammaA).

Mirrors formulas in MassInequalities.nb (CMI Fittings, CMI Inequalities,
Binding Inequalities, sigma block, mass-prediction block).

Usage:
  python3 sensitivity_Cij_ratio.py
"""
from __future__ import annotations

import math
from pathlib import Path

# --- experimental masses (MeV), same as notebook ---
Nu, Delta = 939.6, 1232.0
Xi, Xis = 1314.86, 1531.8
Omega = 1672.45
Xicp, Xics = 2578.7, 2646.16
Omegac, Omegacs = 2695.2, 2765.9
Xibp, Xibs = 5935.0, 5955.33
Sigmabs = 5832.53
Omegab = 6045.2

pim, rhom = 134.98, 775.26
Km, Kms = 497.6, 895.55
Dm, Dms = 1864.84, 2006.85
Bm, Bms = 5279.66, 5324.71
Dsm, Dsms = 1968.35, 2112.2
Bsm, Bsms = 5366.92, 5415.4
etac, etab = 2983.9, 9398.7
Jpsi, Upsilon = 3096.9, 9460.3
Bcm, Bcms = 6274.47, 6332.0

Kmu, Kmd = 493.68, 497.61
Kmsu, Kmsd = 891.67, 895.55
Dmu, Dmd = 1864.84, 1869.66
Dmsu, Dmsd = 2006.85, 2010.26
Bmu, Bmd = 5279.34, 5279.66
Bmsn = 5324.71
piud = 139.57

# eta_s from Witten-like combination (notebook pseudoPhi[1])
phi0 = math.sqrt(Kmu**2 + Kmd**2 - piud**2)

# isospin averages used in Delta M (baryon light sector)
Sigma_suu, Sigma_sdd, Sigma_sud = 1382.8, 1387.2, 1383.7
Xisu, Xisd = 1531.8, 1535.0
Xicsu, Xicsd = 2645.1, 2646.16
Sigmacsuu, Sigmacsdd = 2518.41, 2518.48
Sigmacsud = 2517.4
Xicpu, Xicpd = 2578.2, 2578.7
Omegac_m = 2695.2  # Omega_c
# Xi_c' vs (Omega_c + Sigma_c)/2 uses notebook DeltaMcsn1 — use notebook EXP values from paper when needed

# Bindings from spin-averaged mesons (notebook extraction; baryon B = Bbar/2)
# Use paper Table I values for consistency with published numbers
Bbar = {
    "nn": 0.0,
    "sn": 0.0,
    "ss": -40.6,
    "cn": 0.0,
    "bn": 0.0,
    "cs": -76.0,
    "bs": -91.0,
    "cc": -258.9,
    "bc": -352.0,
    "bb": -566.8,
}
B = {k: v / 2.0 for k, v in Bbar.items()}

HYPF = 16.0 / 3.0 + 16.0  # meson hyperfine denominator


def meson_Cbar():
    return {
        "nn": (rhom - pim) / HYPF,
        "sn": (Kms - Km) / HYPF,
        "ss": (1019.46 - phi0) / HYPF,  # phi mass
        "cn": (Dms - Dm) / HYPF,
        "cs": (Dsms - Dsm) / HYPF,
        "cc": (Jpsi - etac) / HYPF,
        "bn": (Bmsn - Bmu) / HYPF,
        "bs": (Bsms - Bsm) / HYPF,
        "bc": (Bcms - Bcm) / HYPF,
        "bb": (Upsilon - etab) / HYPF,
    }


phi = 1019.46


def delta_M_table():
    """Experimental Delta M for the 14 systems entering sigma (notebook order)."""
    dMsn1 = 0.5 * (
        (Kmsu - 0.5 * (phi + rhom)) + (Kmsd - 0.5 * (phi + rhom))
    )
    dMcn1 = 0.5 * (
        (Dmsu - 0.5 * (Jpsi + rhom)) + (Dmsd - 0.5 * (Jpsi + rhom))
    )
    dMcn0 = 0.5 * (
        (Dmu - 0.5 * (etac + pim)) + (Dmd - 0.5 * (etac + pim))
    )
    dMbn1 = Bmsn - 0.5 * (Upsilon + rhom)
    dMbn0 = 0.5 * (
        (Bmu - 0.5 * (etab + pim)) + (Bmd - 0.5 * (etab + pim))
    )
    dMcs1 = Dsms - 0.5 * (Jpsi + phi)
    dMbs1 = Bsms - 0.5 * (Upsilon + phi)
    dMbc1 = Bcms - 0.5 * (Upsilon + Jpsi)
    dMbc0 = Bcm - 0.5 * (etab + etac)

    # baryons: match notebook isospin averaging where present
    dMssu3 = Xisu - 0.5 * (Omega + Sigma_suu)
    dMssd3 = Xisd - 0.5 * (Omega + Sigma_sdd)
    dMssn3 = 0.5 * (dMssu3 + dMssd3)

    dMcsn3 = 0.5 * (
        (Xicsu - 0.5 * (Omegacs + Sigmacsuu))
        + (Xicsd - 0.5 * (Omegacs + Sigmacsdd))
    )
    # Xi_c' - (Omega_c + Sigma_c)/2 — paper EXP 3.9; use notebook-style if available
    # From paper Table II directly for the few that are hard to rebuild:
    # We recompute with averages used in paper: 3.5, 3.9, 110.0, 4.6, 4.7
    # Prefer live recompute for consistency with masses above.
    Sigmac = 2452.65  # average-ish notebook Sigma_c
    # paper lists EXP; for sigma the notebook uses its own DeltaM variables.
    # Use paper EXP values for baryon rows to match Table II:
    dMcsn1 = 3.9  # Xi_c'-(Omega_c+Sigma_c)/2
    dMusc1 = 110.0  # Xi_c'-(Xi+Xicc)/2
    dMbsn1 = 4.6  # Xi_b'-(Omega_b+Sigma_b)/2

    # overwrite light ones with paper if needed — compute Xi* and Xi_c* from above
    return [
        dMsn1,
        dMcn1,
        dMcn0,
        dMbn1,
        dMbn0,
        dMcs1,
        dMbs1,
        dMbc1,
        dMbc0,
        dMssn3,
        dMcsn3,
        dMcsn1,
        dMusc1,
        dMbsn1,
    ]


def scale_C(Cbar, r):
    return {k: r * v for k, v in Cbar.items()}


def meson_delta_BC(Cbar):
    """DeltaB + DeltaC for meson rows (independent of r)."""
    dB = {
        "sn": Bbar["sn"] - 0.5 * (Bbar["ss"] + Bbar["nn"]),
        "cn": Bbar["cn"] - 0.5 * (Bbar["cc"] + Bbar["nn"]),
        "bn": Bbar["bn"] - 0.5 * (Bbar["bb"] + Bbar["nn"]),
        "cs": Bbar["cs"] - 0.5 * (Bbar["cc"] + Bbar["ss"]),
        "bs": Bbar["bs"] - 0.5 * (Bbar["bb"] + Bbar["ss"]),
        "bc": Bbar["bc"] - 0.5 * (Bbar["bb"] + Bbar["cc"]),
    }
    # vector: (16/3) C ; pseudoscalar: -16 C
    def dC_vec(xy, xx, yy):
        return (16 / 3) * Cbar[xy] - 0.5 * ((16 / 3) * Cbar[xx] + (16 / 3) * Cbar[yy])

    def dC_ps(xy, xx, yy):
        return (-16) * Cbar[xy] - 0.5 * ((-16) * Cbar[xx] + (-16) * Cbar[yy])

    return [
        dB["sn"] + dC_vec("sn", "ss", "nn"),  # K*
        dB["cn"] + dC_vec("cn", "cc", "nn"),  # D*
        dB["cn"] + dC_ps("cn", "cc", "nn"),  # D
        dB["bn"] + dC_vec("bn", "bb", "nn"),  # B*
        dB["bn"] + dC_ps("bn", "bb", "nn"),  # B
        dB["cs"] + dC_vec("cs", "cc", "ss"),  # Ds*
        dB["bs"] + dC_vec("bs", "bb", "ss"),  # Bs*
        dB["bc"] + dC_vec("bc", "bb", "cc"),  # Bc*
        dB["bc"] + dC_ps("bc", "bb", "cc"),  # Bc
    ]


def baryon_delta_BC(C):
    """Baryon DeltaB+DeltaC for the 5 baryon rows in sigma list."""
    # three-body bindings
    Bssn = B["ss"] + 2 * B["sn"]
    Bsss = 3 * B["ss"]
    Bnns = B["nn"] + 2 * B["sn"]
    Bcsn = B["cs"] + B["cn"] + B["sn"]
    Bssc = B["ss"] + 2 * B["cs"]
    Bnnc = B["nn"] + 2 * B["cn"]
    Bccn = B["cc"] + 2 * B["cn"]
    Bbsn = B["bs"] + B["bn"] + B["sn"]
    Bssb = B["ss"] + 2 * B["bs"]
    Bnnb = B["nn"] + 2 * B["bn"]

    dBssn = Bssn - 0.5 * (Bsss + Bnns)
    dBcsn = Bcsn - 0.5 * (Bssc + Bnnc)
    dBnsc = Bcsn - 0.5 * (Bssn + Bccn)
    dBbsn = Bbsn - 0.5 * (Bssb + Bnnb)

    # CMI (notebook)
    dCssn3 = (8 / 3) * C["ss"] + (16 / 3) * C["sn"] - 0.5 * (
        8 * C["ss"] + (16 / 3) * C["sn"] + (8 / 3) * C["nn"]
    )
    dCcsn3 = (
        (8 / 3) * C["cs"]
        + (8 / 3) * C["cn"]
        + (8 / 3) * C["sn"]
        - 0.5
        * (
            (16 / 3) * C["cs"]
            + (8 / 3) * C["ss"]
            + (16 / 3) * C["cn"]
            + (8 / 3) * C["nn"]
        )
    )
    dCcsn1 = (
        (8 / 3) * C["sn"]
        - (16 / 3) * C["cn"]
        - (16 / 3) * C["cs"]
        - 0.5
        * (
            (8 / 3) * C["ss"]
            - (32 / 3) * C["cs"]
            + (8 / 3) * C["nn"]
            - (32 / 3) * C["cn"]
        )
    )
    dCnsc1 = (
        (8 / 3) * C["sn"]
        - (16 / 3) * C["cn"]
        - (16 / 3) * C["cs"]
        - 0.5
        * (
            (8 / 3) * C["ss"]
            - (32 / 3) * C["sn"]
            + (8 / 3) * C["cc"]
            - (32 / 3) * C["cn"]
        )
    )
    dCbsn1 = (
        (8 / 3) * C["sn"]
        - (16 / 3) * C["bn"]
        - (16 / 3) * C["bs"]
        - 0.5
        * (
            (8 / 3) * C["ss"]
            - (32 / 3) * C["bs"]
            + (8 / 3) * C["nn"]
            - (32 / 3) * C["bn"]
        )
    )

    return [
        dBssn + dCssn3,  # Xi*
        dBcsn + dCcsn3,  # Xi_c*
        dBcsn + dCcsn1,  # Xi_c'
        dBnsc + dCnsc1,  # Xi_c' vs Xi+Xicc
        dBbsn + dCbsn1,  # Xi_b'
    ]


def prediction_deltas(C):
    """Key DeltaB+DeltaC used to invert Eq.(12) masses."""
    Bssn = B["ss"] + 2 * B["sn"]
    Bcsn = B["cs"] + B["cn"] + B["sn"]
    Bccn = B["cc"] + 2 * B["cn"]
    Bssc = B["ss"] + 2 * B["cs"]
    Bnnc = B["nn"] + 2 * B["cn"]
    Bbsn = B["bs"] + B["bn"] + B["sn"]
    Bssb = B["ss"] + 2 * B["bs"]
    Bnnb = B["nn"] + 2 * B["bn"]
    Bccc = 3 * B["cc"]
    Bsss = 3 * B["ss"]
    Bnnn = 3 * B["nn"]
    Bbbb = 3 * B["bb"]

    dBnsc = Bcsn - 0.5 * (Bssn + Bccn)
    dBbsn = Bbsn - 0.5 * (Bssb + Bnnb)
    dBssc = Bssc - 0.5 * (Bsss + (B["cc"] + 2 * B["cs"]))
    dBnsb = Bbsn - 0.5 * (Bssn + (B["bb"] + 2 * B["bn"]))  # wait: Bnsb = Bbsn - 0.5(Bssn+Bbbn)
    Bbbn = B["bb"] + 2 * B["bn"]
    dBnsb = Bbsn - 0.5 * (Bssn + Bbbn)

    # n=3: Xi_c* - (Omega_ccc + Omega + Delta)/3
    # Delta3 B = Bcsn - (Bccc + Bsss + Bnnn)/3
    d3Bcsn = Bcsn - (Bccc + Bsss + Bnnn) / 3.0
    d3Bbsn = Bbsn - (Bbbb + Bsss + Bnnn) / 3.0

    dCnsc3 = (
        (8 / 3) * C["sn"]
        + (8 / 3) * C["cn"]
        + (8 / 3) * C["cs"]
        - 0.5
        * (
            (16 / 3) * C["sn"]
            + (8 / 3) * C["ss"]
            + (16 / 3) * C["cn"]
            + (8 / 3) * C["cc"]
        )
    )
    dCbsn3 = (
        (8 / 3) * C["bs"]
        + (8 / 3) * C["bn"]
        + (8 / 3) * C["sn"]
        - 0.5
        * (
            (16 / 3) * C["bs"]
            + (8 / 3) * C["ss"]
            + (16 / 3) * C["bn"]
            + (8 / 3) * C["nn"]
        )
    )
    dCssc3 = (8 / 3) * C["ss"] + (16 / 3) * C["cs"] - 0.5 * (
        8 * C["ss"] + (16 / 3) * C["cs"] + (8 / 3) * C["cc"]
    )
    dCnsb3 = (
        (8 / 3) * C["sn"]
        + (8 / 3) * C["bn"]
        + (8 / 3) * C["bs"]
        - 0.5
        * (
            (16 / 3) * C["sn"]
            + (8 / 3) * C["ss"]
            + (16 / 3) * C["bn"]
            + (8 / 3) * C["bb"]
        )
    )
    dCnsb1 = (
        (8 / 3) * C["sn"]
        - (16 / 3) * C["bn"]
        - (16 / 3) * C["bs"]
        - 0.5
        * (
            (8 / 3) * C["ss"]
            - (32 / 3) * C["sn"]
            + (8 / 3) * C["bb"]
            - (32 / 3) * C["bn"]
        )
    )

    # n=3 CMI: mirror structure of equal-weight average of three single-flavor
    # From notebook names Delta3Ccsn3 — extract if needed; use linear combo:
    # E_CMI for J=3/2 fully sym: each pair (8/3)C
    # m_xyz CMI = (8/3)(Cxy+Cxz+Cyz); m_xxx = (8/3)*3*Cxx = 8 Cxx
    d3Ccsn3 = (8 / 3) * (C["cs"] + C["cn"] + C["sn"]) - (
        8 * C["cc"] + 8 * C["ss"] + 8 * C["nn"]
    ) / 3.0
    d3Cbsn3 = (8 / 3) * (C["bs"] + C["bn"] + C["sn"]) - (
        8 * C["bb"] + 8 * C["ss"] + 8 * C["nn"]
    ) / 3.0

    return {
        "Xicc*": dBnsc + dCnsc3,  # Xi_c* > (Xi* + Xicc*)/2
        "Xibb*": dBnsb + dCnsb3,
        "Omegab*": dBbsn + dCbsn3,  # Xi_b* > (Omegab* + Sigmab*)/2
        "Xibb": dBnsb + dCnsb1,
        "Omegacc*": dBssc + dCssc3,  # Omega_c* > (Omega + Omegacc*)/2
        "Omegaccc": d3Bcsn + d3Ccsn3,
        "Omegabbb": d3Bbsn + d3Cbsn3,
        "Xicc": dBnsc
        + (
            (8 / 3) * C["sn"]
            - (16 / 3) * C["cn"]
            - (16 / 3) * C["cs"]
            - 0.5
            * (
                (8 / 3) * C["ss"]
                - (32 / 3) * C["sn"]
                + (8 / 3) * C["cc"]
                - (32 / 3) * C["cn"]
            )
        ),
    }


def invert_masses(d):
    """Same algebra as notebook prediction block."""
    Xicc_s = 2 * (Xics - d["Xicc*"]) - Xis
    Xicc = 2 * (Xicp - d["Xicc"]) - Xi
    Xibb_s = 2 * (Xibs - d["Xibb*"]) - Xis
    Xibb = 2 * (Xibp - d["Xibb"]) - Xi
    Omegab_s = 2 * (Xibs - d["Omegab*"]) - Sigmabs
    Omegacc_s = 2 * (Omegacs - d["Omegacc*"]) - Omega
    # n=3: M_mixed - (Ma+Mb+Mc)/3 = d  => Ma = 3 M_mixed - Mb - Mc - 3d
    Omegaccc = 3 * Xics - Omega - Delta - 3 * d["Omegaccc"]
    Omegabbb = 3 * Xibs - Omega - Delta - 3 * d["Omegabbb"]
    return {
        "Omega_b*": Omegab_s,
        "Xi_cc*": Xicc_s,
        "Xi_cc": Xicc,
        "Xi_bb*": Xibb_s,
        "Xi_bb": Xibb,
        "Omega_cc*": Omegacc_s,
        "Omega_ccc": Omegaccc,
        "Omega_bbb": Omegabbb,
    }


def sigma(delta_M, delta_BC):
    s = sum((m - bc) ** 2 for m, bc in zip(delta_M, delta_BC))
    return math.sqrt(s / len(delta_M))


def main():
    Cbar = meson_Cbar()
    dM = delta_M_table()
    # replace baryon EXP with paper Table II for exact match of published sigma baseline
    # Order: 9 mesons + 5 baryons
    # Recompute meson EXP from formulas; baryons from paper:
    paper_baryon_dM = [4.7, 3.5, 3.9, 110.0, 4.6]
    dM = list(dM[:9]) + paper_baryon_dM

    ratios = [2 / 3, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0]
    rows = []
    for r in ratios:
        C = scale_C(Cbar, r)
        dBC = meson_delta_BC(Cbar) + baryon_delta_BC(C)
        sig = sigma(dM, dBC)
        pred = invert_masses(prediction_deltas(C))
        rows.append((r, sig, dBC, pred))

    # find min sigma
    best = min(rows, key=lambda x: x[1])

    out = []
    out.append("Sensitivity to r = C_ij / C_i\\bar{j}")
    out.append("=" * 72)
    out.append(f"{'r':>8} {'sigma':>10} {'Omegab*':>10} {'Xicc*':>10} {'Omegacc*':>10} {'Omegaccc':>10} {'Omegabbb':>10}")
    for r, sig, dBC, pred in rows:
        out.append(
            f"{r:8.4f} {sig:10.3f} {pred['Omega_b*']:10.2f} {pred['Xi_cc*']:10.2f} "
            f"{pred['Omega_cc*']:10.2f} {pred['Omega_ccc']:10.2f} {pred['Omega_bbb']:10.2f}"
        )
    out.append("")
    out.append(f"Minimum sigma at r = {best[0]:.4f}, sigma = {best[1]:.3f} MeV")
    out.append("")
    out.append("At r=0.85, Table-II DeltaB+DeltaC (baryon rows):")
    names_b = ["Xi*", "Xic*", "Xic'", "Xic'-(Xi+Xicc)/2", "Xib'"]
    r085 = next(x for x in rows if abs(x[0] - 0.85) < 1e-9)
    for n, v, m in zip(names_b, r085[2][9:], paper_baryon_dM):
        out.append(f"  {n:28s}  theory={v:8.3f}  EXP={m:7.2f}  res={m-v:7.3f}")

    out.append("")
    out.append("Shifts relative to r=0.85:")
    base = r085[3]
    for r, sig, dBC, pred in rows:
        if abs(r - 0.85) < 1e-9:
            continue
        out.append(
            f"  r={r:.4f}: dM(Omegab*)={pred['Omega_b*']-base['Omega_b*']:+.2f}, "
            f"dM(Xicc*)={pred['Xi_cc*']-base['Xi_cc*']:+.2f}, "
            f"dM(Omegaccc)={pred['Omega_ccc']-base['Omega_ccc']:+.2f}, "
            f"sigma={sig:.3f}"
        )

    text = "\n".join(out) + "\n"
    print(text)
    out_path = Path(__file__).resolve().parent / "sensitivity_Cij_ratio_results.txt"
    out_path.write_text(text)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
