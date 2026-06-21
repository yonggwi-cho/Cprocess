"""
Parameter sweep: ion implant energy and dose.

Validates fundamental BCA physics:
- Rp increases monotonically with energy
- dRp < Rp always
- Deposited fraction >= 30%
- Linear dose scaling of peak concentration
"""

import sys
import numpy as np
sys.path.insert(0, "python")
import cprocess as cp

IONS = 20000  # MC ions per run (speed vs. accuracy trade-off)

# ------------------------------------------------------------------
# Helper: run a single MC implant and return key metrics
# ------------------------------------------------------------------

def run_implant(species, energy_kev, dose, ions=IONS, seed=42):
    """Return dict with rp_nm, drp_nm, deposited, peak_conc, peak_z_nm."""
    sim = cp.Simulation()
    # 1 µm x 1 µm x 0.5 µm box, coarse grid for speed
    # Depth: 1 µm for low energies, deeper for high energies to avoid ions exiting
    depth_um = max(1.0, energy_kev * 0.01)  # ~10 nm/keV headroom
    sim.mesh(x=1.0, y=1.0, z=depth_um, nx=8, ny=8, nz=20)
    sim.region("silicon")
    result = sim.implant(species, dose=dose, energy=energy_kev,
                         mc=True, ions=ions, seed=seed)

    # Peak concentration and its depth
    conc = sim.field(species)
    cents = sim.cell_centroids  # µm, shape (n_cells, 3)
    z_um = cents[:, 2]          # depth axis

    if conc.max() > 0:
        idx_peak = int(np.argmax(conc))
        peak_conc = float(conc[idx_peak])
        peak_z_nm = float(z_um[idx_peak]) * 1e3  # µm -> nm
    else:
        peak_conc = 0.0
        peak_z_nm = 0.0

    return {
        "rp_nm":     result.rp_nm,
        "drp_nm":    result.drp_nm,
        "deposited": result.deposited,
        "peak_conc": peak_conc,
        "peak_z_nm": peak_z_nm,
    }


# ------------------------------------------------------------------
# Energy sweep
# ------------------------------------------------------------------

ENERGY_SPECIES = [
    ("B",  [20, 50, 100, 150]),
    ("P",  [30, 60, 100, 150]),
    ("As", [50, 100, 150, 200]),
]
FIXED_DOSE = 1e14

def energy_sweep():
    """Run energy sweep for B, P, As and return list of (row, pass) tuples."""
    rows = []
    print("\n=== Energy Sweep (dose=1e14 cm^-2, ions={}) ===".format(IONS))
    header = f"{'Species':<8}{'E(keV)':>8}{'Rp(nm)':>10}{'dRp(nm)':>10}{'Dep%':>8}{'PkZ(nm)':>10}  Checks"
    print(header)
    print("-" * len(header))

    all_pass = True
    # Collect per-species results for monotonicity check
    for species, energies in ENERGY_SPECIES:
        rp_list = []
        for energy in energies:
            m = run_implant(species, energy, FIXED_DOSE)
            rp_list.append(m["rp_nm"])

            dep_frac = m["deposited"] / IONS

            checks = []
            checks.append(("Rp>0",        m["rp_nm"] > 0))
            checks.append(("dRp<Rp",      m["drp_nm"] < m["rp_nm"]))
            checks.append(("dep>=30%",    dep_frac >= 0.30))

            row_pass = all(v for _, v in checks)
            if not row_pass:
                all_pass = False

            status = "PASS" if row_pass else "FAIL"
            fail_tags = " ".join(k for k, v in checks if not v)
            check_str = status if row_pass else f"FAIL [{fail_tags}]"

            print(f"{species:<8}{energy:>8}{m['rp_nm']:>10.1f}{m['drp_nm']:>10.1f}"
                  f"{dep_frac*100:>7.1f}%{m['peak_z_nm']:>10.1f}  {check_str}")

            rows.append((species, energy, m, row_pass))

        # Monotonicity across energies for this species
        mono_ok = all(rp_list[i+1] > rp_list[i] for i in range(len(rp_list)-1))
        tag = "PASS" if mono_ok else "FAIL"
        if not mono_ok:
            all_pass = False
        print(f"  >> {species} Rp monotonic with energy: {tag}  "
              f"({', '.join(f'{r:.0f}' for r in rp_list)} nm)")

    return all_pass


# ------------------------------------------------------------------
# Dose sweep — linear scaling check
# ------------------------------------------------------------------

DOSE_SPECIES = "B"
DOSE_ENERGY  = 50       # keV
DOSES        = [1e12, 1e13, 1e14, 1e15]

def dose_sweep():
    """Run dose sweep and verify linear peak-concentration scaling."""
    print(f"\n=== Dose Sweep ({DOSE_SPECIES} @ {DOSE_ENERGY} keV, ions={IONS}) ===")
    header = f"{'Dose(cm-2)':<14}{'Rp(nm)':>10}{'PkConc(cm-3)':>16}{'PkConc ratio':>14}  Checks"
    print(header)
    print("-" * len(header))

    results = []
    for dose in DOSES:
        m = run_implant(DOSE_SPECIES, DOSE_ENERGY, dose, seed=99)
        results.append((dose, m))

    all_pass = True
    prev_peak = None
    prev_dose = None
    for dose, m in results:
        ratio_str = "  —  "
        ratio_ok  = True
        if prev_peak is not None and prev_peak > 0 and m["peak_conc"] > 0:
            expected_ratio = dose / prev_dose  # should be ~10
            actual_ratio   = m["peak_conc"] / prev_peak
            # Accept anything between 0.5x and 2x of expected (wide tolerance for MC noise)
            ratio_ok = (actual_ratio > expected_ratio * 0.5) and \
                       (actual_ratio < expected_ratio * 2.0)
            ratio_str = f"{actual_ratio:>10.2f}x"
            if not ratio_ok:
                all_pass = False

        status = "PASS" if ratio_ok else "FAIL"
        print(f"{dose:<14.2e}{m['rp_nm']:>10.1f}{m['peak_conc']:>16.3e}{ratio_str:>14}  {status}")

        prev_peak = m["peak_conc"]
        prev_dose = dose

    return all_pass


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------

def main():
    ok1 = energy_sweep()
    ok2 = dose_sweep()

    print("\n" + "=" * 60)
    if ok1 and ok2:
        print("ALL CHECKS PASSED")
        return 0
    else:
        print("SOME CHECKS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
