#!/usr/bin/env python3
"""Regenerate KGO fixtures for tests. Run from repo root: python3 scripts/generate_kgo_reference.py"""
import numpy as np
from pathlib import Path

from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import calculate_from_energy_spec_array

FIXTURES = Path(__file__).resolve().parent.parent / "tests" / "fixtures"
DOSE_COLS = ["edose", "adose", "dosee", "tn1", "tn2", "tn3", "SEU", "SEL"]


def default_bins_and_altitudes():
    input_energy_bins = 10 ** (0.1 * (np.array(range(1, 52)) - 1) + 1)
    altitudes = np.linspace(0, 60, 60)
    return input_energy_bins, altitudes


def df_to_npz_dict(df):
    out = {"altitude_km": df["altitude (km)"].to_numpy(dtype=np.float64)}
    for c in DOSE_COLS:
        out[c] = df[c].to_numpy(dtype=np.float64)
    return out


def main():
    FIXTURES.mkdir(parents=True, exist_ok=True)
    bins, alts = default_bins_and_altitudes()

    scenarios = [
        ("kgo_energy_flat_proton", np.full(50, 1.0), "proton"),
        ("kgo_energy_flat_alpha", np.full(50, 1.0), "alpha"),
        ("kgo_energy_flat_both", np.full(50, 1.0), "both"),
        ("kgo_energy_singlestep_proton", np.append(1.0, np.full(49, 0.0)), "proton"),
        ("kgo_energy_singlestep_alpha", np.append(1.0, np.full(49, 0.0)), "alpha"),
    ]

    for name, fluxes, species in scenarios:
        df = calculate_from_energy_spec_array(bins, fluxes, alts, particleName=species)
        path = FIXTURES / f"{name}.npz"
        np.savez_compressed(path, **df_to_npz_dict(df))
        print(f"Wrote {path}")

    # Non-flat energy spectra (differential flux varies across bins)
    energy_mid = (bins[1:] + bins[:-1]) / 2.0
    flux_e_minus7 = energy_mid ** (-7.0)
    flux_e_minus2 = energy_mid ** (-2.0)
    for species in ("proton", "alpha", "both"):
        name = f"kgo_energy_powerlaw_em7_{species}"
        df = calculate_from_energy_spec_array(bins, flux_e_minus7, alts, particleName=species)
        path = FIXTURES / f"{name}.npz"
        np.savez_compressed(path, **df_to_npz_dict(df))
        print(f"Wrote {path}")
    df = calculate_from_energy_spec_array(bins, flux_e_minus2, alts, particleName="proton")
    path = FIXTURES / "kgo_energy_powerlaw_em2_proton.npz"
    np.savez_compressed(path, **df_to_npz_dict(df))
    print(f"Wrote {path}")

    # Rigidity: default rigidity bins from energy (same as calculate_from_rigidity_spec)
    from atmosphericRadiationDoseAndFlux.doseAndFluxCalculator import (
        calculate_from_rigidity_spec_array,
    )
    from atmosphericRadiationDoseAndFlux import particle as particle_mod

    p = particle_mod.Particle("proton")
    import ParticleRigidityCalculationTools as PRCT

    rigidity_bins = np.array(
        PRCT.convertParticleEnergyToRigidity(
            bins, particleMassAU=p.atomicMass, particleChargeAU=p.atomicCharge
        )
    )
    rigidity_fluxes = np.full(50, 1.0)
    df_r = calculate_from_rigidity_spec_array(
        rigidity_bins, rigidity_fluxes, [60.0], particleName="proton"
    )
    path_r = FIXTURES / "kgo_rigidity_flat_proton_60km.npz"
    np.savez_compressed(path_r, **df_to_npz_dict(df_r))
    print(f"Wrote {path_r}")

    # Non-flat rigidity spectrum j(R) ∝ R^-7 at bin midpoints (full altitude grid)
    for species in ("proton", "alpha"):
        p = particle_mod.Particle(species)
        rb = np.array(
            PRCT.convertParticleEnergyToRigidity(
                bins, particleMassAU=p.atomicMass, particleChargeAU=p.atomicCharge
            )
        )
        rmid = (rb[1:] + rb[:-1]) / 2.0
        flux_rm7 = rmid ** (-7.0)
        df_nr = calculate_from_rigidity_spec_array(rb, flux_rm7, alts, particleName=species)
        path_nr = FIXTURES / f"kgo_rigidity_powerlaw_rm7_{species}.npz"
        np.savez_compressed(path_nr, **df_to_npz_dict(df_nr))
        print(f"Wrote {path_nr}")


if __name__ == "__main__":
    main()
