import time
import numpy as np
import pandas as pd

from atmosphericRadiationDoseAndFlux import doseAndFluxCalculator as DAF
from atmosphericRadiationDoseAndFlux import particleResponse, units, settings, responseFileParameters


ENERGY_BINS = 10 ** (0.1 * (np.array(range(1, 52)) - 1) + 1)
FLUXES = np.full(50, 1.0)
ALTITUDE = [12.0]
OUTPUT_COLUMNS = ["edose", "adose", "dosee", "tn1", "tn2", "tn3", "SEU", "SEL"]


def legacy_full(input_energy_bins, input_fluxes_mev, altitudes_km, particle_name="proton"):
    particle_for_calculations, input_energy_bins_array, altitudes_km_array, input_fluxes_array = settings.formatInputVariables(
        input_energy_bins, input_fluxes_mev, altitudes_km, particle_name, 0.0
    )
    input_energy_differences = input_energy_bins_array[1:] - input_energy_bins_array[:-1]
    input_fluxes_integrated = DAF.calculate_energy_integrated_flux(input_fluxes_array, input_energy_differences)

    output_df = pd.DataFrame({"altitude (km)": altitudes_km_array})
    for dose_type in particleResponse.fullListOfDoseResponseTypes:
        dose_response_for_particle = particleResponse.fullDoseResponseDict[dose_type](particle_for_calculations, dose_type)
        coordinate_doses = []
        for altitude_km in altitudes_km_array:
            altitude = units.Distance(altitude_km * 1000.0)
            total_dose = dose_response_for_particle.calculateDose(altitude, input_energy_bins_array, input_fluxes_integrated)
            coordinate_doses.append(total_dose)
        output_df[dose_type] = coordinate_doses

    output_df["SEU"] = output_df["tn2"] * 1e-13
    output_df["SEL"] = output_df["tn2"] * 1e-8
    return output_df


def strategy1_cached_responses(input_energy_bins, input_fluxes_mev, altitudes_km, particle_name="proton"):
    particle_for_calculations, input_energy_bins_array, altitudes_km_array, input_fluxes_array = settings.formatInputVariables(
        input_energy_bins, input_fluxes_mev, altitudes_km, particle_name, 0.0
    )
    input_energy_differences = input_energy_bins_array[1:] - input_energy_bins_array[:-1]
    input_fluxes_integrated = np.asarray(
        DAF.calculate_energy_integrated_flux(input_fluxes_array, input_energy_differences), dtype=np.float64
    )
    edose_response, adose_response, dosee_response, neutron_response = DAF._get_cached_response_arrays(
        particle_for_calculations.particleName
    )

    altitude = units.Distance(float(altitudes_km_array[0]) * 1000.0)
    altitude_layer_index, f1 = responseFileParameters.calculate_altitude_layer_params(altitude.meters, altitude.km)
    altitude_index_above = min(137, altitude_layer_index + 1)

    output_df = pd.DataFrame(
        {
            "altitude (km)": altitudes_km_array,
            "edose": [particleResponse.calculate_dose_response(input_fluxes_integrated, edose_response, altitude_layer_index, altitude_index_above, f1)],
            "adose": [particleResponse.calculate_dose_response(input_fluxes_integrated, adose_response, altitude_layer_index, altitude_index_above, f1)],
            "dosee": [particleResponse.calculate_dose_response(input_fluxes_integrated, dosee_response, altitude_layer_index, altitude_index_above, f1)],
            "tn1": [particleResponse.calculate_neutron_response(input_fluxes_integrated, neutron_response, altitude_layer_index, altitude_index_above, f1, 0)],
            "tn2": [particleResponse.calculate_neutron_response(input_fluxes_integrated, neutron_response, altitude_layer_index, altitude_index_above, f1, 50)],
            "tn3": [particleResponse.calculate_neutron_response(input_fluxes_integrated, neutron_response, altitude_layer_index, altitude_index_above, f1, 100)],
        }
    )
    output_df["SEU"] = output_df["tn2"] * 1e-13
    output_df["SEL"] = output_df["tn2"] * 1e-8
    return output_df


def strategy2_current_optimized(input_energy_bins, input_fluxes_mev, altitudes_km, particle_name="proton"):
    return DAF.calculate_from_energy_spec_array(
        input_energy_bins, input_fluxes_mev, altitudes_km, particleName=particle_name
    )


def benchmark(label, function, runs):
    function(ENERGY_BINS, FLUXES, ALTITUDE, "proton")
    start = time.perf_counter()
    for _ in range(runs):
        function(ENERGY_BINS, FLUXES, ALTITUDE, "proton")
    elapsed = time.perf_counter() - start
    milliseconds_per_call = (elapsed * 1000.0) / runs
    print(f"{label}: {milliseconds_per_call:.6f} ms/call over {runs} runs")
    return milliseconds_per_call


def verify_equivalence():
    legacy_output = legacy_full(ENERGY_BINS, FLUXES, ALTITUDE, "proton")
    optimized_output = strategy2_current_optimized(ENERGY_BINS, FLUXES, ALTITUDE, "proton")
    if not np.array_equal(legacy_output[OUTPUT_COLUMNS].to_numpy(), optimized_output[OUTPUT_COLUMNS].to_numpy()):
        differences = np.abs(legacy_output[OUTPUT_COLUMNS].to_numpy() - optimized_output[OUTPUT_COLUMNS].to_numpy())
        raise RuntimeError(f"Output mismatch detected; max abs diff={differences.max()}")
    print("equivalence_check: exact match")


if __name__ == "__main__":
    verify_equivalence()

    legacy_ms = benchmark("legacy_full", legacy_full, runs=1000)
    strategy1_ms = benchmark("strategy1_cached_responses", strategy1_cached_responses, runs=1000)
    strategy2_ms = benchmark("strategy2_current_optimized", strategy2_current_optimized, runs=5000)

    print(f"speedup_strategy1_vs_legacy: {legacy_ms / strategy1_ms:.2f}x")
    print(f"speedup_strategy2_vs_strategy1: {strategy1_ms / strategy2_ms:.2f}x")
    print(f"speedup_total_vs_legacy: {legacy_ms / strategy2_ms:.2f}x")
