from dataclasses import dataclass
import numpy as np


@dataclass
class AppDefaults:
    pulse_sequence: str = "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR36020kHz_30B1_rf10kHz/pulse1000.bruker"
    maximum_amplitude_hz: float = 10000.0
    vector_length: float = 1.0
    time_per_pulse_s: float = 0.5
    inpo_fact: int = 5
    x_expand: float = 0.0
    offset_hz: float = 0.0
    scaling_percent: float = 100.0
    calculation_method: int = 1
    calculation_language: str = "Rust"

    export_directory: str = "/Users/leon/Desktop/Physik/Glaser"

    quality_input_directory: str = "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR_ohne_B1_robustness_20_kHz/UR360_20kHz_noB1_rf10kHz(new_loop_sorted)"
    quality_output_directory: str = "/Users/leon/Desktop/Physik/Glaser/Bachelor_Thesis/images/Quality images"

    initial_vector: np.ndarray = None


def get_default_initial_vector():
    return np.array([0, 0, 1])


DEFAULTS = AppDefaults()