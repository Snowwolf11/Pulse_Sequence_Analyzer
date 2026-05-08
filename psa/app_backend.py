# psa/app_backend.py

from dataclasses import dataclass
import numpy as np
import os
import mpld3


from psa.createCurve import createCurve
from psa.calculate_curveData import (
    angle_between_vectors,
    angle_with_x_axis,
    number_of_pulses,
    calculate_closed_curve_area_app,
    curvature_torsion_3d,
)
from psa.calculate_curve_stability import calculate_curve_stability
from psa.calculate_pulse_sequence_quality_images import calculate_pulse_sequence_quality_images


@dataclass
class AnalyzerParams:
    pulse_sequence: str
    pulse_amplitude_khz: float
    vector_length: float
    time_per_pulse_us: float
    inpo_fact: int
    x_expand: float
    offset_khz: float
    scaling_percent: float
    calculation_method: int
    initial_vector: np.ndarray
    calculation_language: str


@dataclass
class AnalysisResult:
    CM: np.ndarray
    VM: np.ndarray
    PS_mat: np.ndarray

    totalRot: float
    phi: float
    numberOfPulses: int

    Axy: float
    Axz: float
    Ayz: float

    arc_length: np.ndarray
    curvature: np.ndarray
    torsion: np.ndarray

    integrated_curvature: float
    integrated_torsion: float
    integrated_absolut_torsion: float
    avg_curvature: float
    avg_torsion: float

    maximum_amplitude_hz: float
    offset_hz: float
    time_per_pulse_s: float


def analyze_pulse_sequence(params: AnalyzerParams) -> AnalysisResult:
    maximum_amplitude_hz = (
        params.pulse_amplitude_khz
        * params.scaling_percent
        * 1000
        / 100
    )

    offset_hz = params.offset_khz * 1000
    time_per_pulse_s = params.time_per_pulse_us * 1e-6

    CM, VM, PS_mat = createCurve(
        PulseSequence=params.pulse_sequence,
        T=time_per_pulse_s,
        l=params.vector_length,
        maximumAmplitude=maximum_amplitude_hz,
        offset=offset_hz,
        inpoFact=params.inpo_fact,
        xExpand=params.x_expand,
        calculationMethod=params.calculation_method,
        initialVector=params.initial_vector,
        language=params.calculation_language,
    )

    totalRot = angle_between_vectors(VM[0, :], VM[-1, :])
    phi = angle_with_x_axis(VM[-1, :])
    numberOfPulses = number_of_pulses(PS_mat, params.inpo_fact)

    Axy = calculate_closed_curve_area_app(CM[:, [0, 1]], close_curve=False)
    Axz = calculate_closed_curve_area_app(CM[:, [0, 2]], close_curve=False)
    Ayz = calculate_closed_curve_area_app(CM[:, [1, 2]], close_curve=False)

    (
        arc_length,
        curvature,
        torsion,
        integrated_curvature,
        integrated_torsion,
        integrated_absolut_torsion,
        avg_curvature,
        avg_torsion,
    ) = curvature_torsion_3d(CM)

    return AnalysisResult(
        CM=CM,
        VM=VM,
        PS_mat=PS_mat,

        totalRot=totalRot,
        phi=phi,
        numberOfPulses=numberOfPulses,

        Axy=Axy,
        Axz=Axz,
        Ayz=Ayz,

        arc_length=arc_length,
        curvature=curvature,
        torsion=torsion,

        integrated_curvature=integrated_curvature,
        integrated_torsion=integrated_torsion,
        integrated_absolut_torsion=integrated_absolut_torsion,
        avg_curvature=avg_curvature,
        avg_torsion=avg_torsion,

        maximum_amplitude_hz=maximum_amplitude_hz,
        offset_hz=offset_hz,
        time_per_pulse_s=time_per_pulse_s,
    )


@dataclass
class StabilityParams:
    PS_mat: np.ndarray
    time_per_pulse_s: float
    vector_length: float
    maximum_amplitude_hz: float
    scaling_range_percent: float
    offset_range_khz: float
    stability_calculation_method: int
    initial_vector: np.ndarray
    calculation_language: str


def calculate_stability_backend(params: StabilityParams):
    return calculate_curve_stability(
        params.PS_mat,
        params.time_per_pulse_s,
        params.vector_length,
        params.maximum_amplitude_hz,
        scalingRange_percent=params.scaling_range_percent,
        offsetRange_kHz=params.offset_range_khz,
        stabilityCalculationMethod=params.stability_calculation_method,
        initialVector=params.initial_vector,
        language=params.calculation_language,
    )


@dataclass
class QualityImageParams:
    input_directory: str
    range_hz: float
    time_per_pulse_s: float
    maximum_amplitude_hz: float
    initial_vector: np.ndarray
    calculation_type: int
    changing_variable: int
    resolution: float
    calculation_language: str
    progress_callback: object = None

def calculate_quality_image_backend(params: QualityImageParams):
    return calculate_pulse_sequence_quality_images(
        dirname=params.input_directory,
        Range=params.range_hz,
        T=params.time_per_pulse_s,
        Umax=params.maximum_amplitude_hz,
        initialVector=params.initial_vector,
        calcType=params.calculation_type,
        changingVariable=params.changing_variable,
        Resolution=params.resolution,
        language=params.calculation_language,
        progress_callback=params.progress_callback,
    )

def build_curve_data_export_string(
    pulse_sequence,
    totalRot,
    time_per_pulse_s,
    vector_length,
    maximum_amplitude_hz,
    offset_hz,
    inpo_fact,
    x_expand,
    calculation_method,
    initial_vector,
    phi,
    numberOfPulses,
    Axy,
    Axz,
    Ayz,
    integrated_curvature,
    integrated_torsion,
    integrated_absolut_torsion,
):
    return (
        f"Pulse Sequence {pulse_sequence} ({round(1E3 * totalRot) / 1E3}°)\n"
        f"Time per Pulse: {time_per_pulse_s * 10**6} µs\n"
        f"Vector Length: {vector_length}\n"
        f"Maximum Amplitude: {maximum_amplitude_hz / 1000} kHz\n"
        f"Offset: {offset_hz / 1000} kHz\n"
        f"InpoFact: {inpo_fact}\n"
        f"x Expand: {x_expand}\n"
        f"Calculation Method: {calculation_method}\n"
        f"Initial Vector: {initial_vector}\n"
        f"Angle to x-Axis= {phi}°;\n"
        f"number of pulses = {numberOfPulses};\n"
        f"Axy= {Axy};\n"
        f"Axz= {Axz},\n"
        f"Ayz= {Ayz};\n"
        f"Integrated Curvature= {integrated_curvature};\n"
        f"Integrated Torsion= {integrated_torsion};\n"
        f"Integrated abs Torsion= {integrated_absolut_torsion};"
    )


def write_analysis_export(
    output_dir,
    CM,
    VM,
    PS_mat,
    arc_length,
    curvature,
    torsion,
    curve_data_string,
    error_curve_figure=None,
):
    os.mkdir(output_dir)

    np.savetxt(
        os.path.join(output_dir, "Error_Curve.csv"),
        CM,
        delimiter=",",
    )

    np.savetxt(
        os.path.join(output_dir, "Trajectory_Curve.csv"),
        VM,
        delimiter=",",
    )

    np.savetxt(
        os.path.join(output_dir, "Pulse_Sequence.csv"),
        PS_mat,
        delimiter=",",
    )

    np.savetxt(
        os.path.join(output_dir, "curvature-arc_length.csv"),
        np.column_stack((arc_length, curvature)),
        delimiter=",",
    )

    np.savetxt(
        os.path.join(output_dir, "torsion-arc_length.csv"),
        np.column_stack((arc_length, torsion)),
        delimiter=",",
    )

    with open(os.path.join(output_dir, "Curve_data.txt"), "w") as file:
        file.write(curve_data_string)

    if error_curve_figure is not None:
        mpld3.save_html(
            error_curve_figure,
            os.path.join(output_dir, "Error_Curve_plot.html"),
        )