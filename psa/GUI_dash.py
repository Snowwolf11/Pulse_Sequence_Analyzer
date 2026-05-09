#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import traceback
import numpy as np
import matplotlib.pyplot as plt

from dash import Dash, dcc, html, Input, Output, State, no_update
import plotly.graph_objects as go

from psa.app_backend import (
    AnalyzerParams,
    analyze_pulse_sequence,
    StabilityParams,
    calculate_stability_backend,
    QualityImageParams,
    calculate_quality_image_backend,
    build_curve_data_export_string,
    write_analysis_export,
)
from psa.app_config import DEFAULTS, get_default_initial_vector
from psa.documentation import DOCUMENTATION_STR
from psa.pulse_string_to_array import parse_pulse_sequence


APP_STYLE = {
    "fontFamily": "Inter, -apple-system, BlinkMacSystemFont, Segoe UI, sans-serif",
    "background": "#0f172a",
    "color": "#e5e7eb",
    "minHeight": "100vh",
    "padding": "18px",
}

CARD_STYLE = {
    "background": "#111827",
    "border": "1px solid #334155",
    "borderRadius": "14px",
    "padding": "14px",
    "boxShadow": "0 10px 30px rgba(0,0,0,0.25)",
}

INPUT_STYLE = {
    "width": "100%",
    "background": "#020617",
    "color": "#e5e7eb",
    "border": "1px solid #475569",
    "borderRadius": "8px",
    "padding": "7px",
}

BUTTON_STYLE = {
    "background": "#2563eb",
    "color": "white",
    "border": "none",
    "borderRadius": "10px",
    "padding": "10px 14px",
    "cursor": "pointer",
    "fontWeight": "600",
}

SECONDARY_BUTTON_STYLE = {
    **BUTTON_STYLE,
    "background": "#334155",
}

TEXTAREA_STYLE = {
    **INPUT_STYLE,
    "height": "180px",
    "fontFamily": "Menlo, Consolas, monospace",
    "fontSize": "12px",
}


def empty_3d_figure(title):
    fig = go.Figure()
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=45, b=0),
    )
    return fig


def empty_2d_figure(title, x_title="", y_title=""):
    fig = go.Figure()
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        xaxis_title=x_title,
        yaxis_title=y_title,
        margin=dict(l=50, r=20, t=45, b=45),
    )
    return fig


def make_3d_curve_figure(x, y, z, title):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=x,
            y=y,
            z=z,
            mode="lines",
            line=dict(width=5),
            name=title,
        )
    )
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        scene=dict(
            xaxis_title="X",
            yaxis_title="Y",
            zaxis_title="Z",
            aspectmode="data",
        ),
        margin=dict(l=0, r=0, t=45, b=0),
    )
    return fig


def make_2d_curve_figure(x, y, title, x_label, y_label):
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x,
            y=y,
            mode="lines",
            name=title,
        )
    )
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        xaxis_title=x_label,
        yaxis_title=y_label,
        margin=dict(l=50, r=20, t=45, b=45),
    )
    return fig


def make_surface_figure(X, Y, Z, title="Stability"):
    fig = go.Figure()
    fig.add_trace(
        go.Surface(
            x=X,
            y=Y,
            z=Z,
            colorscale="Viridis",
            showscale=True,
        )
    )
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        scene=dict(
            xaxis_title="Offset / Amplitude axis",
            yaxis_title="Sequence axis",
            zaxis_title="Quality / Angle",
            aspectmode="cube",
        ),
        margin=dict(l=0, r=0, t=45, b=0),
    )
    return fig


def make_heatmap_figure(matrix, title="Quality Image"):
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=matrix,
            colorscale="Jet",
            colorbar=dict(title="log10(1-Q)"),
        )
    )
    fig.update_layout(
        title=title,
        template="plotly_dark",
        height=680,
        paper_bgcolor="#111827",
        plot_bgcolor="#111827",
        xaxis_title="Pulse sequence index",
        yaxis_title="Offset / amplitude sample",
        margin=dict(l=50, r=20, t=45, b=45),
    )
    return fig


def build_analyzer_params(
    pulse_sequence,
    pulse_amplitude_khz,
    vector_length,
    time_per_pulse_us,
    inpo_fact,
    x_expand,
    offset_khz,
    scaling_percent,
    calculation_method,
    initial_x,
    initial_y,
    initial_z,
    calculation_language,
):
    return AnalyzerParams(
        pulse_sequence=pulse_sequence,
        pulse_amplitude_khz=float(pulse_amplitude_khz),
        vector_length=float(vector_length),
        time_per_pulse_us=float(time_per_pulse_us),
        inpo_fact=int(inpo_fact),
        x_expand=float(x_expand),
        offset_khz=float(offset_khz),
        scaling_percent=float(scaling_percent),
        calculation_method=int(calculation_method),
        initial_vector=np.array(
            [float(initial_x), float(initial_y), float(initial_z)],
            dtype=float,
        ),
        calculation_language=calculation_language,
    )


def curve_data_string(result):
    return (
        f"Pulse Sequence ({round(1E3 * result.totalRot) / 1E3}°)\n"
        f"Angle to x-Axis= {result.phi}°;\n"
        f"number of pulses = {result.numberOfPulses};\n"
        f"Axy= {result.Axy};\n"
        f"Axz= {result.Axz},\n"
        f"Ayz= {result.Ayz};\n"
        f"Integrated Curvature= {result.integrated_curvature};\n"
        f"Integrated Torsion= {result.integrated_torsion};\n"
        f"Integrated abs Torsion= {result.integrated_absolut_torsion};"
    )


def info_string(update_count, params):
    maximum_amplitude_hz = (
        params.pulse_amplitude_khz * params.scaling_percent * 1000 / 100
    )
    offset_hz = params.offset_khz * 1000
    time_per_pulse_s = params.time_per_pulse_us * 1e-6

    return (
        f"Update Number: {update_count}\n"
        f"Time per Pulse: {time_per_pulse_s * 10**6} µs\n"
        f"Vector Length: {params.vector_length}\n"
        f"Maximum Amplitude: {maximum_amplitude_hz / 1000} kHz\n"
        f"Offset: {offset_hz / 1000} kHz\n"
        f"InpoFact: {params.inpo_fact}\n"
        f"x Expand: {params.x_expand}\n"
        f"Calculation Method: {params.calculation_method}\n"
        f"Initial Vector: {params.initial_vector}\n"
        f"Calculation Language: {params.calculation_language}"
    )


def result_to_store(result):
    return {
        "CM": result.CM.tolist(),
        "VM": result.VM.tolist(),
        "PS_mat": result.PS_mat.tolist(),
        "totalRot": result.totalRot,
        "phi": result.phi,
        "numberOfPulses": result.numberOfPulses,
        "Axy": result.Axy,
        "Axz": result.Axz,
        "Ayz": result.Ayz,
        "arc_length": result.arc_length.tolist(),
        "curvature": result.curvature.tolist(),
        "torsion": result.torsion.tolist(),
        "integrated_curvature": result.integrated_curvature,
        "integrated_torsion": result.integrated_torsion,
        "integrated_absolut_torsion": result.integrated_absolut_torsion,
        "avg_curvature": result.avg_curvature,
        "avg_torsion": result.avg_torsion,
    }


def store_to_arrays(store):
    return {
        "CM": np.array(store["CM"]),
        "VM": np.array(store["VM"]),
        "PS_mat": np.array(store["PS_mat"]),
        "arc_length": np.array(store["arc_length"]),
        "curvature": np.array(store["curvature"]),
        "torsion": np.array(store["torsion"]),
    }


def next_numbered_dir(base_dir, prefix):
    count = 0
    while os.path.exists(os.path.join(base_dir, f"{prefix}{count}")):
        count += 1
    return os.path.join(base_dir, f"{prefix}{count}")


def labeled_input(label, component):
    return html.Div(
        [
            html.Label(label, style={"fontSize": "13px", "color": "#cbd5e1"}),
            component,
        ],
        style={"marginBottom": "10px"},
    )


def build_layout():
    initial_vector = get_default_initial_vector()

    return html.Div(
        style=APP_STYLE,
        children=[
            dcc.Store(id="analysis-store"),
            dcc.Store(id="update-count-store", data=0),

            html.Div(
                [
                    html.H1(
                        "Pulse Sequence Analyzer",
                        style={"margin": "0", "fontSize": "28px"},
                    ),
                    html.Div(
                        "Dash/Plotly interface — Tkinter fallback remains available.",
                        style={"color": "#94a3b8", "marginTop": "5px"},
                    ),
                ],
                style={"marginBottom": "18px"},
            ),

            html.Div(
                [
                    html.Div(
                        style={**CARD_STYLE, "width": "360px", "flexShrink": 0},
                        children=[
                            html.H3("Inputs", style={"marginTop": 0}),

                            labeled_input(
                                "Pulse Sequence Path",
                                dcc.Input(
                                    id="pulse-sequence",
                                    type="text",
                                    value=DEFAULTS.pulse_sequence,
                                    style=INPUT_STYLE,
                                ),
                            ),

                            html.Div(
                                [
                                    html.Button(
                                        "Refresh / Calculate",
                                        id="refresh-button",
                                        n_clicks=0,
                                        style=BUTTON_STYLE,
                                    ),
                                    html.Button(
                                        "Reset Defaults",
                                        id="reset-button",
                                        n_clicks=0,
                                        style={
                                            **SECONDARY_BUTTON_STYLE,
                                            "marginLeft": "8px",
                                        },
                                    ),
                                ],
                                style={"marginBottom": "16px"},
                            ),

                            html.H4("Pulse Parameters"),
                            labeled_input(
                                "Pulse Amplitude [kHz]",
                                dcc.Input(
                                    id="pulse-amplitude",
                                    type="number",
                                    value=10,
                                    style=INPUT_STYLE,
                                ),
                            ),
                            labeled_input(
                                "Vector Length",
                                dcc.Input(
                                    id="vector-length",
                                    type="number",
                                    value=1,
                                    style=INPUT_STYLE,
                                ),
                            ),
                            labeled_input(
                                "Time/pulse [µs]",
                                dcc.Input(
                                    id="time-per-pulse",
                                    type="number",
                                    value=0.5,
                                    style=INPUT_STYLE,
                                ),
                            ),
                            labeled_input(
                                "Inpofact",
                                dcc.Input(
                                    id="inpo-fact",
                                    type="number",
                                    value=5,
                                    step=1,
                                    style=INPUT_STYLE,
                                ),
                            ),
                            labeled_input(
                                "x_Expand",
                                dcc.Input(
                                    id="x-expand",
                                    type="number",
                                    value=0,
                                    style=INPUT_STYLE,
                                ),
                            ),

                            html.H4("Live Controls"),
                            html.Label("Amplitude / Scaling [%]"),
                            dcc.Slider(
                                id="scaling-slider",
                                min=0,
                                max=200,
                                step=1,
                                value=100,
                                marks={0: "0", 100: "100", 200: "200"},
                                updatemode="mouseup",
                            ),
                            dcc.Input(
                                id="scaling",
                                type="number",
                                value=100,
                                style={**INPUT_STYLE, "marginTop": "6px"},
                            ),

                            html.Label(
                                "Offset [kHz]",
                                style={"display": "block", "marginTop": "12px"},
                            ),
                            dcc.Slider(
                                id="offset-slider",
                                min=-20,
                                max=20,
                                step=0.1,
                                value=0,
                                marks={-20: "-20", 0: "0", 20: "20"},
                                updatemode="mouseup",
                            ),
                            dcc.Input(
                                id="offset",
                                type="number",
                                value=0,
                                style={**INPUT_STYLE, "marginTop": "6px"},
                            ),

                            html.H4("Settings"),
                            labeled_input(
                                "Calculation Method",
                                dcc.RadioItems(
                                    id="calculation-method",
                                    options=[
                                        {
                                            "label": "Rotation matrices",
                                            "value": 1,
                                        },
                                        {"label": "Helices", "value": 2},
                                    ],
                                    value=DEFAULTS.calculation_method,
                                    style={"color": "#e5e7eb"},
                                ),
                            ),
                            labeled_input(
                                "Calculation Language",
                                dcc.RadioItems(
                                    id="calculation-language",
                                    options=[
                                        {"label": "Rust", "value": "Rust"},
                                        {"label": "Python", "value": "Python"},
                                    ],
                                    value=DEFAULTS.calculation_language,
                                    style={"color": "#e5e7eb"},
                                ),
                            ),

                            html.H4("Initial Vector"),
                            html.Div(
                                [
                                    dcc.Input(
                                        id="initial-x",
                                        type="number",
                                        value=float(initial_vector[0]),
                                        style={**INPUT_STYLE, "width": "31%"},
                                    ),
                                    dcc.Input(
                                        id="initial-y",
                                        type="number",
                                        value=float(initial_vector[1]),
                                        style={
                                            **INPUT_STYLE,
                                            "width": "31%",
                                            "marginLeft": "3%",
                                        },
                                    ),
                                    dcc.Input(
                                        id="initial-z",
                                        type="number",
                                        value=float(initial_vector[2]),
                                        style={
                                            **INPUT_STYLE,
                                            "width": "31%",
                                            "marginLeft": "3%",
                                        },
                                    ),
                                ]
                            ),
                        ],
                    ),

                    html.Div(
                        style={"flex": 1, "minWidth": 0},
                        children=[
                            dcc.Tabs(
                                id="tabs",
                                value="tab-error",
                                children=[
                                    dcc.Tab(
                                        label="Error Curve",
                                        value="tab-error",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Graph(
                                                    id="error-curve-graph",
                                                    figure=empty_3d_figure(
                                                        "Error Curve"
                                                    ),
                                                ),
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Trajectory",
                                        value="tab-trajectory",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Graph(
                                                    id="trajectory-graph",
                                                    figure=empty_3d_figure(
                                                        "Trajectory"
                                                    ),
                                                ),
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Curvature",
                                        value="tab-curvature",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Graph(
                                                    id="curvature-graph",
                                                    figure=empty_2d_figure(
                                                        "Curvature",
                                                        "arc_length",
                                                        "curvature",
                                                    ),
                                                ),
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Torsion",
                                        value="tab-torsion",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Graph(
                                                    id="torsion-graph",
                                                    figure=empty_2d_figure(
                                                        "Torsion",
                                                        "arc_length",
                                                        "torsion",
                                                    ),
                                                ),
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Stability",
                                        value="tab-stability",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=[
                                                    html.Div(
                                                        [
                                                            dcc.RadioItems(
                                                                id="stability-method",
                                                                options=[
                                                                    {
                                                                        "label": "PP",
                                                                        "value": 1,
                                                                    },
                                                                    {
                                                                        "label": "UR",
                                                                        "value": 2,
                                                                    },
                                                                    {
                                                                        "label": "Angle",
                                                                        "value": 3,
                                                                    },
                                                                ],
                                                                value=1,
                                                                inline=True,
                                                            ),
                                                            dcc.Input(
                                                                id="stability-amplitude-range",
                                                                type="number",
                                                                value=20,
                                                                placeholder="Amplitude range %",
                                                                style={
                                                                    **INPUT_STYLE,
                                                                    "width": "190px",
                                                                    "marginLeft": "16px",
                                                                },
                                                            ),
                                                            dcc.Input(
                                                                id="stability-offset-range",
                                                                type="number",
                                                                value=20,
                                                                placeholder="Offset range kHz",
                                                                style={
                                                                    **INPUT_STYLE,
                                                                    "width": "170px",
                                                                    "marginLeft": "8px",
                                                                },
                                                            ),
                                                            html.Button(
                                                                "Calculate Stability",
                                                                id="stability-button",
                                                                n_clicks=0,
                                                                style={
                                                                    **BUTTON_STYLE,
                                                                    "marginLeft": "8px",
                                                                },
                                                            ),
                                                        ],
                                                        style={
                                                            "display": "flex",
                                                            "alignItems": "center",
                                                            "gap": "8px",
                                                            "marginBottom": "12px",
                                                        },
                                                    ),
                                                    dcc.Loading(
                                                        dcc.Graph(
                                                            id="stability-graph",
                                                            figure=empty_3d_figure(
                                                                "Stability"
                                                            ),
                                                        )
                                                    ),
                                                ],
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Pulse Sequence",
                                        value="tab-ps",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Textarea(
                                                    id="ps-matrix-text",
                                                    value="Pulse sequence will display here",
                                                    readOnly=True,
                                                    style={
                                                        **TEXTAREA_STYLE,
                                                        "height": "650px",
                                                    },
                                                ),
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Options",
                                        value="tab-options",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=[
                                                    html.H3("Export Data"),
                                                    labeled_input(
                                                        "Export Directory",
                                                        dcc.Input(
                                                            id="export-directory",
                                                            type="text",
                                                            value=DEFAULTS.export_directory,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    html.Button(
                                                        "Export Current Data",
                                                        id="export-button",
                                                        n_clicks=0,
                                                        style=BUTTON_STYLE,
                                                    ),
                                                    html.Button(
                                                        "Export Directory Data",
                                                        id="export-dir-button",
                                                        n_clicks=0,
                                                        style={
                                                            **SECONDARY_BUTTON_STYLE,
                                                            "marginLeft": "8px",
                                                        },
                                                    ),
                                                    html.Div(
                                                        id="export-status",
                                                        style={
                                                            "marginTop": "10px",
                                                            "color": "#93c5fd",
                                                        },
                                                    ),

                                                    html.Hr(),

                                                    html.H3("Create Quality Images"),
                                                    html.Div(
                                                        [
                                                            html.Div(
                                                                [
                                                                    html.Label(
                                                                        "Changing Variable"
                                                                    ),
                                                                    dcc.RadioItems(
                                                                        id="quality-changing-variable",
                                                                        options=[
                                                                            {
                                                                                "label": "Amplitude",
                                                                                "value": 1,
                                                                            },
                                                                            {
                                                                                "label": "Offset",
                                                                                "value": 2,
                                                                            },
                                                                        ],
                                                                        value=1,
                                                                    ),
                                                                ],
                                                                style={
                                                                    "width": "48%",
                                                                    "display": "inline-block",
                                                                },
                                                            ),
                                                            html.Div(
                                                                [
                                                                    html.Label(
                                                                        "Calculation Type"
                                                                    ),
                                                                    dcc.RadioItems(
                                                                        id="quality-calc-type",
                                                                        options=[
                                                                            {
                                                                                "label": "SS",
                                                                                "value": 1,
                                                                            },
                                                                            {
                                                                                "label": "UR",
                                                                                "value": 2,
                                                                            },
                                                                        ],
                                                                        value=1,
                                                                    ),
                                                                ],
                                                                style={
                                                                    "width": "48%",
                                                                    "display": "inline-block",
                                                                    "marginLeft": "4%",
                                                                },
                                                            ),
                                                        ]
                                                    ),
                                                    labeled_input(
                                                        "Range [kHz]",
                                                        dcc.Input(
                                                            id="quality-range",
                                                            type="number",
                                                            value=40,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    labeled_input(
                                                        "Maximum Amplitude [kHz]",
                                                        dcc.Input(
                                                            id="quality-umax",
                                                            type="number",
                                                            value=10,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    labeled_input(
                                                        "Resolution [p/kHz]",
                                                        dcc.Input(
                                                            id="quality-resolution",
                                                            type="number",
                                                            value=5,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    labeled_input(
                                                        "Time/pulse [µs]",
                                                        dcc.Input(
                                                            id="quality-time-per-pulse",
                                                            type="number",
                                                            value=0.5,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    labeled_input(
                                                        "Input Directory",
                                                        dcc.Input(
                                                            id="quality-input-directory",
                                                            type="text",
                                                            value=DEFAULTS.quality_input_directory,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    labeled_input(
                                                        "Output Directory",
                                                        dcc.Input(
                                                            id="quality-output-directory",
                                                            type="text",
                                                            value=DEFAULTS.quality_output_directory,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    html.Button(
                                                        "Calculate Quality Image",
                                                        id="quality-button",
                                                        n_clicks=0,
                                                        style=BUTTON_STYLE,
                                                    ),
                                                    html.Div(
                                                        id="quality-status",
                                                        style={
                                                            "marginTop": "10px",
                                                            "color": "#93c5fd",
                                                        },
                                                    ),
                                                    dcc.Loading(
                                                        dcc.Graph(
                                                            id="quality-graph",
                                                            figure=empty_2d_figure(
                                                                "Quality Image"
                                                            ),
                                                        )
                                                    ),
                                                ],
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Edit Pulse Sequence",
                                        value="tab-edit",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=[
                                                    html.H3("Pulse Sequence Editor"),
                                                    labeled_input(
                                                        "Directory Path",
                                                        dcc.Input(
                                                            id="edit-directory",
                                                            type="text",
                                                            value=DEFAULTS.export_directory,
                                                            style=INPUT_STYLE,
                                                        ),
                                                    ),
                                                    dcc.Textarea(
                                                        id="edit-text",
                                                        value=(
                                                            "R(B1 amplitude [%], Phase [°], Duration "
                                                            "(multiple of T sub pulse) [µs]);\n"
                                                            "PR(B1 amplitude [%], starting Phase [°], ending Phase [°], "
                                                            "Duration (multiple of T sub pulse) [µs]); eg:\n"
                                                            "PR(100, 0, 360, 90);\n"
                                                            "R(99.6, 180, 50.5);\n"
                                                            "R(100, 0.1, 60);"
                                                        ),
                                                        style={
                                                            **TEXTAREA_STYLE,
                                                            "height": "420px",
                                                        },
                                                    ),
                                                    html.Button(
                                                        "Save and Use Edited Pulse Sequence",
                                                        id="edit-save-button",
                                                        n_clicks=0,
                                                        style={
                                                            **BUTTON_STYLE,
                                                            "marginTop": "10px",
                                                        },
                                                    ),
                                                    html.Div(
                                                        id="edit-status",
                                                        style={
                                                            "marginTop": "10px",
                                                            "color": "#93c5fd",
                                                        },
                                                    ),
                                                ],
                                            )
                                        ],
                                    ),
                                    dcc.Tab(
                                        label="Documentation",
                                        value="tab-doc",
                                        children=[
                                            html.Div(
                                                style=CARD_STYLE,
                                                children=dcc.Markdown(
                                                    DOCUMENTATION_STR.replace(
                                                        "•", "-"
                                                    ),
                                                    style={
                                                        "whiteSpace": "pre-wrap",
                                                        "lineHeight": "1.45",
                                                    },
                                                ),
                                            )
                                        ],
                                    ),
                                ],
                                colors={
                                    "border": "#334155",
                                    "primary": "#2563eb",
                                    "background": "#111827",
                                },
                            ),

                            html.Div(
                                style={
                                    "display": "grid",
                                    "gridTemplateColumns": "1fr 1fr",
                                    "gap": "12px",
                                    "marginTop": "12px",
                                },
                                children=[
                                    html.Div(
                                        style=CARD_STYLE,
                                        children=[
                                            html.H3(
                                                "Curve Data",
                                                style={"marginTop": 0},
                                            ),
                                            dcc.Textarea(
                                                id="curve-data-text",
                                                value="Curve data will appear here.",
                                                readOnly=True,
                                                style=TEXTAREA_STYLE,
                                            ),
                                        ],
                                    ),
                                    html.Div(
                                        style=CARD_STYLE,
                                        children=[
                                            html.H3(
                                                "Infos and Errors",
                                                style={"marginTop": 0},
                                            ),
                                            dcc.Textarea(
                                                id="info-error-text",
                                                value="Status messages will appear here.",
                                                readOnly=True,
                                                style=TEXTAREA_STYLE,
                                            ),
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),
                ],
                style={
                    "display": "flex",
                    "gap": "14px",
                    "alignItems": "flex-start",
                },
            ),
        ],
    )


def register_callbacks(app):
    @app.callback(
        Output("scaling", "value"),
        Input("scaling-slider", "value"),
        prevent_initial_call=True,
    )
    def sync_scaling_from_slider(value):
        return value

    @app.callback(
        Output("scaling-slider", "value"),
        Input("scaling", "value"),
        prevent_initial_call=True,
    )
    def sync_scaling_to_slider(value):
        if value is None:
            return no_update
        return value

    @app.callback(
        Output("offset", "value"),
        Input("offset-slider", "value"),
        prevent_initial_call=True,
    )
    def sync_offset_from_slider(value):
        return value

    @app.callback(
        Output("offset-slider", "value"),
        Input("offset", "value"),
        prevent_initial_call=True,
    )
    def sync_offset_to_slider(value):
        if value is None:
            return no_update
        return value

    @app.callback(
        Output("pulse-sequence", "value"),
        Output("pulse-amplitude", "value"),
        Output("vector-length", "value"),
        Output("time-per-pulse", "value"),
        Output("inpo-fact", "value"),
        Output("x-expand", "value"),
        Output("offset", "value", allow_duplicate=True),
        Output("scaling", "value", allow_duplicate=True),
        Output("calculation-method", "value"),
        Output("calculation-language", "value"),
        Output("initial-x", "value"),
        Output("initial-y", "value"),
        Output("initial-z", "value"),
        Input("reset-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_defaults(_):
        v = get_default_initial_vector()
        return (
            DEFAULTS.pulse_sequence,
            10,
            DEFAULTS.vector_length,
            DEFAULTS.time_per_pulse_s,
            DEFAULTS.inpo_fact,
            DEFAULTS.x_expand,
            DEFAULTS.offset_hz / 1000,
            DEFAULTS.scaling_percent,
            DEFAULTS.calculation_method,
            DEFAULTS.calculation_language,
            float(v[0]),
            float(v[1]),
            float(v[2]),
        )

    @app.callback(
        Output("analysis-store", "data"),
        Output("update-count-store", "data"),
        Output("error-curve-graph", "figure"),
        Output("trajectory-graph", "figure"),
        Output("curvature-graph", "figure"),
        Output("torsion-graph", "figure"),
        Output("curve-data-text", "value"),
        Output("info-error-text", "value"),
        Output("ps-matrix-text", "value"),
        Input("refresh-button", "n_clicks"),
        Input("offset", "value"),
        Input("scaling", "value"),
        State("update-count-store", "data"),
        State("pulse-sequence", "value"),
        State("pulse-amplitude", "value"),
        State("vector-length", "value"),
        State("time-per-pulse", "value"),
        State("inpo-fact", "value"),
        State("x-expand", "value"),
        State("calculation-method", "value"),
        State("initial-x", "value"),
        State("initial-y", "value"),
        State("initial-z", "value"),
        State("calculation-language", "value"),
    )
    def update_analysis(
        refresh_clicks,
        offset_khz,
        scaling_percent,
        update_count,
        pulse_sequence,
        pulse_amplitude_khz,
        vector_length,
        time_per_pulse_us,
        inpo_fact,
        x_expand,
        calculation_method,
        initial_x,
        initial_y,
        initial_z,
        calculation_language,
    ):
        try:
            update_count = int(update_count or 0) + 1

            params = build_analyzer_params(
                pulse_sequence=pulse_sequence,
                pulse_amplitude_khz=pulse_amplitude_khz,
                vector_length=vector_length,
                time_per_pulse_us=time_per_pulse_us,
                inpo_fact=inpo_fact,
                x_expand=x_expand,
                offset_khz=offset_khz,
                scaling_percent=scaling_percent,
                calculation_method=calculation_method,
                initial_x=initial_x,
                initial_y=initial_y,
                initial_z=initial_z,
                calculation_language=calculation_language,
            )

            result = analyze_pulse_sequence(params)

            CM = result.CM
            VM = result.VM
            expanded_x = CM[:, 0] + np.linspace(
                0,
                params.x_expand,
                np.shape(CM)[0],
            )

            error_fig = make_3d_curve_figure(
                expanded_x,
                CM[:, 1],
                CM[:, 2],
                "Error Curve",
            )
            trajectory_fig = make_3d_curve_figure(
                VM[:, 0],
                VM[:, 1],
                VM[:, 2],
                "Trajectory",
            )
            curvature_fig = make_2d_curve_figure(
                result.arc_length,
                result.curvature,
                "Curvature",
                "arc_length",
                "curvature",
            )
            torsion_fig = make_2d_curve_figure(
                result.arc_length,
                result.torsion,
                "Torsion",
                "arc_length",
                "torsion",
            )

            return (
                result_to_store(result),
                update_count,
                error_fig,
                trajectory_fig,
                curvature_fig,
                torsion_fig,
                curve_data_string(result),
                info_string(update_count, params),
                str(result.PS_mat),
            )

        except Exception:
            err = traceback.format_exc()
            return (
                no_update,
                update_count,
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                err,
                no_update,
            )

    @app.callback(
        Output("stability-graph", "figure"),
        Output("info-error-text", "value", allow_duplicate=True),
        Input("stability-button", "n_clicks"),
        State("analysis-store", "data"),
        State("time-per-pulse", "value"),
        State("vector-length", "value"),
        State("pulse-amplitude", "value"),
        State("scaling", "value"),
        State("stability-amplitude-range", "value"),
        State("stability-offset-range", "value"),
        State("stability-method", "value"),
        State("initial-x", "value"),
        State("initial-y", "value"),
        State("initial-z", "value"),
        State("calculation-language", "value"),
        prevent_initial_call=True,
    )
    def calculate_stability_dash(
        n_clicks,
        store,
        time_per_pulse_us,
        vector_length,
        pulse_amplitude_khz,
        scaling_percent,
        stability_amplitude_range,
        stability_offset_range,
        stability_method,
        initial_x,
        initial_y,
        initial_z,
        calculation_language,
    ):
        if not store:
            return no_update, "Please run Refresh / Calculate first."

        try:
            arrays = store_to_arrays(store)
            PS_mat = arrays["PS_mat"]

            maximum_amplitude_hz = (
                float(pulse_amplitude_khz) * float(scaling_percent) * 1000 / 100
            )

            stability_params = StabilityParams(
                PS_mat=PS_mat,
                time_per_pulse_s=float(time_per_pulse_us) * 1e-6,
                vector_length=float(vector_length),
                maximum_amplitude_hz=maximum_amplitude_hz,
                scaling_range_percent=float(stability_amplitude_range),
                offset_range_khz=float(stability_offset_range),
                stability_calculation_method=int(stability_method),
                initial_vector=np.array(
                    [float(initial_x), float(initial_y), float(initial_z)]
                ),
                calculation_language=calculation_language,
            )

            X, Y, Z, quality, axis, angle = calculate_stability_backend(
                stability_params
            )

            msg = (
                f"Pulse Sequence stability\n"
                f"Offset range: {stability_offset_range} kHz\n"
                f"Amplitude range: {stability_amplitude_range} %\n"
                f"rotation axis: {axis}\n"
                f"rotation angle: {angle} °\n"
                f"average quality: {quality}\n"
            )

            return make_surface_figure(X, Y, Z, "Stability"), msg

        except Exception:
            return no_update, traceback.format_exc()

    @app.callback(
        Output("quality-graph", "figure"),
        Output("quality-status", "children"),
        Output("info-error-text", "value", allow_duplicate=True),
        Input("quality-button", "n_clicks"),
        State("quality-input-directory", "value"),
        State("quality-output-directory", "value"),
        State("quality-range", "value"),
        State("quality-time-per-pulse", "value"),
        State("quality-umax", "value"),
        State("initial-x", "value"),
        State("initial-y", "value"),
        State("initial-z", "value"),
        State("quality-calc-type", "value"),
        State("quality-changing-variable", "value"),
        State("quality-resolution", "value"),
        State("calculation-language", "value"),
        prevent_initial_call=True,
    )
    def calculate_quality_dash(
        n_clicks,
        input_directory,
        output_directory,
        quality_range,
        quality_time_per_pulse,
        quality_umax,
        initial_x,
        initial_y,
        initial_z,
        quality_calc_type,
        quality_changing_variable,
        quality_resolution,
        calculation_language,
    ):
        try:
            progress_messages = []

            def progress_callback(message):
                progress_messages.append(message)

            params = QualityImageParams(
                input_directory=input_directory,
                range_hz=float(quality_range) * 1000,
                time_per_pulse_s=float(quality_time_per_pulse) * 1e-6,
                maximum_amplitude_hz=float(quality_umax) * 1000,
                initial_vector=np.array(
                    [float(initial_x), float(initial_y), float(initial_z)]
                ),
                calculation_type=int(quality_calc_type),
                changing_variable=int(quality_changing_variable),
                resolution=float(quality_resolution),
                calculation_language=calculation_language,
                progress_callback=progress_callback,
            )

            matrix = calculate_quality_image_backend(params)

            output_dir = next_numbered_dir(output_directory, "QualityImage")
            os.mkdir(output_dir)

            np.savetxt(
                os.path.join(output_dir, "QualityMatrix.csv"),
                matrix,
                delimiter=",",
            )
            plt.imsave(
                os.path.join(output_dir, "QualityImage.png"),
                matrix,
                cmap="jet",
            )

            status = f"Quality Image exported to: {output_dir}"
            info = "\n".join(progress_messages[-10:] + [status])

            return make_heatmap_figure(matrix, "Quality Image"), status, info

        except Exception:
            err = traceback.format_exc()
            return no_update, err, err

    @app.callback(
        Output("export-status", "children"),
        Output("info-error-text", "value", allow_duplicate=True),
        Input("export-button", "n_clicks"),
        State("analysis-store", "data"),
        State("export-directory", "value"),
        State("pulse-sequence", "value"),
        State("time-per-pulse", "value"),
        State("vector-length", "value"),
        State("pulse-amplitude", "value"),
        State("scaling", "value"),
        State("offset", "value"),
        State("inpo-fact", "value"),
        State("x-expand", "value"),
        State("calculation-method", "value"),
        State("initial-x", "value"),
        State("initial-y", "value"),
        State("initial-z", "value"),
        State("error-curve-graph", "figure"),
        prevent_initial_call=True,
    )
    def export_current_data(
        n_clicks,
        store,
        export_directory,
        pulse_sequence,
        time_per_pulse_us,
        vector_length,
        pulse_amplitude_khz,
        scaling_percent,
        offset_khz,
        inpo_fact,
        x_expand,
        calculation_method,
        initial_x,
        initial_y,
        initial_z,
        error_curve_figure,
    ):
        if not store:
            return "Please run Refresh / Calculate first.", no_update

        try:
            arrays = store_to_arrays(store)

            output_dir = next_numbered_dir(
                export_directory,
                "PulseSequenceAnalyzer_Data",
            )

            maximum_amplitude_hz = (
                float(pulse_amplitude_khz) * float(scaling_percent) * 1000 / 100
            )
            offset_hz = float(offset_khz) * 1000
            time_per_pulse_s = float(time_per_pulse_us) * 1e-6

            export_string = build_curve_data_export_string(
                pulse_sequence=pulse_sequence,
                totalRot=store["totalRot"],
                time_per_pulse_s=time_per_pulse_s,
                vector_length=float(vector_length),
                maximum_amplitude_hz=maximum_amplitude_hz,
                offset_hz=offset_hz,
                inpo_fact=int(inpo_fact),
                x_expand=float(x_expand),
                calculation_method=int(calculation_method),
                initial_vector=np.array(
                    [float(initial_x), float(initial_y), float(initial_z)]
                ),
                phi=store["phi"],
                numberOfPulses=store["numberOfPulses"],
                Axy=store["Axy"],
                Axz=store["Axz"],
                Ayz=store["Ayz"],
                integrated_curvature=store["integrated_curvature"],
                integrated_torsion=store["integrated_torsion"],
                integrated_absolut_torsion=store[
                    "integrated_absolut_torsion"
                ],
            )

            write_analysis_export(
                output_dir=output_dir,
                CM=arrays["CM"],
                VM=arrays["VM"],
                PS_mat=arrays["PS_mat"],
                arc_length=arrays["arc_length"],
                curvature=arrays["curvature"],
                torsion=arrays["torsion"],
                curve_data_string=export_string,
                error_curve_figure=None,
            )

            if error_curve_figure is not None:
                go.Figure(error_curve_figure).write_html(
                    os.path.join(output_dir, "Error_Curve_plot.html")
                )

            status = f"Data exported to: {output_dir}"
            return status, status

        except Exception:
            err = traceback.format_exc()
            return err, err

    @app.callback(
        Output("export-status", "children", allow_duplicate=True),
        Output("info-error-text", "value", allow_duplicate=True),
        Input("export-dir-button", "n_clicks"),
        State("export-directory", "value"),
        State("pulse-sequence", "value"),
        State("pulse-amplitude", "value"),
        State("vector-length", "value"),
        State("time-per-pulse", "value"),
        State("inpo-fact", "value"),
        State("x-expand", "value"),
        State("offset", "value"),
        State("scaling", "value"),
        State("calculation-method", "value"),
        State("initial-x", "value"),
        State("initial-y", "value"),
        State("initial-z", "value"),
        State("calculation-language", "value"),
        prevent_initial_call=True,
    )
    def export_directory_data(
        n_clicks,
        export_directory,
        active_pulse_sequence,
        pulse_amplitude_khz,
        vector_length,
        time_per_pulse_us,
        inpo_fact,
        x_expand,
        offset_khz,
        scaling_percent,
        calculation_method,
        initial_x,
        initial_y,
        initial_z,
        calculation_language,
    ):
        try:
            base_output_dir = next_numbered_dir(
                export_directory,
                "PulseSequenceAnalyzer_Dir_Data",
            )
            os.mkdir(base_output_dir)

            dirname = os.path.dirname(active_pulse_sequence)
            files = sorted(os.listdir(dirname))
            ps_names = [f for f in files if "bruker" in f]

            for ps_name in ps_names:
                ps_path = os.path.join(dirname, ps_name)

                params = build_analyzer_params(
                    pulse_sequence=ps_path,
                    pulse_amplitude_khz=pulse_amplitude_khz,
                    vector_length=vector_length,
                    time_per_pulse_us=time_per_pulse_us,
                    inpo_fact=inpo_fact,
                    x_expand=x_expand,
                    offset_khz=offset_khz,
                    scaling_percent=scaling_percent,
                    calculation_method=calculation_method,
                    initial_x=initial_x,
                    initial_y=initial_y,
                    initial_z=initial_z,
                    calculation_language=calculation_language,
                )

                result = analyze_pulse_sequence(params)
                output_dir = os.path.join(
                    base_output_dir,
                    ps_name.replace(".bruker", ""),
                )

                export_string = build_curve_data_export_string(
                    pulse_sequence=ps_path,
                    totalRot=result.totalRot,
                    time_per_pulse_s=params.time_per_pulse_us * 1e-6,
                    vector_length=params.vector_length,
                    maximum_amplitude_hz=(
                        params.pulse_amplitude_khz
                        * params.scaling_percent
                        * 1000
                        / 100
                    ),
                    offset_hz=params.offset_khz * 1000,
                    inpo_fact=params.inpo_fact,
                    x_expand=params.x_expand,
                    calculation_method=params.calculation_method,
                    initial_vector=params.initial_vector,
                    phi=result.phi,
                    numberOfPulses=result.numberOfPulses,
                    Axy=result.Axy,
                    Axz=result.Axz,
                    Ayz=result.Ayz,
                    integrated_curvature=result.integrated_curvature,
                    integrated_torsion=result.integrated_torsion,
                    integrated_absolut_torsion=result.integrated_absolut_torsion,
                )

                write_analysis_export(
                    output_dir=output_dir,
                    CM=result.CM,
                    VM=result.VM,
                    PS_mat=result.PS_mat,
                    arc_length=result.arc_length,
                    curvature=result.curvature,
                    torsion=result.torsion,
                    curve_data_string=export_string,
                    error_curve_figure=None,
                )

            status = f"Directory data exported to: {base_output_dir}"
            return status, status

        except Exception:
            err = traceback.format_exc()
            return err, err

    @app.callback(
        Output("edit-status", "children"),
        Output("pulse-sequence", "value", allow_duplicate=True),
        Input("edit-save-button", "n_clicks"),
        State("edit-directory", "value"),
        State("edit-text", "value"),
        State("time-per-pulse", "value"),
        prevent_initial_call=True,
    )
    def save_edited_pulse_sequence(
        n_clicks,
        edit_directory,
        edit_text,
        time_per_pulse_us,
    ):
        try:
            ps_array = parse_pulse_sequence(
                edit_text.strip(),
                subpulse_dt_us=float(time_per_pulse_us),
            )

            os.makedirs(edit_directory, exist_ok=True)
            output_path = os.path.join(edit_directory, "Edited_Pulse_Sequence.csv")
            np.savetxt(output_path, ps_array, delimiter=",")

            return f"Saved edited pulse sequence to: {output_path}", output_path

        except Exception:
            return traceback.format_exc(), no_update


def create_app():
    app = Dash(__name__, suppress_callback_exceptions=True)
    app.title = "Pulse Sequence Analyzer"
    app.layout = build_layout()
    register_callbacks(app)
    return app


def main():
    app = create_app()
    app.run(debug=True, host="127.0.0.1", port=8050)