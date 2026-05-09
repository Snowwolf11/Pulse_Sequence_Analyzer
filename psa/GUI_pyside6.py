#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import traceback
import numpy as np
import matplotlib

matplotlib.use("QtAgg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSlider,
    QSpinBox,
    QDoubleSpinBox,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

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
from psa.documentation import DOCUMENTATION_STR
from psa.app_config import DEFAULTS, get_default_initial_vector
from psa.pulse_string_to_array import parse_pulse_sequence

np.set_printoptions(threshold=np.inf, linewidth=100000)


DARK_BG = "#0f172a"
PANEL_BG = "#111827"
PANEL_2 = "#1e293b"
TEXT_FG = "#e5e7eb"
MUTED_FG = "#94a3b8"
ACCENT = "#2563eb"
ENTRY_BG = "#020617"
LINE_COLOR = "#38bdf8"


def app_stylesheet():
    return f"""
    QWidget {{
        background-color: {DARK_BG};
        color: {TEXT_FG};
        font-size: 13px;
    }}

    QMainWindow {{
        background-color: {DARK_BG};
    }}

    QGroupBox {{
        border: 1px solid #334155;
        border-radius: 10px;
        margin-top: 12px;
        padding: 10px;
        background-color: {PANEL_BG};
        font-weight: bold;
    }}

    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 12px;
        padding: 0 4px;
    }}

    QLineEdit, QTextEdit, QSpinBox, QDoubleSpinBox {{
        background-color: {ENTRY_BG};
        color: {TEXT_FG};
        border: 1px solid #475569;
        border-radius: 6px;
        padding: 5px;
        selection-background-color: {ACCENT};
    }}

    QPushButton {{
        background-color: {PANEL_2};
        color: {TEXT_FG};
        border: 1px solid #475569;
        border-radius: 8px;
        padding: 7px 10px;
        font-weight: 600;
    }}

    QPushButton:hover {{
        background-color: #334155;
    }}

    QPushButton:pressed {{
        background-color: {ACCENT};
    }}

    QTabWidget::pane {{
        border: 1px solid #334155;
        border-radius: 8px;
        background-color: {PANEL_BG};
    }}

    QTabBar::tab {{
        background: {PANEL_2};
        color: {TEXT_FG};
        padding: 9px 14px;
        border-top-left-radius: 8px;
        border-top-right-radius: 8px;
        margin-right: 2px;
    }}

    QTabBar::tab:selected {{
        background: {ACCENT};
        color: white;
    }}

    QRadioButton, QCheckBox {{
        color: {TEXT_FG};
        spacing: 8px;
    }}

    QSlider::groove:horizontal {{
        border: 1px solid #334155;
        height: 6px;
        background: {PANEL_2};
        border-radius: 3px;
    }}

    QSlider::handle:horizontal {{
        background: {ACCENT};
        border: 1px solid #93c5fd;
        width: 16px;
        margin: -6px 0;
        border-radius: 8px;
    }}
    """


class MplPlot(QWidget):
    def __init__(self, projection=None):
        super().__init__()
        self.figure = Figure(figsize=(5, 4))
        self.figure.patch.set_facecolor(PANEL_BG)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        if projection == "3d":
            self.ax = self.figure.add_subplot(111, projection="3d")
        else:
            self.ax = self.figure.add_subplot(111)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

    def draw(self):
        self.canvas.draw_idle()


class SettingsDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent_app = parent
        self.setWindowTitle("Settings")
        self.setMinimumWidth(360)

        layout = QVBoxLayout(self)

        vector_group = QGroupBox("Initial Vector")
        vector_form = QFormLayout(vector_group)

        self.x_spin = QDoubleSpinBox()
        self.y_spin = QDoubleSpinBox()
        self.z_spin = QDoubleSpinBox()
        for spin in (self.x_spin, self.y_spin, self.z_spin):
            spin.setRange(-1000000, 1000000)
            spin.setDecimals(8)

        self.x_spin.setValue(float(parent.initialVector[0]))
        self.y_spin.setValue(float(parent.initialVector[1]))
        self.z_spin.setValue(float(parent.initialVector[2]))

        vector_form.addRow("x:", self.x_spin)
        vector_form.addRow("y:", self.y_spin)
        vector_form.addRow("z:", self.z_spin)
        layout.addWidget(vector_group)

        method_group = QGroupBox("Calculation Method")
        method_layout = QHBoxLayout(method_group)
        self.method_group = QButtonGroup(self)
        self.rotation_radio = QRadioButton("Rotation matrices")
        self.helix_radio = QRadioButton("Helices")
        self.method_group.addButton(self.rotation_radio, 1)
        self.method_group.addButton(self.helix_radio, 2)
        method_layout.addWidget(self.rotation_radio)
        method_layout.addWidget(self.helix_radio)
        if parent.calculationMethod == 1:
            self.rotation_radio.setChecked(True)
        else:
            self.helix_radio.setChecked(True)
        layout.addWidget(method_group)

        language_group = QGroupBox("Calculation Language")
        language_layout = QHBoxLayout(language_group)
        self.language_group = QButtonGroup(self)
        self.rust_radio = QRadioButton("Rust")
        self.python_radio = QRadioButton("Python")
        self.language_group.addButton(self.rust_radio)
        self.language_group.addButton(self.python_radio)
        language_layout.addWidget(self.rust_radio)
        language_layout.addWidget(self.python_radio)
        if parent.calculation_language == "Rust":
            self.rust_radio.setChecked(True)
        else:
            self.python_radio.setChecked(True)
        layout.addWidget(language_group)

        self.display_values_checkbox = QCheckBox("Display values")
        self.display_values_checkbox.setChecked(parent.disp_values)
        layout.addWidget(self.display_values_checkbox)

        buttons = QHBoxLayout()
        save_button = QPushButton("Save")
        close_button = QPushButton("Close")
        reset_button = QPushButton("Reset Analyzer")
        save_button.clicked.connect(self.save)
        close_button.clicked.connect(self.close)
        reset_button.clicked.connect(self.reset)
        buttons.addWidget(save_button)
        buttons.addWidget(reset_button)
        buttons.addWidget(close_button)
        layout.addLayout(buttons)

    def save(self):
        self.parent_app.initialVector = np.array(
            [self.x_spin.value(), self.y_spin.value(), self.z_spin.value()]
        )
        self.parent_app.calculationMethod = self.method_group.checkedId()
        self.parent_app.calculation_language = (
            "Rust" if self.rust_radio.isChecked() else "Python"
        )
        self.parent_app.disp_values = self.display_values_checkbox.isChecked()
        self.accept()

    def reset(self):
        self.parent_app.reset_to_defaults()
        self.accept()


class EditPulseSequenceDialog(QDialog):
    def __init__(self, parent):
        super().__init__(parent)
        self.parent_app = parent
        self.setWindowTitle("Edit Pulse Sequence")
        self.resize(900, 700)

        layout = QVBoxLayout(self)

        self.directory_edit = QLineEdit(DEFAULTS.export_directory)
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_directory)

        dir_layout = QHBoxLayout()
        dir_layout.addWidget(QLabel("Directory Path:"))
        dir_layout.addWidget(self.directory_edit)
        dir_layout.addWidget(browse_button)
        layout.addLayout(dir_layout)

        self.editor = QTextEdit()
        self.editor.setPlainText(
            "R(B1 amplitude [%], Phase [°], Duration (multiple of T sub pulse) [µs]);\n"
            "PR(B1 amplitude [%], starting Phase [°], ending Phase [°], "
            "Duration (multiple of T sub pulse) [µs]); eg:\n"
            "PR(100, 0, 360, 90);\n"
            "R(99.6, 180, 50.5);\n"
            "R(100, 0.1, 60);"
        )
        layout.addWidget(self.editor)

        button_layout = QHBoxLayout()
        save_button = QPushButton("Save and Plot")
        close_button = QPushButton("Close")
        save_button.clicked.connect(self.save_and_plot)
        close_button.clicked.connect(self.close)
        button_layout.addStretch()
        button_layout.addWidget(save_button)
        button_layout.addWidget(close_button)
        layout.addLayout(button_layout)

    def browse_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Directory")
        if directory:
            self.directory_edit.setText(directory)

    def save_and_plot(self):
        try:
            directory = self.directory_edit.text()
            os.makedirs(directory, exist_ok=True)
            output_path = os.path.join(directory, "Edited_Pulse_Sequence.csv")

            ps_array = parse_pulse_sequence(
                self.editor.toPlainText().strip(),
                subpulse_dt_us=self.parent_app.T * 1e6,
            )
            np.savetxt(output_path, ps_array, delimiter=",")

            self.parent_app.pulse_sequence_edit.setText(output_path)
            self.parent_app.update_analysis()
            self.parent_app.append_info(f"Edited pulse sequence saved to: {output_path}")

        except Exception:
            QMessageBox.critical(self, "Error", traceback.format_exc())


class PulseSequenceAnalyzerQt(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Pulse Sequence Analyzer — PySide6")
        self.resize(1500, 950)

        self.PS = DEFAULTS.pulse_sequence
        self.maximumAmplitude = DEFAULTS.maximum_amplitude_hz
        self.Vector_Length = DEFAULTS.vector_length
        self.T = DEFAULTS.time_per_pulse_s
        self.InpoFact = DEFAULTS.inpo_fact
        self.x_Expand = DEFAULTS.x_expand
        self.Offset = DEFAULTS.offset_hz
        self.Scaling = DEFAULTS.scaling_percent
        self.calculationMethod = DEFAULTS.calculation_method
        self.initialVector = get_default_initial_vector()
        self.calculation_language = DEFAULTS.calculation_language

        self.disp_values = False
        self.update_count = 0
        self.export_count = 0
        self.export_dir_count = 0
        self.latest_result = None

        self.play_amplitude_direction = 1
        self.play_offset_direction = 1

        self.play_amplitude_timer = QTimer(self)
        self.play_amplitude_timer.timeout.connect(self.step_play_amplitude)

        self.play_offset_timer = QTimer(self)
        self.play_offset_timer.timeout.connect(self.step_play_offset)   

        self._build_ui()

    def keyPressEvent(self, event):
        if event.key() in (Qt.Key_Return, Qt.Key_Enter):
            self.update_analysis()
            event.accept()
            return

        super().keyPressEvent(event)
        
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QHBoxLayout(central)

        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFixedWidth(420)
        left_panel = QWidget()
        left_scroll.setWidget(left_panel)
        left_layout = QVBoxLayout(left_panel)

        top_buttons = QHBoxLayout()
        menu_button = QPushButton("Menu")
        doc_button = QPushButton("Doc")
        refresh_button = QPushButton("Refresh")
        menu_button.clicked.connect(self.open_settings)
        doc_button.clicked.connect(self.show_documentation)
        refresh_button.clicked.connect(self.update_analysis)
        top_buttons.addWidget(menu_button)
        top_buttons.addWidget(doc_button)
        top_buttons.addWidget(refresh_button)
        left_layout.addLayout(top_buttons)

        ps_group = QGroupBox("Pulse Sequence")
        ps_layout = QGridLayout(ps_group)
        self.pulse_sequence_edit = QLineEdit(self.PS)
        browse_ps_button = QPushButton("Browse")
        edit_ps_button = QPushButton("Edit")
        browse_ps_button.clicked.connect(self.browse_pulse_sequence)
        edit_ps_button.clicked.connect(self.open_editor)

        ps_layout.addWidget(self.pulse_sequence_edit, 0, 0, 1, 3)
        ps_layout.addWidget(browse_ps_button, 1, 1)
        ps_layout.addWidget(edit_ps_button, 1, 2)
        left_layout.addWidget(ps_group)

        param_group = QGroupBox("Parameters")
        form = QFormLayout(param_group)

        self.pulse_amplitude_spin = QDoubleSpinBox()
        self.pulse_amplitude_spin.setRange(-1e9, 1e9)
        self.pulse_amplitude_spin.setDecimals(6)
        self.pulse_amplitude_spin.setValue(10)

        self.vector_length_spin = QDoubleSpinBox()
        self.vector_length_spin.setRange(-1e9, 1e9)
        self.vector_length_spin.setDecimals(6)
        self.vector_length_spin.setValue(1)

        self.time_per_pulse_spin = QDoubleSpinBox()
        self.time_per_pulse_spin.setRange(0, 1e9)
        self.time_per_pulse_spin.setDecimals(6)
        self.time_per_pulse_spin.setValue(0.5)

        self.inpo_fact_spin = QSpinBox()
        self.inpo_fact_spin.setRange(1, 1000000)
        self.inpo_fact_spin.setValue(5)

        self.x_expand_spin = QDoubleSpinBox()
        self.x_expand_spin.setRange(-1e9, 1e9)
        self.x_expand_spin.setDecimals(6)
        self.x_expand_spin.setValue(0)

        self.offset_spin = QDoubleSpinBox()
        self.offset_spin.setRange(-1e9, 1e9)
        self.offset_spin.setDecimals(6)
        self.offset_spin.setValue(0)

        self.scaling_spin = QDoubleSpinBox()
        self.scaling_spin.setRange(0, 1e9)
        self.scaling_spin.setDecimals(6)
        self.scaling_spin.setValue(100)

        form.addRow("Pulse Amplitude [kHz]:", self.pulse_amplitude_spin)
        form.addRow("Vector Length:", self.vector_length_spin)
        form.addRow("Time/pulse [µs]:", self.time_per_pulse_spin)
        form.addRow("Inpofact:", self.inpo_fact_spin)
        form.addRow("x_Expand:", self.x_expand_spin)
        form.addRow("Offset [kHz]:", self.offset_spin)
        form.addRow("Scaling [%]:", self.scaling_spin)

        left_layout.addWidget(param_group)

        slider_group = QGroupBox("Live Controls")
        slider_layout = QVBoxLayout(slider_group)

        self.amplitude_slider = QSlider(Qt.Horizontal)
        self.amplitude_slider.setRange(0, 200)
        self.amplitude_slider.setValue(100)
        self.amplitude_slider.sliderReleased.connect(self.apply_amplitude_slider)

        self.offset_slider = QSlider(Qt.Horizontal)
        self.offset_slider.setRange(-200, 200)
        self.offset_slider.setValue(0)
        self.offset_slider.sliderReleased.connect(self.apply_offset_slider)

        self.play_amp_button = QPushButton("Play Amplitude")
        self.play_off_button = QPushButton("Play Offset")

        self.play_amp_button.clicked.connect(self.toggle_play_amplitude)
        self.play_off_button.clicked.connect(self.toggle_play_offset)

        slider_layout.addWidget(QLabel("Amplitude [%]"))
        slider_layout.addWidget(self.amplitude_slider)
        slider_layout.addWidget(self.play_amp_button)
        slider_layout.addWidget(QLabel("Offset [kHz]"))
        slider_layout.addWidget(self.offset_slider)
        slider_layout.addWidget(self.play_off_button)

        left_layout.addWidget(slider_group)

        self.curve_data_text = QTextEdit()
        self.curve_data_text.setReadOnly(True)
        self.curve_data_text.setMinimumHeight(190)
        curve_group = QGroupBox("Curve Data")
        curve_layout = QVBoxLayout(curve_group)
        curve_layout.addWidget(self.curve_data_text)
        left_layout.addWidget(curve_group)

        self.info_error_text = QTextEdit()
        self.info_error_text.setReadOnly(True)
        self.info_error_text.setMinimumHeight(190)
        info_group = QGroupBox("Infos and Errors")
        info_layout = QVBoxLayout(info_group)
        info_layout.addWidget(self.info_error_text)
        left_layout.addWidget(info_group)

        left_layout.addStretch()

        root.addWidget(left_scroll)

        self.tabs = QTabWidget()
        root.addWidget(self.tabs, stretch=1)

        self.error_plot = MplPlot(projection="3d")
        self.trajectory_plot = MplPlot(projection="3d")
        self.curvature_plot = MplPlot()
        self.torsion_plot = MplPlot()

        self.tabs.addTab(self.error_plot, "Error-Curve")
        self.tabs.addTab(self.trajectory_plot, "Trajectory")
        self.tabs.addTab(self.curvature_plot, "Curvature")
        self.tabs.addTab(self.torsion_plot, "Torsion")

        self._build_stability_tab()
        self._build_pulse_sequence_tab()
        self._build_options_tab()

    def _build_stability_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        controls = QHBoxLayout()
        self.stability_group = QButtonGroup(self)
        self.stability_pp = QRadioButton("PP")
        self.stability_ur = QRadioButton("UR")
        self.stability_angle = QRadioButton("Angle")
        self.stability_group.addButton(self.stability_pp, 1)
        self.stability_group.addButton(self.stability_ur, 2)
        self.stability_group.addButton(self.stability_angle, 3)
        self.stability_pp.setChecked(True)

        self.stability_amplitude_range = QDoubleSpinBox()
        self.stability_amplitude_range.setRange(0, 1e9)
        self.stability_amplitude_range.setValue(20)

        self.stability_offset_range = QDoubleSpinBox()
        self.stability_offset_range.setRange(0, 1e9)
        self.stability_offset_range.setValue(20)

        calc_button = QPushButton("Calculate")
        calc_button.clicked.connect(self.calculate_stability)

        controls.addWidget(self.stability_pp)
        controls.addWidget(self.stability_ur)
        controls.addWidget(self.stability_angle)
        controls.addWidget(QLabel("Amplitude range [%]:"))
        controls.addWidget(self.stability_amplitude_range)
        controls.addWidget(QLabel("Offset range [kHz]:"))
        controls.addWidget(self.stability_offset_range)
        controls.addWidget(calc_button)

        self.stability_plot = MplPlot(projection="3d")

        layout.addLayout(controls)
        layout.addWidget(self.stability_plot)

        self.tabs.addTab(tab, "Stability")

    def _build_pulse_sequence_tab(self):
        self.ps_matrix_text = QTextEdit()
        self.ps_matrix_text.setReadOnly(True)
        self.ps_matrix_text.setPlainText("Pulse sequence will display here.")
        self.tabs.addTab(self.ps_matrix_text, "Pulse_Sequence")

    def _build_options_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        export_group = QGroupBox("Export Data")
        export_layout = QGridLayout(export_group)

        self.export_directory_edit = QLineEdit(DEFAULTS.export_directory)
        export_browse = QPushButton("Browse")
        export_button = QPushButton("Export Data")
        export_dir_button = QPushButton("Export Dir Data")

        export_browse.clicked.connect(self.browse_export_directory)
        export_button.clicked.connect(self.export_current_data)
        export_dir_button.clicked.connect(self.export_directory_data)

        export_layout.addWidget(QLabel("Path:"), 0, 0)
        export_layout.addWidget(self.export_directory_edit, 0, 1)
        export_layout.addWidget(export_browse, 0, 2)
        export_layout.addWidget(export_button, 0, 3)
        export_layout.addWidget(export_dir_button, 0, 4)

        layout.addWidget(export_group)

        quality_group = QGroupBox("Create Quality Images")
        q_layout = QGridLayout(quality_group)

        self.quality_variable_group = QButtonGroup(self)
        self.quality_amp_radio = QRadioButton("Amplitude")
        self.quality_offset_radio = QRadioButton("Offset")
        self.quality_variable_group.addButton(self.quality_amp_radio, 1)
        self.quality_variable_group.addButton(self.quality_offset_radio, 2)
        self.quality_amp_radio.setChecked(True)

        self.quality_type_group = QButtonGroup(self)
        self.quality_ss_radio = QRadioButton("SS")
        self.quality_ur_radio = QRadioButton("UR")
        self.quality_type_group.addButton(self.quality_ss_radio, 1)
        self.quality_type_group.addButton(self.quality_ur_radio, 2)
        self.quality_ss_radio.setChecked(True)

        self.quality_range = QDoubleSpinBox()
        self.quality_range.setRange(0, 1e9)
        self.quality_range.setValue(40)

        self.quality_umax = QDoubleSpinBox()
        self.quality_umax.setRange(0, 1e9)
        self.quality_umax.setValue(10)

        self.quality_resolution = QDoubleSpinBox()
        self.quality_resolution.setRange(0, 1e9)
        self.quality_resolution.setValue(5)

        self.quality_time_per_pulse = QDoubleSpinBox()
        self.quality_time_per_pulse.setRange(0, 1e9)
        self.quality_time_per_pulse.setDecimals(6)
        self.quality_time_per_pulse.setValue(0.5)

        self.quality_input_dir = QLineEdit(DEFAULTS.quality_input_directory)
        self.quality_output_dir = QLineEdit(DEFAULTS.quality_output_directory)

        input_browse = QPushButton("Browse")
        output_browse = QPushButton("Browse")
        calc_quality = QPushButton("Calculate")
        input_browse.clicked.connect(self.browse_quality_input)
        output_browse.clicked.connect(self.browse_quality_output)
        calc_quality.clicked.connect(self.calculate_quality_images)

        row = 0
        q_layout.addWidget(QLabel("Changing Variable:"), row, 0)
        q_layout.addWidget(self.quality_amp_radio, row, 1)
        q_layout.addWidget(self.quality_offset_radio, row, 2)
        q_layout.addWidget(QLabel("Calculation Type:"), row, 3)
        q_layout.addWidget(self.quality_ss_radio, row, 4)
        q_layout.addWidget(self.quality_ur_radio, row, 5)

        row += 1
        q_layout.addWidget(QLabel("Range [kHz]:"), row, 0)
        q_layout.addWidget(self.quality_range, row, 1)
        q_layout.addWidget(QLabel("Maximum Amplitude [kHz]:"), row, 2)
        q_layout.addWidget(self.quality_umax, row, 3)
        q_layout.addWidget(QLabel("Resolution [p/kHz]:"), row, 4)
        q_layout.addWidget(self.quality_resolution, row, 5)

        row += 1
        q_layout.addWidget(QLabel("Time/pulse [µs]:"), row, 0)
        q_layout.addWidget(self.quality_time_per_pulse, row, 1)

        row += 1
        q_layout.addWidget(QLabel("Input Directory:"), row, 0)
        q_layout.addWidget(self.quality_input_dir, row, 1, 1, 4)
        q_layout.addWidget(input_browse, row, 5)

        row += 1
        q_layout.addWidget(QLabel("Output Directory:"), row, 0)
        q_layout.addWidget(self.quality_output_dir, row, 1, 1, 4)
        q_layout.addWidget(output_browse, row, 5)

        row += 1
        q_layout.addWidget(calc_quality, row, 5)

        layout.addWidget(quality_group)

        self.quality_plot = MplPlot()
        layout.addWidget(self.quality_plot)

        self.tabs.addTab(tab, "options")

    def browse_pulse_sequence(self):
        path, _ = QFileDialog.getOpenFileName(self, "Choose Pulse Sequence")
        if path:
            self.pulse_sequence_edit.setText(path)

    def browse_export_directory(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Export Directory")
        if directory:
            self.export_directory_edit.setText(directory)

    def browse_quality_input(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Input Directory")
        if directory:
            self.quality_input_dir.setText(directory)

    def browse_quality_output(self):
        directory = QFileDialog.getExistingDirectory(self, "Choose Output Directory")
        if directory:
            self.quality_output_dir.setText(directory)

    def open_settings(self):
        SettingsDialog(self).exec()

    def open_editor(self):
        EditPulseSequenceDialog(self).exec()

    def show_documentation(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Documentation")
        dialog.resize(950, 700)
        layout = QVBoxLayout(dialog)
        text = QTextEdit()
        text.setReadOnly(True)
        text.setPlainText(DOCUMENTATION_STR)
        layout.addWidget(text)
        close = QPushButton("Close")
        close.clicked.connect(dialog.close)
        layout.addWidget(close)
        dialog.exec()

    def toggle_play_amplitude(self):
        if self.play_amplitude_timer.isActive():
            self.play_amplitude_timer.stop()
            self.play_amp_button.setText("Play Amplitude")
            return

        self.play_amp_button.setText("Stop Amplitude")
        self.play_amplitude_timer.start(250)
    

    def toggle_play_offset(self):
        if self.play_offset_timer.isActive():
            self.play_offset_timer.stop()
            self.play_off_button.setText("Play Offset")
            return

        self.play_off_button.setText("Stop Offset")
        self.play_offset_timer.start(250)


    def step_play_amplitude(self):
        step = 2
        value = self.amplitude_slider.value()
        new_value = value + self.play_amplitude_direction * step

        if new_value >= self.amplitude_slider.maximum():
            new_value = self.amplitude_slider.maximum()
            self.play_amplitude_direction = -1

        elif new_value <= self.amplitude_slider.minimum():
            new_value = self.amplitude_slider.minimum()
            self.play_amplitude_direction = 1

        self.amplitude_slider.setValue(new_value)
        self.scaling_spin.setValue(float(new_value))
        self.update_analysis()


    def step_play_offset(self):
        step = 2
        value = self.offset_slider.value()
        new_value = value + self.play_offset_direction * step

        if new_value >= self.offset_slider.maximum():
            new_value = self.offset_slider.maximum()
            self.play_offset_direction = -1

        elif new_value <= self.offset_slider.minimum():
            new_value = self.offset_slider.minimum()
            self.play_offset_direction = 1

        self.offset_slider.setValue(new_value)
        self.offset_spin.setValue(float(new_value) / 10)
        self.update_analysis()

    def reset_to_defaults(self):
        self.PS = DEFAULTS.pulse_sequence
        self.maximumAmplitude = DEFAULTS.maximum_amplitude_hz
        self.Vector_Length = DEFAULTS.vector_length
        self.T = DEFAULTS.time_per_pulse_s
        self.InpoFact = DEFAULTS.inpo_fact
        self.x_Expand = DEFAULTS.x_expand
        self.Offset = DEFAULTS.offset_hz
        self.Scaling = DEFAULTS.scaling_percent
        self.calculationMethod = DEFAULTS.calculation_method
        self.initialVector = get_default_initial_vector()
        self.calculation_language = DEFAULTS.calculation_language

        self.pulse_sequence_edit.setText(self.PS)
        self.pulse_amplitude_spin.setValue(10)
        self.vector_length_spin.setValue(self.Vector_Length)
        self.time_per_pulse_spin.setValue(self.T)
        self.inpo_fact_spin.setValue(self.InpoFact)
        self.x_expand_spin.setValue(self.x_Expand)
        self.offset_spin.setValue(self.Offset / 1000)
        self.scaling_spin.setValue(self.Scaling)

    def apply_amplitude_slider(self):
        self.scaling_spin.setValue(float(self.amplitude_slider.value()))
        self.update_analysis()

    def apply_offset_slider(self):
        self.offset_spin.setValue(float(self.offset_slider.value()) / 10)
        self.update_analysis()

    def get_params(self):
        return AnalyzerParams(
            pulse_sequence=self.pulse_sequence_edit.text(),
            pulse_amplitude_khz=self.pulse_amplitude_spin.value(),
            vector_length=self.vector_length_spin.value(),
            time_per_pulse_us=self.time_per_pulse_spin.value(),
            inpo_fact=self.inpo_fact_spin.value(),
            x_expand=self.x_expand_spin.value(),
            offset_khz=self.offset_spin.value(),
            scaling_percent=self.scaling_spin.value(),
            calculation_method=self.calculationMethod,
            initial_vector=self.initialVector,
            calculation_language=self.calculation_language,
        )

    def update_values_from_params(self, params):
        self.maximumAmplitude = (
            params.pulse_amplitude_khz * params.scaling_percent * 1000 / 100
        )
        self.Offset = params.offset_khz * 1000
        self.T = params.time_per_pulse_us * 1e-6
        self.PS = params.pulse_sequence
        self.Vector_Length = params.vector_length
        self.InpoFact = params.inpo_fact
        self.x_Expand = params.x_expand
        self.Scaling = params.scaling_percent

        self.amplitude_slider.setValue(int(round(self.Scaling)))
        self.offset_slider.setValue(int(round(self.Offset / 100)))

    def update_analysis(self):
        try:
            self.update_count += 1
            params = self.get_params()
            self.update_values_from_params(params)

            result = analyze_pulse_sequence(params)
            self.latest_result = result

            CM = result.CM
            VM = result.VM
            expanded_x = CM[:, 0] + np.linspace(0, self.x_Expand, np.shape(CM)[0])

            self.plot_3d(
                self.error_plot,
                expanded_x,
                CM[:, 1],
                CM[:, 2],
                "Error-Kurve",
            )
            self.plot_3d(
                self.trajectory_plot,
                VM[:, 0],
                VM[:, 1],
                VM[:, 2],
                "Trajectory",
            )
            self.plot_2d(
                self.curvature_plot,
                result.arc_length,
                result.curvature,
                "Curvature",
                "arc_length",
                "curvature",
            )
            self.plot_2d(
                self.torsion_plot,
                result.arc_length,
                result.torsion,
                "Torsion",
                "arc_length",
                "torsion",
            )

            self.curve_data_text.setPlainText(self.curve_data_string(result))
            self.info_error_text.setPlainText(self.info_error_string(params))
            self.ps_matrix_text.setPlainText(str(result.PS_mat))

            return result

        except Exception:
            self.info_error_text.setPlainText(traceback.format_exc())
            return None

    def scientific_formatter(self, value, pos):
        if np.isclose(value, 0.0):
            return "0"
        magnitude = int(np.floor(np.log10(abs(value))))
        formatted_value = value / (10 ** magnitude)
        if magnitude >= 0:
            return f"{formatted_value:.1f}e{magnitude}"
        return f"{formatted_value:.1f}e-{abs(magnitude)}"

    def style_axes_3d(self, ax):
        ax.set_facecolor(PANEL_BG)
        ax.figure.patch.set_facecolor(PANEL_BG)
        ax.xaxis.pane.set_facecolor((0.07, 0.09, 0.13, 1.0))
        ax.yaxis.pane.set_facecolor((0.07, 0.09, 0.13, 1.0))
        ax.zaxis.pane.set_facecolor((0.07, 0.09, 0.13, 1.0))
        ax.tick_params(colors=TEXT_FG)
        ax.xaxis.label.set_color(TEXT_FG)
        ax.yaxis.label.set_color(TEXT_FG)
        ax.zaxis.label.set_color(TEXT_FG)
        ax.title.set_color(TEXT_FG)

    def style_axes_2d(self, ax):
        ax.set_facecolor(PANEL_BG)
        ax.figure.patch.set_facecolor(PANEL_BG)
        ax.tick_params(colors=TEXT_FG)
        ax.xaxis.label.set_color(TEXT_FG)
        ax.yaxis.label.set_color(TEXT_FG)
        ax.title.set_color(TEXT_FG)
        ax.grid(True, alpha=0.25)

    def plot_3d(self, plot, x, y, z, title):
        ax = plot.ax
        ax.clear()
        self.style_axes_3d(ax)
        ax.plot3D(x, y, z, color=LINE_COLOR, linewidth=1.5)
        ax.set_title(title)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.zaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter))
        ax.zaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter))

        try:
            ax.set_box_aspect(
                (
                    abs(max(x) - min(x)) or 1,
                    abs(max(y) - min(y)) or 1,
                    abs(max(z) - min(z)) or 1,
                )
            )
        except Exception:
            pass

        plot.draw()

    def plot_2d(self, plot, x, y, title, x_label, y_label):
        ax = plot.ax
        ax.clear()
        self.style_axes_2d(ax)
        ax.plot(x, y, color=LINE_COLOR, linewidth=1.5)
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        plot.draw()

    def curve_data_string(self, result):
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

    def info_error_string(self, params):
        return (
            f"Update Number: {self.update_count}\n"
            f"Time per Pulse: {self.T * 10**6} µs\n"
            f"Vector Length: {self.Vector_Length}\n"
            f"Maximum Amplitude: {self.maximumAmplitude / 1000} kHz\n"
            f"Offset: {self.Offset / 1000} kHz\n"
            f"InpoFact: {self.InpoFact}\n"
            f"x Expand: {self.x_Expand}\n"
            f"Calculation Method: {self.calculationMethod}\n"
            f"Initial Vector: {self.initialVector}\n"
            f"Calculation Language: {self.calculation_language}"
        )

    def append_info(self, message):
        current = self.info_error_text.toPlainText()
        self.info_error_text.setPlainText(f"{message}\n\n{current}")

    def calculate_stability(self):
        try:
            self.inpo_fact_spin.setValue(1)
            result = self.update_analysis()
            if result is None:
                return

            stability_params = StabilityParams(
                PS_mat=result.PS_mat,
                time_per_pulse_s=self.T,
                vector_length=self.Vector_Length,
                maximum_amplitude_hz=self.maximumAmplitude,
                scaling_range_percent=self.stability_amplitude_range.value(),
                offset_range_kHz=self.stability_offset_range.value(),
                stability_calculation_method=self.stability_group.checkedId(),
                initial_vector=self.initialVector,
                calculation_language=self.calculation_language,
            )

        except TypeError:
            stability_params = StabilityParams(
                PS_mat=self.latest_result.PS_mat,
                time_per_pulse_s=self.T,
                vector_length=self.Vector_Length,
                maximum_amplitude_hz=self.maximumAmplitude,
                scaling_range_percent=self.stability_amplitude_range.value(),
                offset_range_khz=self.stability_offset_range.value(),
                stability_calculation_method=self.stability_group.checkedId(),
                initial_vector=self.initialVector,
                calculation_language=self.calculation_language,
            )

        try:
            X, Y, Z, quality, axis, angle = calculate_stability_backend(
                stability_params
            )

            ax = self.stability_plot.ax
            ax.clear()
            self.style_axes_3d(ax)
            ax.plot_surface(X, Y, Z, cmap="viridis", edgecolor="black", linewidth=0.1)
            ax.set_title("Stability")
            self.stability_plot.draw()

            msg = (
                f"Pulse Sequence stability\n"
                f"Offset range: {self.stability_offset_range.value()} kHz\n"
                f"Amplitude range: {self.stability_amplitude_range.value()} %\n"
                f"rotation axis: {axis}\n"
                f"rotation angle: {angle} °\n"
                f"average quality: {quality}\n"
            )
            self.append_info(msg)

        except Exception:
            self.info_error_text.setPlainText(traceback.format_exc())

    def calculate_quality_images(self):
        try:
            messages = []

            def progress_callback(message):
                messages.append(message)
                self.append_info(message)
                QApplication.processEvents()

            params = QualityImageParams(
                input_directory=self.quality_input_dir.text(),
                range_hz=self.quality_range.value() * 1000,
                time_per_pulse_s=self.quality_time_per_pulse.value() * 1e-6,
                maximum_amplitude_hz=self.quality_umax.value() * 1000,
                initial_vector=self.initialVector,
                calculation_type=self.quality_type_group.checkedId(),
                changing_variable=self.quality_variable_group.checkedId(),
                resolution=self.quality_resolution.value(),
                calculation_language=self.calculation_language,
                progress_callback=progress_callback,
            )

            matrix = calculate_quality_image_backend(params)

            output_dir = self.next_numbered_dir(
                self.quality_output_dir.text(),
                "QualityImage",
            )
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

            self.quality_plot.figure.clear()
            self.quality_plot.ax = self.quality_plot.figure.add_subplot(111)

            ax = self.quality_plot.ax
            self.style_axes_2d(ax)

            im = ax.imshow(matrix, cmap="jet", aspect="auto")
            ax.set_title("Quality Image")
            ax.set_xlabel("Pulse sequence index")
            ax.set_ylabel("Offset / amplitude sample")

            cbar = self.quality_plot.figure.colorbar(im, ax=ax)
            cbar.ax.tick_params(colors=TEXT_FG)
            cbar.outline.set_edgecolor(TEXT_FG)

            self.quality_plot.draw()

            self.append_info(f"Quality Image exported to: {output_dir}")

        except Exception:
            self.info_error_text.setPlainText(traceback.format_exc())

    def export_current_data(self):
        try:
            result = self.latest_result or self.update_analysis()
            if result is None:
                return

            output_dir = self.next_numbered_dir(
                self.export_directory_edit.text(),
                "PulseSequenceAnalyzer_Data",
            )

            export_string = build_curve_data_export_string(
                pulse_sequence=self.PS,
                totalRot=result.totalRot,
                time_per_pulse_s=self.T,
                vector_length=self.Vector_Length,
                maximum_amplitude_hz=self.maximumAmplitude,
                offset_hz=self.Offset,
                inpo_fact=self.InpoFact,
                x_expand=self.x_Expand,
                calculation_method=self.calculationMethod,
                initial_vector=self.initialVector,
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
                error_curve_figure=self.error_plot.figure,
            )

            self.append_info(f"Data exported to: {output_dir}")

        except Exception:
            self.info_error_text.setPlainText(traceback.format_exc())

    def export_directory_data(self):
        try:
            base_output_dir = self.next_numbered_dir(
                self.export_directory_edit.text(),
                "PulseSequenceAnalyzer_Dir_Data",
            )
            os.mkdir(base_output_dir)

            dirname = os.path.dirname(self.pulse_sequence_edit.text())
            files = sorted(os.listdir(dirname))
            ps_names = [f for f in files if "bruker" in f]

            original_path = self.pulse_sequence_edit.text()

            for ps_name in ps_names:
                ps_path = os.path.join(dirname, ps_name)
                self.pulse_sequence_edit.setText(ps_path)
                result = self.update_analysis()
                if result is None:
                    continue

                output_dir = os.path.join(base_output_dir, ps_name.replace(".bruker", ""))

                export_string = build_curve_data_export_string(
                    pulse_sequence=ps_path,
                    totalRot=result.totalRot,
                    time_per_pulse_s=self.T,
                    vector_length=self.Vector_Length,
                    maximum_amplitude_hz=self.maximumAmplitude,
                    offset_hz=self.Offset,
                    inpo_fact=self.InpoFact,
                    x_expand=self.x_Expand,
                    calculation_method=self.calculationMethod,
                    initial_vector=self.initialVector,
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
                    error_curve_figure=self.error_plot.figure,
                )

                QApplication.processEvents()

            self.pulse_sequence_edit.setText(original_path)
            self.append_info(f"Directory data exported to: {base_output_dir}")

        except Exception:
            self.info_error_text.setPlainText(traceback.format_exc())

    def next_numbered_dir(self, base_dir, prefix):
        count = 0
        while os.path.exists(os.path.join(base_dir, f"{prefix}{count}")):
            count += 1
        return os.path.join(base_dir, f"{prefix}{count}")


def main():
    app = QApplication.instance() or QApplication(sys.argv)
    app.setStyleSheet(app_stylesheet())
    window = PulseSequenceAnalyzerQt()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()