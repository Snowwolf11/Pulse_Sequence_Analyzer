#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 14:25:41 2024

@author: leon
"""

#most important:
######TODO: Jemand schicken und funktioniert einfach!!
######TODO: Error handling
######TODO: plots verbessern: Achsenbeschriftung, datapoints, zoom,...!!!
######TODO: updates bei quality Images kommen erst nachdem komplett fertig gerechnet, threading for GUI?
######TODO: play buttons
#not so immediate
######TODO: Check all calculated values and plots
######TODO: Documentation 
######TODO: Profiling
######TODO: Plane projections, visibility and colour of curve
######TODO: plots exportieren
######TODO: UR stability plot only for 0° & 360° Pulse !!!


from psa.createCurve import *
from psa.calculate_curveData import *
from psa.calculate_curve_stability import *
from psa.calculate_pulse_sequence_quality_images import *

import cProfile	#just for profiling
import pstats	#just for profiling

import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import messagebox
from tkinter import simpledialog
from tkinter import filedialog
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.ticker as ticker
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import mplcursors
import os
import mpld3

########################################################################
class PulseSequenceAnalyzerApp:
    def __init__(self, master):
        self.master = master
        master.title("Pulse Sequence Analyses")
        master.bind('<Return>', self.handle_enter)
        self.master.pack_propagate(0)  # Prevent resizing based on content
        self.master.grid_rowconfigure(0, weight=2)  # Allow row 0 to expand when window is resized
        self.master.grid_rowconfigure(1, weight=1)  # Allow row 1 to expand when window is resized
        self.master.grid_rowconfigure(2, weight=2)  # Allow row 2 to expand when window is resized
        self.master.grid_rowconfigure(3, weight=2)  # Allow row 3 to expand when window is resized
        self.master.grid_rowconfigure(4, weight=2)  # Allow row 4 to expand when window is resized
        self.master.grid_rowconfigure(5, weight=2)  # Allow row 5 to expand when window is resized
        self.master.grid_rowconfigure(6, weight=1)  # Allow row 6 to expand when window is resized
        self.master.grid_rowconfigure(7, weight=1)  # Allow row 7 to expand when window is resized
        self.master.grid_columnconfigure(0, weight=1)  # Allow column 0 to expand when window is resized
        self.master.grid_columnconfigure(1, weight=2)  # Allow column 1 to expand when window is resized
        self.master.grid_columnconfigure(2, weight=2)  # Allow column 2 to expand when window is resized
        
        self.PS = "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR36020kHz_30B1_rf10kHz/pulse1000.bruker"
        self.maximumAmplitude = 10000.0
        self.Vector_Length = 1.0
        self.T = 0.5
        self.InpoFact = 5
        self.x_Expand = 0.0
        self.Offset = 0.0
        self.Scaling = 100.0
        self.calculationMethod = 1
        self.initialVector = np.array([0,0,1])
        
        self.documentation_str = """Documentation:
Start window:
Insert complete path of pulse sequence (e.g. /Users/leon/Documents/MATLAB/Pulssequenzen/BIBOP_sorted_20kHz_noB1_rf10kHz/pulse0015.bruker)
 and the correct parameter values into the edit fields to the left. Press refresh to load the AHT-Curve, its Curvature, Torsion and the Path (projectory of the vector end). Use the Sliders or the associated edit fields to change the offset or the pulse Amplitude and click play to automatically change the value.
Menu:
choose the desired calculation Method and the initial Condition (Bloch Vector)
Curve Window:
Plot of the AHT-Curve and its xy-plane projection. Updates through pressing the refresh or  play button or changing the offset or scaling value.
Curvature Window:
Plot of the Curvature along the Curve. Updates through pressing the refresh or play button or changing the offset or scaling value.
Torsion Window:
Plot of the Torsion along the Curve. Updates through pressing the refresh or play button or changing the offset or scaling value. The Torsion plot is only accurate for inpoFact = 1.
Path Window:
Plot of the Trajectory of the Vector Ends. Updates through pressing the refresh or play button or changing the offset or scaling value.
Stability window:
Select if the pulse sequence is a state to state or an universal rotation pulse. Then choose the offset range (in kHz, 0kHz is in center; e.g. 20kHz->-10kHz - 10kHz)   and pulse amplitude range (in %, 100 is in center; e.g. 20%->90% - 110%)  and click calculate. The plot displays the accuracy/quality (/rotation angle if angle is selected) of the pulse sequence as function of the offset and amplitude.
PS Window:
Displays pulse sequence as text or as plot.
Options Window:
1.Animation
Animate a video of AHT-Curves with offset or pulse amplitude as changing parameter in .avi format. Insert a complete path as saving location (e.g. /Users/leon/Documents/MATLAB/PulseImages/001)
2.createQualityImages
Create an Image which displays the stability of pulse sequences with increasing lengths in dependence of the offset per color coding. To do so insert the path of the folder which contains the pulse sequences into the directory field and the name the Image shall be saved under in the 'save under' field. After the Image is calculated it will also open in a figure, which can also be saved. The parameter Offset Range defines in which Offset Range (in kHz, 0kHz is in center; e.g. 20kHz->-10kHz - 10kHz) the curve quality shall be calculated. The field 'class' sets if the pulses are UR or SS pulses. Resolution sets how many pixels per kHz Range per pulse shall be calculated. Pulse Angle isn't yet functional. The calculate Image feature till now just works with all State to State pulses and 0° /360° UR pulses. With the default settings calculating an Image from a folder with 150 pulse sequences can take up to half an hour."""

        self.user_interaction = False
        self.disp_values = False       
        self.update_count = 0
        self.export_count = 0
        self.export_dir_count = 0
        self.qualityImages_export_count = 0
        self.play_amplitude_bool = False
        self.play_offset_bool = False
        
        np.set_printoptions(threshold=np.inf)
               
        # Main Frame
        #self.master = tk.Frame(master)
        #self.master.grid(row = 0, column =0)
        
        # tabs
        self.tabControl = ttk.Notebook(self.master)
        self.tabControl.grid(row = 0, column = 2, pady=15, rowspan = 5,sticky="NSEW")
        
		# Enable interactive mode
        #plt.ion()
        # Enable data picking
        #mplcursors.cursor(hover=True)
        # Axes for plotting AHT-Curve
        self.fig1 = plt.figure(figsize=(5, 3))
        self.ax1 = self.fig1.add_subplot(111, projection='3d')
        self.ax1.set_xlabel('X', fontsize=8)
        self.ax1.set_ylabel('Y', fontsize=8)
        self.ax1.set_zlabel('Z', fontsize=8)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.tabControl)
        self.canvas1.get_tk_widget().grid(row = 0, column = 0,sticky="NSEW")
        self.tabControl.add(self.canvas1.get_tk_widget(), text ='AHT-Curve') 
        
        # Axes for plotting Trajectory
        self.fig2 = plt.figure(figsize=(5, 3))
        self.ax2 = self.fig2.add_subplot(111, projection='3d')
        self.ax2.set_xlabel('X', fontsize=8)
        self.ax2.set_ylabel('Y', fontsize=8)
        self.ax2.set_zlabel('Z', fontsize=8)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.tabControl)
        self.canvas2.get_tk_widget().grid(row = 0, column = 0,sticky="NSEW")
        self.tabControl.add(self.canvas2.get_tk_widget(), text ='Trajectory') 
        
        # Axes for plotting Curvature
        self.fig3 = plt.figure(figsize=(5, 3))
        self.ax3 = self.fig3.add_subplot(111)
        self.ax3.set_xlabel('Curvature', fontsize=8)
        self.ax3.set_ylabel('Arclength', fontsize=8)
        self.canvas3 = FigureCanvasTkAgg(self.fig3, master=self.tabControl)
        self.canvas3.get_tk_widget().grid(row = 0, column = 0,sticky="NSEW")
        self.tabControl.add(self.canvas3.get_tk_widget(), text ='Curvature') 
        
        # Axes for plotting Torsion
        self.fig4 = plt.figure(figsize=(5, 3))
        self.ax4 = self.fig4.add_subplot(111)
        self.ax4.set_xlabel('Torsion', fontsize=8)
        self.ax4.set_ylabel('Arclength', fontsize=8)
        self.canvas4 = FigureCanvasTkAgg(self.fig4, master=self.tabControl)
        self.canvas4.get_tk_widget().grid(row = 0, column = 0,sticky="NSEW")
        self.tabControl.add(self.canvas4.get_tk_widget(), text ='Torsion') 
        
        # Frame for stability calculation
        self.stability_frame = tk.Frame(master=self.tabControl)
        self.stability_frame.grid(row = 0, column =0)
        self.tabControl.add(self.stability_frame, text ='Stability') 
        self.stability_frame.grid_rowconfigure(0, weight=0) 
        self.stability_frame.grid_rowconfigure(1, weight=1) 
        self.stability_frame.grid_columnconfigure(0, weight=1)  
        self.stability_frame.grid_columnconfigure(1, weight=1)  
        self.stability_frame.grid_columnconfigure(2, weight=1)
        self.stability_frame.grid_columnconfigure(3, weight=1)  
        self.stability_frame.grid_columnconfigure(4, weight=1)  
			#Radio buttons 
        self.selected_stability_option = tk.StringVar()
        self.stability_options_frame = tk.Frame(self.stability_frame)
        self.stability_options_frame.grid(row=0, column=0, columnspan=2)
        self.stability_option1 = tk.Radiobutton(self.stability_options_frame, text="PP", variable=self.selected_stability_option, value=1)
        self.stability_option1.grid(row=0, column=0)
        self.stability_option2 = tk.Radiobutton(self.stability_options_frame, text="UR", variable=self.selected_stability_option, value=2)
        self.stability_option2.grid(row=0, column=1)
        self.stability_option3 = tk.Radiobutton(self.stability_options_frame, text="Angle", variable=self.selected_stability_option, value=3)
        self.stability_option3.grid(row=0, column=2)
        self.selected_stability_option.set(1)
			#Amplitude Range
        self.stability_amplitude_range_label =	tk.Label(self.stability_frame, text = "Amplitude range (%):")
        self.stability_amplitude_range_label.grid(row = 0, column = 2, sticky = "W")
        self.stability_amplitude_range_text = tk.Entry(self.stability_frame, width = 8)
        self.stability_amplitude_range_text.grid(row = 0, column = 2, sticky = "E")	
        self.stability_amplitude_range_text.insert(0, "20")	
			#Offset Range
        self.stability_offset_range_label =	tk.Label(self.stability_frame, text = "Offset range (kHz):")
        self.stability_offset_range_label.grid(row = 0, column = 3, sticky = "W")
        self.stability_offset_range_text = tk.Entry(self.stability_frame, width = 8)
        self.stability_offset_range_text.grid(row = 0, column = 3, sticky = "E")
        self.stability_offset_range_text.insert(0, "20")
			#calculate Button
        self.calculate_stability_button = tk.Button(self.stability_frame, text="Calculate", command=self.calculate_stability)
        self.calculate_stability_button.grid(row = 0, column = 4, sticky = "E")
			#Axes for surface plot
        self.fig5 = plt.figure(figsize=(5, 3))
        self.ax5 = self.fig5.add_subplot(111, projection='3d')
        self.canvas5 = FigureCanvasTkAgg(self.fig5, master=self.stability_frame)
        self.canvas5.get_tk_widget().grid(row = 1, column = 0, columnspan = 5, sticky="NSEW")
        	
        # Frame for pulse sequence
        self.PS_frame = tk.Frame(master=self.tabControl)
        self.PS_frame.grid(row = 0, column =0)
        self.tabControl.add(self.PS_frame, text ='Pulse_Sequence') 
        self.PS_frame.grid_rowconfigure(0, weight=1) 
        self.PS_frame.grid_columnconfigure(0, weight=1)  
        self.PS_mat_text_area = scrolledtext.ScrolledText(self.PS_frame)
        self.PS_mat_text_area.insert(0.0, "Pulse sequence will display here")
        self.PS_mat_text_area.configure(state = "disabled")
        self.PS_mat_text_area.grid(row =0, column = 0, sticky = "NSEW") 
        
        # Frame for other option
        self.options_frame = tk.Frame(master=self.tabControl)
        self.options_frame.grid(row = 0, column =0)
        self.tabControl.add(self.options_frame, text ='options') 
        self.options_frame.grid_rowconfigure(0, weight = 1)
        self.options_frame.grid_rowconfigure(1, weight = 1)
        self.options_frame.grid_rowconfigure(2, weight = 1)
        self.options_frame.grid_rowconfigure(3, weight = 1)
        self.options_frame.grid_rowconfigure(4, weight = 1)
        self.options_frame.grid_rowconfigure(5, weight = 1)
        self.options_frame.grid_rowconfigure(6, weight = 1)
        self.options_frame.grid_rowconfigure(7, weight = 1)
        self.options_frame.grid_rowconfigure(8, weight = 1)
        self.options_frame.grid_rowconfigure(9, weight = 1)
        self.options_frame.grid_columnconfigure(0, weight = 1)
        self.options_frame.grid_columnconfigure(1, weight = 1)
        self.options_frame.grid_columnconfigure(2, weight = 1)
        self.options_frame.grid_columnconfigure(3, weight = 1)
        self.options_frame.grid_columnconfigure(4, weight = 1)
        self.options_frame.grid_columnconfigure(5, weight = 1)   
        	#Export option
        self.export_data_label = tk.Label(self.options_frame, text="Export Data")
        self.export_data_label.grid(row = 0, column = 2, padx = 10, pady = 10, sticky = "E")
        self.export_data_path_label = tk.Label(self.options_frame, text="Path:")
        self.export_data_path_label.grid(row = 1, column = 0, padx = 10, pady = 10)
        self.export_data_path_text =  tk.Entry(self.options_frame, width = 40)
        self.export_data_path_text.grid(row = 1, column= 1, padx = 10, pady = 10)
        self.export_data_path_text.insert(0, "/Users/leon/Desktop/Physik/Glaser")
        self.export_browse_directory_button =  tk.Button(self.options_frame,  text="Browse", command=self.export_browse_directory)
        self.export_browse_directory_button.grid(row = 1, column = 2, padx = 10, pady = 10)
        self.export_data_path_button =  tk.Button(self.options_frame,  text="Export Data", command=self.export_data_interface)
        self.export_data_path_button.grid(row = 1, column = 3, padx = 10, pady = 10)    
        self.export_data_path_button =  tk.Button(self.options_frame,  text="Export Dir Data", command=self.export_dir_data_interface)  #export data of all the pulses nin the same directory as the active pulse
        self.export_data_path_button.grid(row = 1, column = 4, padx = 10, pady = 10)    
			######## Quality Images stuff
        self.createQualityImages_label = tk.Label(self.options_frame, text = "create Quality Images")
        self.createQualityImages_label.grid(row = 2, column = 2, sticky = "E")
			# Radio buttons changing Variable
        self.qualityImages_changingVariable_label = tk.Label(self.options_frame, text = "changing Variable:")
        self.qualityImages_changingVariable_label.grid(row = 3, column = 0, sticky = "E")	
        self.selected_qualityImages_changingVariable = tk.StringVar()
        self.qualityImages_changingVariable_frame = tk.Frame(self.options_frame)
        self.qualityImages_changingVariable_frame.grid(row=3, column=1, columnspan=2, sticky = "W")
        self.qualityImages_changingVariable1 = tk.Radiobutton(self.qualityImages_changingVariable_frame, text="Amplitude", variable=self.selected_qualityImages_changingVariable, value=1)
        self.qualityImages_changingVariable1.grid(row=0, column=0)
        self.qualityImages_changingVariable2 = tk.Radiobutton(self.qualityImages_changingVariable_frame, text="Offset", variable=self.selected_qualityImages_changingVariable, value=2)
        self.qualityImages_changingVariable2.grid(row=0, column=1)
        self.selected_qualityImages_changingVariable.set(1)
			# Radio buttons calc Type
        self.qualityImages_calcType_label = tk.Label(self.options_frame, text = "calculation Type:")
        self.qualityImages_calcType_label.grid(row = 3, column = 3, sticky = "E")	
        self.selected_qualityImages_calcType = tk.StringVar()
        self.qualityImages_calcType_frame = tk.Frame(self.options_frame)
        self.qualityImages_calcType_frame.grid(row=3, column=4, columnspan=2, sticky = "W")
        self.qualityImages_calcType1 = tk.Radiobutton(self.qualityImages_calcType_frame, text="SS", variable=self.selected_qualityImages_calcType, value=1)
        self.qualityImages_calcType1.grid(row=0, column=0)
        self.qualityImages_calcType2 = tk.Radiobutton(self.qualityImages_calcType_frame, text="UR", variable=self.selected_qualityImages_calcType, value=2)
        self.qualityImages_calcType2.grid(row=0, column=1)
        self.selected_qualityImages_calcType.set(1)
			# Other parameters (Range, Umax, Resolution, Time/pulse, 
        self.qualityImages_Range_label = tk.Label(self.options_frame, text = "Range [kHz]:")
        self.qualityImages_Range_label.grid(row = 4, column = 0, sticky = "E")
        self.qualityImages_Range_text = tk.Entry(self.options_frame, width = 15)
        self.qualityImages_Range_text.grid(row = 4, column = 1, sticky = "W")
        self.qualityImages_Range_text.insert(0, "40")
        self.qualityImages_Umax_label = tk.Label(self.options_frame, text = "maximum Amplitude [kHz]:")
        self.qualityImages_Umax_label.grid(row = 4, column = 2, sticky = "E")
        self.qualityImages_Umax_text = tk.Entry(self.options_frame, width = 15)
        self.qualityImages_Umax_text.grid(row = 4, column = 3, sticky = "W")
        self.qualityImages_Umax_text.insert(0, "10")
        self.qualityImages_Resolution_label = tk.Label(self.options_frame, text = "Resolution [p/kHz]:")
        self.qualityImages_Resolution_label.grid(row = 4, column = 4, sticky = "E")
        self.qualityImages_Resolution_text = tk.Entry(self.options_frame, width = 15)
        self.qualityImages_Resolution_text.grid(row = 4, column = 5, sticky = "W")
        self.qualityImages_Resolution_text.insert(0, "5")
        self.qualityImages_TimeperPuls_label = tk.Label(self.options_frame, text = "Time/puls [µs]:")
        self.qualityImages_TimeperPuls_label.grid(row = 5, column = 0, sticky = "E")
        self.qualityImages_TimeperPuls_text = tk.Entry(self.options_frame, width = 15)
        self.qualityImages_TimeperPuls_text.grid(row = 5, column = 1, sticky = "W")
        self.qualityImages_TimeperPuls_text.insert(0, "0.5")
			# choose  directory
        self.qualityImages_input_dir_label = tk.Label(self.options_frame, text="Input Directory:")
        self.qualityImages_input_dir_label.grid(row = 6, column = 0, sticky = "E")
        self.qualityImages_input_dir_text =  tk.Entry(self.options_frame, width = 40)
        self.qualityImages_input_dir_text.grid(row = 6, column= 1, columnspan = 2, sticky ="W")
        self.qualityImages_input_dir_text.insert(0, "/Users/leon/Desktop/Physik/Glaser/Analyse_und_Visualisierung_von_robusten_Kontrollpulsen/Pulssequenzen/UR_Pulse/UR_ohne_B1_robustness_20_kHz/UR360_20kHz_noB1_rf10kHz(new_loop_sorted)")
        self.qualityImages_input_dir_button =  tk.Button(self.options_frame,  text="Browse", command=self.qualityImages_input_browse_directory)
        self.qualityImages_input_dir_button.grid(row = 6, column = 2)        
        self.qualityImages_output_dir_label = tk.Label(self.options_frame, text="Output Directory:")
        self.qualityImages_output_dir_label.grid(row = 7, column = 0, sticky = "E")
        self.qualityImages_output_dir_text =  tk.Entry(self.options_frame, width = 40)
        self.qualityImages_output_dir_text.grid(row = 7, column= 1, columnspan = 2, sticky ="W")
        self.qualityImages_output_dir_text.insert(0, "/Users/leon/Desktop/Physik/Glaser/Bachelor_Thesis/images/Quality images")
        self.qualityImages_output_dir_button =  tk.Button(self.options_frame,  text="Browse", command=self.qualityImages_output_browse_directory)
        self.qualityImages_output_dir_button.grid(row = 7, column = 2)   
			# calculate Button
        self.qualityImages_calculate_button =  tk.Button(self.options_frame,  text="Calculate", command=self.calc_qualityImages)
        self.qualityImages_calculate_button.grid(row = 8, column = 3, stick = "E")    	
        
        # Menu Button
        self.menu_button = tk.Button(self.master, text="Menu", command=self.open_menu)
        self.menu_button.grid(row = 0, column = 0, padx=20, pady=10,sticky="NW")
        
        # Documentation Button
        self.doc_button = tk.Button(self.master, text="Doc", command=self.doc)
        self.doc_button.grid(row=0, column = 1, padx=2, pady=10, sticky="NW")

        # Text Edit Field for Pulse Sequence
        self.pulse_sequence_label = tk.Label(self.master, text="Pulse Sequence:")
        self.pulse_sequence_label.grid(row = 1, column = 0, padx=10, pady=5)
        self.pulse_sequence_text = tk.Entry(self.master, width=20)
        self.pulse_sequence_text.insert(0, self.PS) 
        self.pulse_sequence_text.grid(row = 1, column = 1, padx=5, pady=5, sticky = "W")
        self.PS_browse_directory_button =  tk.Button(self.master,  text="Browse", command=self.PS_browse_file)
        self.PS_browse_directory_button.grid(row = 1, column = 1, padx = 10, pady = 10, sticky = "E")      

        # Parameter Adjustments
        self.param_frame = tk.Frame(self.master)
        self.param_frame.grid(row=2, column = 0, columnspan =2, padx=10, pady=10)
        self.param_frame.grid_rowconfigure(0, weight=1)  # Allow row 0 to expand when window is resized
        self.param_frame.grid_rowconfigure(1, weight=1)  # Allow row 1 to expand when window is resized
        self.param_frame.grid_rowconfigure(2, weight=1)  # Allow row 2 to expand when window is resized
        self.param_frame.grid_rowconfigure(3, weight=1)  # Allow row 3 to expand when window is resized
        self.param_frame.grid_rowconfigure(4, weight=1)  # Allow row 4 to expand when window is resized
        self.param_frame.grid_rowconfigure(5, weight=1)  # Allow row 5 to expand when window is resized
        self.param_frame.grid_rowconfigure(6, weight=1)  # Allow row 6 to expand when window is resized
        self.param_labels = ["Pulse Amplitude [kHz]:", "Vector Length:", "Time/pulse [µs]:", "Inpofact:", "x_Expand:", "Offset [kHz]:", "Scaling %:"]
        self.param_initials = ["10", "1", "0.5", "5", "0", "0", "100"]
        self.param_entries = []
        for i, label in enumerate(self.param_labels):
            lbl = tk.Label(self.param_frame, text=label)
            lbl.grid(row = i, column = 0, pady=2)
            entry = tk.Entry(self.param_frame, width=10)
            entry.insert(0, self.param_initials[i])  
            entry.grid(row = i, column = 1, pady=2)
            self.param_entries.append(entry)
	

        # Curve Data Text Area
        self.curve_data_label = tk.Label(self.master, text="Curve Data:")
        self.curve_data_label.grid(row=3, column = 0, columnspan = 2, padx=10, pady=5)
        self.curve_data_text = scrolledtext.ScrolledText(self.master, width=45, height=15)
        self.curve_data_text.configure(state ='disabled')
        self.curve_data_text.grid(row=4, column = 0, columnspan = 2, padx=10, pady=4, sticky="NSEW")

        # Infos and Errors Text Area
        self.info_error_label = tk.Label(self.master, text="Infos and Errors:")
        self.info_error_label.grid(row=5, column = 0, columnspan = 2, padx=10, pady=5)
        self.info_error_text = scrolledtext.ScrolledText(self.master, width=60, height=8)
        self.info_error_text.configure(state ='disabled')
        self.info_error_text.grid(row=6, column = 0, columnspan = 2, padx=10, pady=4, sticky="NSEW")

        # Refresh Button
        self.refresh_button = tk.Button(self.master, text="Refresh", command=self.refresh)
        self.refresh_button.grid(row=0, column = 1, padx=2, pady=10, sticky="NE")

        # Sliders for Amplitude and Offset
        self.slider_frame = tk.Frame(self.master)
        self.slider_frame.grid(row=5, column = 2, rowspan = 3, padx=10, pady=2)
        self.amplitude_label = tk.Label(self.slider_frame, text="Amplitude %:")
        self.amplitude_label.grid(row=0, column = 0, padx=5, pady=5)
        self.amplitude_slider = tk.Scale(self.slider_frame, length = 600, width = 20, from_=0, to=200, orient=tk.HORIZONTAL, command=self.on_amplitudeSlider_change)
        #self.amplitude_slider = tk.Scale(self.slider_frame, length = 600, width = 20, from_=0, to=200, orient=tk.HORIZONTAL)
        self.amplitude_slider.set(self.Scaling)
        self.amplitude_slider.grid(row=0, column = 1, padx=5, pady=5)
        self.amplitude_slider_play_button = tk.Button(self.slider_frame, text = "play", command=self.play_amplitude)
        self.amplitude_slider_play_button.grid(row = 0, column = 2, padx=10, pady=5)
        self.offset_label = tk.Label(self.slider_frame, text="Offset kHz:")
        self.offset_label.grid(row=1, column = 0, padx=5, pady=5)
        self.offset_slider = tk.Scale(self.slider_frame, length = 600, width = 20, from_=-20, to=20, orient=tk.HORIZONTAL, command=self.on_offsetSlider_change)
        #self.offset_slider = tk.Scale(self.slider_frame, length = 600, width = 20, from_=-20, to=20, orient=tk.HORIZONTAL)
        self.offset_slider.set(self.Offset/1000)
        self.offset_slider.grid(row=1, column = 1, padx=5, pady=5)
        self.offset_slider_play_button = tk.Button(self.slider_frame, text = "play", command=self.play_offset)
        self.offset_slider_play_button.grid(row = 1, column = 2, padx=10, pady=5)

        
	#########################
 
    def open_menu(self):
        menu_window = MenuWindow(self.master, self)

    def play_amplitude(self):
        self.play_amplitude_bool = not self.play_amplitude_bool
        print("Not implemeted")
        #while self.play_amplitude_bool:
            #print("Not implemeted")

    def play_offset(self):
        self.play_offset_bool = not self.play_offset_bool
        print("Not implemeted")

    def calculate_stability(self):
        self.param_entries[3].delete(0,1000)
        self.param_entries[3].insert(0,"1")
        CM, VM, PS_mat, totalRot, phi, numberOfPulses, Axy, Axz, Ayz, arc_length, curvature, torsion, integrated_curvature, integrated_torsion, integrated_absolut_torsion, avg_curvature, avg_torsion = self.update()
        #profiler = cProfile.Profile()
        #profiler.enable()
        X,Y,Z,quality = calculate_curve_stability(PS_mat, 
                                                  self.T, 
                                                  self.Vector_Length, 
                                                  self.maximumAmplitude, 
                                                  scalingRange_percent=float(self.stability_amplitude_range_text.get()), 
												  offsetRange_kHz=float(self.stability_offset_range_text.get()), 
                                                  stabilityCalculationMethod = int(self.selected_stability_option.get()), 
												  initialVector = self.initialVector)
        #profiler.disable()
        self.ax5.clear()
        self.ax5.plot_surface(X, Y, Z, cmap='viridis', edgecolor='black', linewidth=0.1)
        self.canvas5.draw()
        self.text_area_set(text_area = self.info_error_text, text_str = "Quality: " + str(quality), reset_bool = 0)
        #stats = pstats.Stats(profiler).sort_stats('cumulative')
        #stats.print_stats()
        #print("Calculate stability function not yet fully implemented!!")
		 
    def doc(self):
        doc_window = tk.Toplevel(self.master)
        doc_window.title("Documentation")
        doc_window.pack_propagate(0)  # Prevent resizing based on content
        doc_window.grid_rowconfigure(0, weight=1) 
        doc_window.grid_columnconfigure(0, weight=1) 
        doc_text_area = scrolledtext.ScrolledText(doc_window, width=100, height=40)
        doc_text_area.grid(row=0, column=0, sticky ="NSEW")
        doc_text_area.insert(0.0, self.documentation_str)
        self.curve_data_text.configure(state ='disabled')
        
    def PS_browse_file(self):
        # Open file dialog to choose directory
        chosen_file = filedialog.askopenfilename()      
        # If directory chosen, update the Entry widget with the chosen directory
        if chosen_file:
            self.pulse_sequence_text.delete(0, tk.END)  # Clear any existing text
            self.pulse_sequence_text.insert(tk.END, chosen_file)
            
    def export_browse_directory(self):
        # Open file dialog to choose directory
        chosen_directory = filedialog.askdirectory()      
        # If directory chosen, update the Entry widget with the chosen directory
        if chosen_directory:
            self.export_data_path_text.delete(0, tk.END)  # Clear any existing text
            self.export_data_path_text.insert(tk.END, chosen_directory)
   
    def qualityImages_input_browse_directory(self):
        # Open file dialog to choose directory
        chosen_directory = filedialog.askdirectory()      
        # If directory chosen, update the Entry widget with the chosen directory
        if chosen_directory:
            self.qualityImages_input_dir_text.delete(0, tk.END)  # Clear any existing text
            self.qualityImages_input_dir_text.insert(tk.END, chosen_directory)

    def qualityImages_output_browse_directory(self):
        # Open file dialog to choose directory
        chosen_directory = filedialog.askdirectory()      
        # If directory chosen, update the Entry widget with the chosen directory
        if chosen_directory:
            self.qualityImages_output_dir_text.delete(0, tk.END)  # Clear any existing text
            self.qualityImages_output_dir_text.insert(tk.END, chosen_directory)
 
    def calc_qualityImages(self):
        PSQM = calculate_pulse_sequence_quality_images(self,
                                                       self.qualityImages_input_dir_text.get(), 
                                                       float(self.qualityImages_Range_text.get())*1000, 
                                                       float(self.qualityImages_TimeperPuls_text.get())*(10**(-6)), 
													   float(self.qualityImages_Umax_text.get())*1000, 
                                                       self.initialVector, 
													   int(self.selected_qualityImages_calcType.get()), 
                                                       int(self.selected_qualityImages_changingVariable.get()),
													   float(self.qualityImages_Resolution_text.get()))
        while os.path.exists(self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count)):
            self.export_count += 1
        os.mkdir(self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count))
        np.savetxt(self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count)+"/QualityMatrix.csv", PSQM, delimiter=",")
        #fig, ax = plt.subplots()
        #cax = ax.imshow(PSQM, cmap='jet')
        #cbar = fig.colorbar(cax)
        #plt.savefig(self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count)+"/QualityImage.png")
        plt.imsave(self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count)+"/QualityImage.png", PSQM, cmap='jet')
        calc_qualityImages_str  = "Quality Image exported to: "+  self.qualityImages_output_dir_text.get()+"/QualityImage"+str(self.export_count)
        self.text_area_set(text_area = self.info_error_text, text_str = calc_qualityImages_str, reset_bool = 0)
		                                     
    def export_data_interface(self):
        while os.path.exists(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)):
            self.export_count += 1
        info_window = tk.Toplevel(self.master)
        info_window.title("Info Window")
        info_window_label = tk.Label(info_window, text= "Data will be exported to directory: " + self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count))
        info_window_label.grid(row = 0, column = 0, columnspan = 2)
        info_window_ok_button = tk.Button(info_window, text="OK", command=lambda: self.info_window_ok_button_pressed(info_window, 0))
        info_window_ok_button.grid(row = 1, column = 0)
        info_window_cancel_button = tk.Button(info_window, text="cancel", command=lambda: self.cancel_export(info_window))
        info_window_cancel_button.grid(row = 1, column = 1)

    def export_dir_data_interface(self):
        while os.path.exists(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)):
            self.export_dir_count += 1
        info_window = tk.Toplevel(self.master)
        info_window.title("Info Window")
        info_window_label = tk.Label(info_window, text= "Directory Data will be exported to directory: " + self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count))
        info_window_label.grid(row = 0, column = 0, columnspan = 2)
        info_window_ok_button = tk.Button(info_window, text="OK", command=lambda: self.info_window_ok_button_pressed(info_window, 1))
        info_window_ok_button.grid(row = 1, column = 0)
        info_window_cancel_button = tk.Button(info_window, text="cancel", command=lambda: self.cancel_export(info_window))
        info_window_cancel_button.grid(row = 1, column = 1)

    def info_window_ok_button_pressed(self, window, type):
        window.destroy()
        if(type == 0):          #just export this pulse
            self.export_data(type=type)
        if(type == 1):          #export all sequences from directory
            os.mkdir(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count))
            dirname = os.path.dirname(self.PS)
            files = os.listdir(dirname)
            sorted_files = sorted(files)
            PSnames = [f for f in sorted_files if 'bruker' in f]
            PSnumber = len(PSnames)
            for n in range(PSnumber):
                PS = os.path.join(dirname, PSnames[n])
                self.pulse_sequence_text.delete(0, tk.END)  # Clear any existing text
                self.pulse_sequence_text.insert(tk.END, PS)
                self.export_data(type=type, exp_str = PSnames[n].replace(".bruker",""))
        self.export_count = 0
        self.export_dir_count = 0

    def export_data(self,type,exp_str="empty"):
        CM, VM, PS_mat, totalRot, phi, numberOfPulses, Axy, Axz, Ayz, arc_length, curvature, torsion, integrated_curvature, integrated_torsion, integrated_absolut_torsion, avg_curvature, avg_torsion = self.update()
        export_string = (f"Pulse Sequence {self.PS} ({round(1E3 * totalRot) / 1E3}°)\n"
	    f"Time per Pulse: {self.T * 10**6} µs\n"
		f"Vector Length: {self.Vector_Length}\n"
		f"Maximum Amplitude: {self.maximumAmplitude / 1000} kHz\n"
		f"Offset: {self.Offset / 1000} kHz\n"
		f"InpoFact: {self.InpoFact}\n"
		f"x Expand: {self.x_Expand}\n"
		f"Calculation Method: {self.calculationMethod}\n"
		f"Initial Vector: {self.initialVector}\n"
		f"Angle to x-Axis= {phi}°;\n"
		f"number of pulses = {numberOfPulses};\n"
		f"Axy= {Axy};\n"
		f"Axz= {Axz},\n"
		f"Ayz= {Ayz};\n"
		f"Integrated Curvature= {integrated_curvature};\n"
		f"Integrated Torsion= {integrated_torsion};\n"
		f"Integrated abs Torsion= {integrated_absolut_torsion};")
        if(type == 0):
            os.mkdir(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count))
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/AHT_Curve.csv", CM, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/Trajectory_Curve.csv", VM, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/Pulse_Sequence.csv", PS_mat, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/curvature-arc_length.csv", np.column_stack((arc_length, curvature)), delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/torsion-arc_length.csv", np.column_stack((arc_length, torsion)), delimiter=",")
            with open(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/Curve_data.txt", 'w') as file:
                file.write(export_string)
            mpld3.save_html(self.fig1, self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Data"+str(self.export_count)+"/AHT_Curve_plot.html")	#######TODO
            export_data_interface_str = "data exported to: "+  self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_count)
            self.text_area_set(text_area = self.info_error_text, text_str = export_data_interface_str, reset_bool = 0)
        if(type == 1):
            os.mkdir(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str)
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/AHT_Curve.csv", CM, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/Trajectory_Curve.csv", VM, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/Pulse_Sequence.csv", PS_mat, delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/curvature-arc_length.csv", np.column_stack((arc_length, curvature)), delimiter=",")
            np.savetxt(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/torsion-arc_length.csv", np.column_stack((arc_length, torsion)), delimiter=",")
            with open(self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/Curve_data.txt", 'w') as file:
                file.write(export_string)
            mpld3.save_html(self.fig1, self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str+"/AHT_Curve_plot.html")	#######TODO
            export_data_interface_str = "data exported to: "+  self.export_data_path_text.get()+"/PulseSequenceAnalyzer_Dir_Data"+str(self.export_dir_count)+"/"+exp_str
            self.text_area_set(text_area = self.info_error_text, text_str = export_data_interface_str, reset_bool = 0)
        
    def cancel_export(self, window):
        window.destroy()
        self.export_count = 0
        self.export_dir_count = 0
        
    def refresh(self):
        # Placeholder for refresh functionality
        self.update()
        
    def on_amplitudeSlider_change(self, amplitudeSlider_Value):
        if self.user_interaction:
            self.param_entries[6].delete(0,1000)
            self.param_entries[6].insert(0, amplitudeSlider_Value) 
            self.update()
        self.user_interaction = True
        
    def on_offsetSlider_change(self, offsetSlider_Value):
        if self.user_interaction:
            self.param_entries[5].delete(0,1000)
            self.param_entries[5].insert(0, offsetSlider_Value) 
            self.update()
        self.user_interaction = True    	    	
				
    def handle_enter(self, event):
        # Diese Methode wird aufgerufen, wenn die Eingabetaste gedrückt wird
        self.update()
    
    def scientific_formatter(self, value, pos):
        if np.isclose(value, 0.0):
            return "0"
        magnitude = int(np.floor(np.log10(abs(value))))
        formatted_value = value / (10 ** magnitude)
        if magnitude >= 0:
            return f'{formatted_value:.1f}e{magnitude}'
        else:
            return f'{formatted_value:.1f}e-{abs(magnitude)}'
        
    def plot3DCurve(self, CurveX, CurveY, CurveZ, axes, canvas, Title=""):
        axes.clear()
        axes.plot3D(CurveX, CurveY, CurveZ, color='blue', linewidth=0.9)
        axes.grid(False)
        axes.tick_params(axis='both', which='major', labelsize=8)  # Set font size of tick labels
        axes.tick_params(axis='both', which='minor', labelsize=6)  # Set font size of minor tick labels
        axes.xaxis.set_tick_params(rotation=45)  # Rotate x-axis tick labels by 45 degrees
        axes.yaxis.set_tick_params(rotation=-45)  # Rotate y-axis tick labels by -45 degrees
        axes.zaxis.set_tick_params(rotation=30)  # Rotate z-axis tick labels by 30 degrees
        axes.xaxis.set_major_locator(plt.MaxNLocator(5))  # Set maximum number of ticks for x-axis
        axes.yaxis.set_major_locator(plt.MaxNLocator(5))  # Set maximum number of ticks for y-axis
        axes.zaxis.set_major_locator(plt.MaxNLocator(5))  # Set maximum number of ticks for z-axis
        axes.xaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter)) # Set the tick formatter for x axes
        axes.yaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter)) # Set the tick formatter for y axes
        axes.zaxis.set_major_formatter(ticker.FuncFormatter(self.scientific_formatter)) # Set the tick formatter for z axes
        axes.set_title(Title, fontsize=8)
        axes.set_xlabel('X', fontsize=8)
        axes.set_ylabel('Y', fontsize=8)
        axes.set_zlabel('Z', fontsize=8)
        #axes.set_box_aspect((1, 1, 1))  
        axes.set_box_aspect((abs(max(CurveX)-min(CurveX)), abs(max(CurveY)-min(CurveY)), abs(max(CurveZ)-min(CurveZ))))
        if self.disp_values:
           mplcursors.cursor(hover=True)
        canvas.draw()
        
    def plot2DCurve(self, CurveX, CurveY, axes, canvas, xLabel, yLabel, Title=""):
        axes.clear()
        axes.plot(CurveX, CurveY, color='blue', linewidth=0.8)
        axes.grid(True)
        axes.tick_params(axis='both', which='major', labelsize=8)  # Set font size of tick labels
        axes.tick_params(axis='both', which='minor', labelsize=6)  # Set font size of minor tick labels
        axes.set_title(Title)
        axes.set_xlabel(xLabel, fontsize=8)
        axes.set_ylabel(yLabel, fontsize=8)
        #ax.set_aspect('equal')
        #axes.set_box_aspect((abs(max(CurveX)-min(CurveX)), abs(max(CurveY)-min(CurveY))))
        if self.disp_values:
           mplcursors.cursor(hover=False)
        canvas.draw()        

    def get_info_error_str(self):
        info_error_str = (
    f"Update Number: {self.update_count}\n"
    f"Time per Pulse: {self.T * 10**6} µs\n"
    f"Vector Length: {self.Vector_Length}\n"
    f"Maximum Amplitude: {self.maximumAmplitude / 1000} kHz\n"
    f"Offset: {self.Offset / 1000} kHz\n"
    f"InpoFact: {self.InpoFact}\n"
    f"x Expand: {self.x_Expand}\n"
    f"Calculation Method: {self.calculationMethod}\n"
    f"Initial Vector: {self.initialVector}")
        return info_error_str
        
    def get_curve_data_str(self, totalRot, phi, numberOfPulses, Axy, Axz, Ayz,CurvInt, TorsInt, TorsInt_abs):
        
        curve_data_str = (
    f"Pulse Sequence ({round(1E3 * totalRot) / 1E3}°)\n"
    f"Angle to x-Axis= {phi}°;\n"
    f"number of pulses = {numberOfPulses};\n"
    f"Axy= {Axy};\n"
    f"Axz= {Axz},\n"
    f"Ayz= {Ayz};\n"
    f"Integrated Curvature= {CurvInt};\n"
    f"Integrated Torsion= {TorsInt};\n"
    f"Integrated abs Torsion= {TorsInt_abs};")
        return curve_data_str
        
    def text_area_set(self, text_area, text_str, reset_bool):
        text_area.configure(state ='normal') 
        if reset_bool:
           text_area.delete(1.0,tk.END)
        text_area.insert(1.0, "\n############\n") 
        text_area.insert(1.0, text_str) 
        text_area.configure(state ='disabled')

		
    def updateValues(self):
        PS_str = self.pulse_sequence_text.get()
        Pulse_Amplitude_str= self.param_entries[0].get()
        Vector_Length_str= self.param_entries[1].get()
        T_str= self.param_entries[2].get()
        InpoFact_str= self.param_entries[3].get()
        x_Expand_str= self.param_entries[4].get()
        Offset_str= self.param_entries[5].get()
        Scaling_str= self.param_entries[6].get()
       
        Pulse_Amplitude_Numeric= float(Pulse_Amplitude_str)
        Vector_Length_Numeric= float(Vector_Length_str)
        T_Numeric= float(T_str)
        InpoFact_Numeric= int(InpoFact_str)
        x_Expand_Numeric= float(x_Expand_str)
        Offset_Numeric= float(Offset_str)
        Scaling_Numeric= float(Scaling_str)
        
        self.maximumAmplitude = Pulse_Amplitude_Numeric*Scaling_Numeric*(10**3)/100
        self.Offset = Offset_Numeric*(10**3)
        self.T = T_Numeric*(10**-6)
        self.PS = PS_str
        self.Vector_Length = Vector_Length_Numeric
        self.InpoFact = InpoFact_Numeric
        self.x_Expand = x_Expand_Numeric
        self.Scaling = Scaling_Numeric
        
        self.user_interaction = False
        self.amplitude_slider.set(self.Scaling)
        self.offset_slider.set(self.Offset/1000)

		
    def update(self):
        #print(self.disp_values)
        self.update_count += 1
        self.updateValues()
        CM, VM, PS_mat = createCurve(PulseSequence=self.PS,T=self.T,l=self.Vector_Length,maximumAmplitude = self.maximumAmplitude,offset = self.Offset, inpoFact = self.InpoFact, xExpand = self.x_Expand, calculationMethod=self.calculationMethod, initialVector=self.initialVector)
        self.plot3DCurve(CurveX = CM[:,0]+np.linspace(0,self.x_Expand,np.shape(CM)[0]).transpose(), CurveY = CM[:,1], CurveZ = CM[:,2], axes=self.ax1, canvas = self.canvas1, Title="AHT-Kurve")
        self.plot3DCurve(CurveX = VM[:,0], CurveY = VM[:,1], CurveZ = VM[:,2], axes=self.ax2, canvas = self.canvas2, Title="Trajectory")
        totalRot = angle_between_vectors(VM[0,:], VM[-1,:])
        phi = angle_with_x_axis(VM[-1,:])
        numberOfPulses = number_of_pulses(PS_mat, self.InpoFact)
        Axy = calculate_closed_curve_area_app(CM[:,[0,1]], close_curve=False)
        Axz = calculate_closed_curve_area_app(CM[:,[0,2]], close_curve=False)
        Ayz = calculate_closed_curve_area_app(CM[:,[1,2]], close_curve=False)
        arc_length, curvature, torsion, integrated_curvature, integrated_torsion, integrated_absolut_torsion, avg_curvature, avg_torsion = curvature_torsion_3d(CM)
        self.plot2DCurve(CurveX = arc_length, CurveY = curvature, axes = self.ax3, canvas = self.canvas3, xLabel = "arc_length", yLabel = "curvature", Title = "Curvature")
        self.plot2DCurve(CurveX = arc_length, CurveY = torsion, axes = self.ax4, canvas = self.canvas4, xLabel = "arc_length", yLabel = "torsion", Title = "Torsion")
        self.text_area_set(self.curve_data_text, self.get_curve_data_str(totalRot, phi, numberOfPulses, Axy, Axz, Ayz,CurvInt=integrated_curvature, TorsInt=integrated_torsion, TorsInt_abs=integrated_absolut_torsion), True)
        self.text_area_set(self.info_error_text, self.get_info_error_str(), False)
        self.PS_mat_text_area.configure(state = "normal")
        self.PS_mat_text_area.delete(1.0,tk.END)
        self.PS_mat_text_area.insert(0.0,str(PS_mat))
        self.PS_mat_text_area.configure(state = "disabled")
        return CM, VM, PS_mat, totalRot, phi, numberOfPulses, Axy, Axz, Ayz, arc_length, curvature, torsion, integrated_curvature, integrated_torsion, integrated_absolut_torsion, avg_curvature, avg_torsion

########################################################################
class MenuWindow:
    def __init__(self, master, PS_analyzer):
        self.master = master
        self.PS_analyzer = PS_analyzer
        self.menu_window = tk.Toplevel(master)
        self.menu_window.title("Menu")

        self.disp_values = tk.BooleanVar(self.master)
        self.disp_values.set(self.PS_analyzer.disp_values)
        
        # initialVector
        label_initialVec = tk.Label(self.menu_window, text="InitialVector:")
        label_initialVec.grid(row=1, column=0)
        
        self.label_x_initial = tk.Label(self.menu_window, text="x:")
        self.label_x_initial.grid(row=2, column=0)

        self.entry_x_initial = tk.Entry(self.menu_window)
        self.entry_x_initial.insert(0, self.PS_analyzer.initialVector[0]) 
        self.entry_x_initial.grid(row=2, column=1)
        
        self.label_y_initial = tk.Label(self.menu_window, text="y:")
        self.label_y_initial.grid(row=3, column=0)

        self.entry_y_initial = tk.Entry(self.menu_window)
        self.entry_y_initial.insert(0, self.PS_analyzer.initialVector[1]) 
        self.entry_y_initial.grid(row=3, column=1)
        
        self.label_z_initial = tk.Label(self.menu_window, text="z:")
        self.label_z_initial.grid(row=4, column=0)

        self.entry_z_initial = tk.Entry(self.menu_window)
        self.entry_z_initial.insert(0, self.PS_analyzer.initialVector[2]) 
        self.entry_z_initial.grid(row=4, column=1)

        # Radiobuttons mit calculationMethod
        self.label_initialVec = tk.Label(self.menu_window, text="Calculation Method:")
        self.label_initialVec.grid(row=5, column=0)
        
        self.selected_option = tk.StringVar()
        options_frame = tk.Frame(self.menu_window)
        options_frame.grid(row=6, column=0, columnspan=2)

        option1 = tk.Radiobutton(options_frame, text="Rotation matrices", variable=self.selected_option, value=1)
        option1.grid(row=6, column=0)

        option3 = tk.Radiobutton(options_frame, text="Helices", variable=self.selected_option, value=2)
        option3.grid(row=6, column=2)
        
        self.selected_option.set(self.PS_analyzer.calculationMethod)

		# Checkbox for Display values
        self.display_values_checkbox = tk.Checkbutton(self.menu_window, text="Display values", variable=self.disp_values)
        self.display_values_checkbox.grid(row=7, column=0, columnspan=2)

		# Werte sichern	
        self.save_button = tk.Button(self.menu_window, text="save", command=self.save_menu)
        self.save_button.grid(row=8, column=0, columnspan=2, pady=10)
        
        # Schließen des Menüfensters
        self.close_button = tk.Button(self.menu_window, text="Close", command=self.close_menu)
        self.close_button.grid(row=8, column=2, columnspan=2, pady=10)
		
		#reset the app to initial values
        self.reset_button = tk.Button(self.menu_window, text="Reset Analyzer", command=self.reset_analyzer)
        self.reset_button.grid(row=7, column=2, columnspan=1, pady=5)

    def reset_analyzer(self):
        # Reinitialize the PulseSequenceAnalyzer object
        self.PS_analyzer = PulseSequenceAnalyzerApp(self.master)
        self.close_menu()
        
    def save_menu(self):
        self.PS_analyzer.initialVector = np.array([float(self.entry_x_initial.get()), float(self.entry_y_initial.get()), float(self.entry_z_initial.get())])
        selected_option = self.selected_option.get()
        self.PS_analyzer.calculationMethod = int(selected_option)
        #print(self.disp_values)
        self.PS_analyzer.disp_values = self.disp_values.get()
        self.menu_window.destroy()
		
    def close_menu(self):
        self.menu_window.destroy()
        
########################################################################        
def main():
    #plt.ion()
    root = tk.Tk()
    app = PulseSequenceAnalyzerApp(root)
    root.mainloop()


