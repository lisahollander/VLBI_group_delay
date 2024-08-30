import os
import tkinter as tk
from tkinter import filedialog
from tkinter import END
from tkinter import messagebox
from tkinter import ttk
from tkinter.filedialog import asksaveasfile

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Circle
import numpy as np
import logging

from read_file import FileData
from source_model import SourceModel, SourceModelAnalytical
from plot import get_plot_dictionaries_gui


"""
Class for managing text output in GUI
"""
class TextHandler(logging.Handler):

    def __init__(self, text_widget):
        logging.Handler.__init__(self)
        self.text_widget = text_widget

    def emit(self, record):
        msg = self.format(record)
        def append():
            self.text_widget.configure(state='normal')
            self.text_widget.insert(tk.END, msg + '\n')
            self.text_widget.configure(state='disabled')
            self.text_widget.yview(tk.END)  # Autoscroll to the end
        self.text_widget.after(0, append)

class GUI(tk.Frame):

    def __init__(self, root):
        tk.Frame.__init__(self, root)
        self.root = root
        self.root.title("Source Model")
        self.root.config(bg="white")

        # Configure grid for expanding
        self.root.grid_rowconfigure(1, weight=1)  
        self.root.grid_columnconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)

        self.plot_dict_mdl = None 
        self.plot_dict_anl = None 
        self.plot_dict_clean = None

        # Create top frame
        self.top_frame = tk.Frame(self.root, width=600, height=50, bg='lightgray')
        self.top_frame.grid(row=0, column=0, columnspan=10, padx=10, pady=5, sticky='ew')

        # Create left_frame
        self.left_frame = tk.Frame(self.root, bg='lightgray')
        self.left_frame.grid(row=1, column=0, padx=10, pady=5, sticky='nsew')

        # Create right_frame
        self.right_frame = tk.Frame(self.root, bg='lightgray')
        self.right_frame.grid(row=1, column=1, padx=10, pady=5, sticky='nsew')

        # Configure left_frame grid
        self.left_frame.grid_rowconfigure(2, weight=1)
        self.left_frame.grid_columnconfigure(0, weight=1)
        self.left_frame.grid_columnconfigure(1, weight=1)
        self.left_frame.grid_columnconfigure(2, weight=1)
        self.left_frame.grid_columnconfigure(4, weight=1)

        # Configure right_frame grid
        self.right_frame.grid_rowconfigure(0, weight=0)
        self.right_frame.grid_rowconfigure(1, weight=1)
        self.right_frame.grid_columnconfigure(0, weight=1)

        # Add radio buttons for file selection
        self.file_option = tk.StringVar(value="local")
        self.local_radio = tk.Radiobutton(self.top_frame, text="Local File", variable=self.file_option, value="local", bg='lightgray', command=self.update_entry_state)
        self.local_radio.grid(row=0, column=0, padx=5, pady=0, sticky='nw')

        self.repo_radio = tk.Radiobutton(self.top_frame, text="Repo Path", variable=self.file_option, value="repo", bg='lightgray', command=self.update_entry_state)
        self.repo_radio.grid(row=0, column=1, padx=5, pady=0, sticky='nw')

        self.third_option_radio = tk.Radiobutton(self.top_frame, text="Third Option", variable=self.file_option, value="third", bg='lightgray', command=self.update_entry_state)
        self.third_option_radio.grid(row=0, column=2, padx=5, pady=0, sticky='nw')

        # Add label, entry, and buttons to left_frame
        label_file = tk.Label(self.top_frame, text="Choose a file:", bg='lightgray')
        label_file.grid(row=1, column=0, padx=5, pady=5, sticky='nw')

        # Text entry box for file path
        self.entry = tk.Entry(self.top_frame, width=30)
        self.entry.grid(row=1, column=1, padx=5, pady=5, sticky='nw')

        # Browse button
        import_button = tk.Button(self.top_frame, text="Browse", command=self.import_file)
        import_button.grid(row=1, column=2, padx=5, pady=5, sticky='nw')

        # Run button
        run_button = ttk.Button(self.top_frame, text="Run", command=self.Run)
        run_button.grid(row=1, column=3, padx=0, pady=5, sticky='nw')

        # Frame for plotting inside right_frame
        self.plot_frame = tk.Frame(self.right_frame, bg='white')
        self.plot_frame.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')

        # Initialize figure object
        self.f = Figure(figsize=(6, 8), dpi=100)

        # Initialize canvas for figures
        self.canvas = FigureCanvasTkAgg(self.f, master=self.plot_frame)
        self.canvas.get_tk_widget().pack(side='right', fill='both', expand=1)

        stop_button = ttk.Button(self.top_frame, text="Stop run", command=self.stop_run)
        stop_button.grid(row=1, column=4, padx=5, pady=5, sticky='W')

        save_button = ttk.Button(self.left_frame, text="Save log", command=self.save_log)
        save_button.grid(row=10, column=0, padx=5, pady=5, sticky='W')

        # Frame for buttons inside right_frame
        self.button_frame = tk.Frame(self.top_frame, bg='lightgray')
        self.button_frame.grid(row=0, column=5, padx=10, pady=5, sticky='ew')

        # Label above buttons
        label_result = tk.Label(self.button_frame, text="View results:", bg='lightgray')
        label_result.grid(row=2, column=5, padx=5, pady=5, sticky='s')

        # Add Analytical button
        analytical_button = ttk.Button(self.button_frame, text="Analytical", command=self.plot_analytical)
        analytical_button.grid(row=2, column=6, padx=5, pady=5, sticky='s')

        # Add Model button
        model_button = ttk.Button(self.button_frame, text="Model", command=self.plot_model)
        model_button.grid(row=2, column=7, padx=5, pady=5, sticky='s')

        # Add a Text widget for log messages to left_frame
        self.log_text = tk.Text(self.left_frame, wrap='word', state='disabled', bg='white', height=10)
        self.log_text.grid(row=2, column=0, columnspan=5, padx=5, pady=5, sticky='nsew')

        # Set up logging
        text_handler = TextHandler(self.log_text)
        logging.basicConfig(format='%(message)s', level=logging.INFO)
        self.logger = logging.getLogger('QuasarModelLog')
        self.logger.setLevel(logging.INFO)
        self.logger.addHandler(text_handler)

    """
    Function for running GUI
    """
    def Run(self):
        filePath = self.entry.get()

        self.logger.info(" #-#-#-#-#-#-#-#-# RUNNING FILE #-#-#-#-#-#-#-#-# ")
        self.logger.info(f"Running analysis on file: {filePath}")
        try:

            file_data = FileData()
            image = file_data.get_image_from_path(filePath) 
                
            source = SourceModel()
            org, mdl, anl, mdl2, anl2, anlDerivative,residuals,_ = source.process(image,file_data)

            self.plot_dict_mdl, self.plot_dict_anl = get_plot_dictionaries_gui(file_data, org, mdl2, anl2, anlDerivative,residuals)

            self.logger.info("Analysis completed successfully.")
            self.logger.info("Figures generated and displayed.")

            self.plot_model()

        except FileNotFoundError:
            self.logger.error("File not found.")
            messagebox.showerror("Error", "File not found.")
        except IOError:
            self.logger.error("Could not open the file.")
            messagebox.showerror("Error", "Could not open the file.")

    """
    Function for creating plots for gui
    """
    def make_figures(self, plot_dict):

        if not plot_dict or not plot_dict['data']:
            return
    
        self.f.clf()                            #Clear current figure
        num_images = len(plot_dict['data'])     #Number of images 
        num_cols = 2                            # Fixed number of columns (2)
        num_rows = (num_images + num_cols - 1) // num_cols  # Number of rows based on total images

        axes = self.f.subplots(num_rows, num_cols)
        axes = axes.flatten() if isinstance(axes, np.ndarray) else [axes]

        # Iterate over the data and titles in plot_dict to plot on the axes
        for i, (data, title,extent,cmap, circle,x_axis, y_axis) in enumerate(zip(plot_dict['data'], plot_dict['titles'],plot_dict['extent'],
                                                                   plot_dict['cmap'], plot_dict['circle'],plot_dict['x_axis'],plot_dict['y_axis'])):
            ax = axes[i]
            im = ax.imshow(data, cmap=cmap,extent = extent, origin='lower')  # Plot the image
        
            if circle:
                circ = Circle((plot_dict['circle_dict']['x'], plot_dict['circle_dict']['y']), plot_dict['circle_dict']['r'], facecolor='None', edgecolor='w', lw=1)
                ax.add_patch(circ)
            
            ax.set_xlabel(x_axis)
            ax.set_ylabel(y_axis)
            ax.set_title(title)
            self.f.colorbar(im, ax=ax)  # Add a colorbar to each subplot

        # Hide any unused subplots
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')
    
        # Update the canvas to display the plots
        self.f.tight_layout(pad=3.0)
        self.canvas.draw()
    
    def update_entry_state(self):
        file_type = self.file_option.get()
        if file_type == "local":
            self.entry.config(state='normal')
        elif file_type == "repo":
            self.entry.config(state='readonly')  # Make entry read-only for repo path

   
    def import_file(self):
        if self.file_option.get() == "local":
            file_path = filedialog.askopenfilename(title="Select a file",filetypes=[("Fits files","*.fits"),("All files","*.*")])
            if file_path:
                self.entry.delete(0,tk.END)
                self.entry.insert(0,file_path)
        elif self.file_option.get() == "repo":
            script_dir = os.path.dirname(os.path.abspath(__file__))
            folder_path = os.path.join(script_dir, "fits")
            print(f"Repository path: {folder_path}")  # Debugging line to check the path
            if os.path.isdir(folder_path):
                file_path = filedialog.askopenfilename(initialdir=folder_path,title="Select file", filetypes=[("Fits files","*.fits"),("All files","*.*")])
                if file_path:
                    self.entry.config(state='normal')
                    self.entry.delete(0,tk.END)
                    self.entry.insert(0,file_path)
                    self.entry.config(state='disabled')
            else:
                messagebox.showerror("Error",f"Directory does not exist: {folder_path}")

    """
    Functions for events when pressing buttons in GUI
    """
    def plot_analytical(self):
        self.make_figures(self.plot_dict_anl)

    def plot_model(self):
        self.make_figures(self.plot_dict_mdl)

    def stop_run(self):
        self.root.destroy() 
        
    def save_log(self):
        file = asksaveasfile(initialfile='Untitled.txt',
                                     defaultextension=".txt", filetypes=[("All Files", "*.*"),
                                                                         ("Text Documents", "*.txt")])
        try:
            log_content = self.log_text.get('1.0', tk.END)  # Get all text from the Text widget
            file.write(log_content)  # Write the text to the file
            file.close()
        except AttributeError:
            pass  # if user closes the browse before choosing a place to save their file
