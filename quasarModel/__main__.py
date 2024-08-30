import typer
from source_model import SourceModel
from gui import *
from tools import Header, model_several_sources, set_logger, check_if_test_source
from constants import help_real, help_directory
from plot import get_plot_dictionaries,plot_results
from clean import * 
import tkinter as tk

def main(source_path: str = typer.Option(None, '-r', clamp=True, help=help_real),
         dir_path: str = typer.Option(None, '-d', clamp=True, help=help_directory)):

    if source_path:
        set_logger(log_path='QuasarModelLog', console=True)
        source_path = check_if_test_source(source_path)

        header = Header()
        image, _ = header.get_image_from_path_2(source_path)

        source = SourceModel()

        org, mdl, anl, mdl2, anl2, anlDerivative,residuals, _= source.process(image,header)

        plot_dict = get_plot_dictionaries(header, org, mdl2, anl2, anlDerivative,residuals)
        
        plot_results(plot_dict)

    elif dir_path:
        set_logger(log_path='QuasarModelLog', console=True)
        model_several_sources(dir_path)

    else:
        root = tk.Tk()
        GUI(root)
        root.mainloop()

if __name__ == "__main__":
    typer.run(main)
