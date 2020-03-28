import tkinter as tk
from tkinter import *
from tkinter import messagebox, scrolledtext
from tkinter.filedialog import askopenfilename, askdirectory
from modules.analysis import *
import sys, os

#from matplotlib.backends.backend_tkagg import (
#    FigureCanvasTkAgg, NavigationToolbar2Tk)
#from matplotlib.backend_bases import key_press_handler
#from matplotlib.figure import Figure

#import numpy as np


#####################################################################################################
#   This is going to be where we write all the functions involving windows, buttons within it, etc  #
#   Also stuff like getting the file path to the input .txt or .csv, and path to an output file     #
#####################################################################################################

window = tk.Tk()

window.title('Peptide-Analysis Tool')
window.geometry('646x80')

#trying to change the icon on the window to the one found in the modules folder, can't figure it out

# get path to input file, sequences
inp_file = Label(window, text="Path to .txt file containing sequences: ")
inp_file.grid(row=0, column=0)

inp_path_box = Entry(window, width=30)
inp_path_box.grid(row=0,column=1)

def get_inp_path():
    fname = askopenfilename()
    if inp_path_box != '':
        inp_path_box.delete(0,END)
    inp_path_box.insert(END,fname)

inp_path_button = Button(window,text='Click to select file',command=get_inp_path, width=30)
inp_path_button.grid(row=0,column=2)


# get path to output folder
out_dir = Label(window, text="Path to output directory: ")
out_dir.grid(row=1, column=0)

out_path_box = Entry(window, width=30)
out_path_box.grid(row=1,column=1)


def get_out_path():
    fname = askdirectory()
    if out_path_box != '':
        out_path_box.delete(0,END)
    out_path_box.insert(END,fname)

out_path_button = Button(window,text='Click to select file',command=get_out_path, width=30)
out_path_button.grid(row=1,column=2)

def get_inp_filename():
    filename = inp_path_box.get() 
    if not filename.endswith('.txt'):
        return "error, incorrect input filetype"
    else:
        return filename

def analyze(): 
    fname = get_inp_filename()
    outpath = out_path_box.get()
    
    try:
        check_input_file(fname)
    except:
        messagebox.showinfo("Error!!!","Input file is not correct format. See README.txt for details.")
        window.destroy()
        return None
    
    main(fname,outpath)

analysis_button = Button(window,text='Analyze',command=analyze,width=30)
analysis_button.grid(row=2,column=1)    
