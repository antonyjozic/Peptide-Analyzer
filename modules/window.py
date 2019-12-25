import tkinter as tk
from tkinter import *
from tkinter import messagebox, scrolledtext
from tkinter.filedialog import askopenfilename, askdirectory

import sys
import os


#####################################################################################################
#   This is going to be where we write all the functions involving windows, buttons within it, etc  #
#   Also stuff like getting the file path to the input .txt or .csv, and path to an output file     #
#####################################################################################################

window = tk.Tk()

window.title('Peptide-Analysis Tool')
window.geometry('800x600')

inp_file = Label(window, text="Path to .txt file containing sequences:")
inp_file.grid(row=0, column=0, sticky=W)

inp_path_box = Entry(window, width=30)
inp_path_box.grid(row=0,column=1)

def get_inp_path():
    fname = askopenfilename()
    if inp_path_box != '':
        inp_path_box.delete(0,END)
    inp_path_box.insert(END,fname)

inp_path_button = Button(window,text='Click to select file',command=get_inp_path(), width=30)
inp_path_button.grid(row=0,column=2)

#         #
# TESTING #
#         #

window.mainloop()
