"""Utility methods"""
from copy import deepcopy
from tkinter import Canvas
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import Menu
from tkinter import filedialog
#
import sys
from pathlib import Path  
from PIL import ImageGrab
#
from Bio.SeqRecord import SeqRecord
from wrappers import MySeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location

import gl
from model import Model


def addSequencesHandler(self)->Seq:
    loadModel(False, append=True)
    self.root.title("OpenBio "+Model.modelInstance.loadedFileName)
    drawCanvas(self.canvas) 

def denaturate( root):
  messagebox.showerror("TBD", f"denaturate") 
def anealPrimers( root):
  messagebox.showinfo("TBD", f"anealPrimers") 
def elongate( root):
  messagebox.showinfo("TBD", f"elongate") 
