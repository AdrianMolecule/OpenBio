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
from util import *
from main import UiApp
import gl


def addPrimerHandler(self)->Seq:
    seqRecList, filePath= loadFile( )
    for newRecord in seqRecList:
        newRecord.annotations=None
        newRecord.features=None
        myRecord:MySeqRecord=MySeqRecord(newRecord, True,fiveTo3=True,primer=True)
        len=len(myRecord.seq) 
        minLen=gl.prefs.get_preference_value("minPrimerLength")
        maxLen=gl.prefs.get_preference_value("maxPrimerLength")
        if len< minLen:
             messagebox.showerror("Invalid Primer", f"primer length: {len} is smaller than minimal primar length {minLen}") 
             return
        if len> maxLen:
             messagebox.showerror("Invalid Primer", f"primer length: {len} is larger than maximum primar length {maxLen}") 
             return
        myRecord.isPrimer=True
        myRecord.singleStranded=True
        myRecord.fiveTo3=True
        myRecord.hybridizedTo=None
        Model.modelInstance.sequenceRecordList.append(newRecord)
    Model.modelInstance.dumpModel("in main")
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))


def denaturate( uiApp:UiApp):
  messagebox.showerror("TBD", f"denaturate") 
def anealPrimers( uiApp:UiApp):
  messagebox.showinfo("TBD", f"anealPrimers") 
def elongate( uiApp:UiApp):
  messagebox.showinfo("TBD", f"elongate") 
