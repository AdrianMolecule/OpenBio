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
from typing import Type  
from PIL import ImageGrab
#
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from wrappers import MySeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location

import gl
from model import Model
from util import *
from main import UiApp
from checkPrimerBinding import findPrimerOverlaps

def addPrimerHandler(canvas:Canvas)->Seq:
    seqRecList, filePath= loadFile( )
    if len(seqRecList)>1:
        messagebox.showerror("Too many sequences", f" A primer should contain ony one sequence and this one contains {len(seqRecList)}") 
        return None
    newRecord = seqRecList[0]
    if not newRecord.annotations.get("molecule_type")=="ss-DNA":
        messagebox.showerror("Not a primer candidate", f" A primer should be single stranded but this record does not have molecule_type =ss-DNA") 
        return None
    mandatoryFeature:SeqFeature=SeqFeature(SimpleLocation(0, len(newRecord.seq), strand=1), type="primer", id="Adrian")
    if newRecord.description !="":
        mandatoryFeature.qualifiers.update({'label': newRecord.description})# this is what is shown
    else:        
        mandatoryFeature.qualifiers.update({'label': 'A Primer'})# this is what is shown
    newRecord.features.append(mandatoryFeature)   
    #todo add a feature label?                     
    # feature.qualifiers.get("label")[0]

    myRecord:MySeqRecord=MySeqRecord(newRecord, True,fiveTo3=True,primer=True)
    leng=len(myRecord.seq) 
    minLen=gl.prefs.get_preference_value("minPrimerOverlapLength")
    maxLen=gl.prefs.get_preference_value("maxPrimerLength")
    if leng< minLen:
            messagebox.showerror("Invalid Primer", f"primer length: {len} is smaller than minimal primar length {minLen}") 
            return
    if leng> maxLen:
            messagebox.showerror("Invalid Primer", f"primer length: {len} is larger than maximum primar length {maxLen}") 
            return
    myRecord.primer=True
    myRecord.singleStranded=True
    myRecord.fiveTo3=True
    myRecord.hybridizedTo=None
    Model.modelInstance.sequenceRecordList.append(myRecord)
    Model.modelInstance.dumpModel("in main")
    drawCanvas(canvas)
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def denaturate( uiApp:UiApp):
  messagebox.showerror("TBD", f"denaturate") 

def anealPrimers( uiApp:UiApp):
    minOverlapLength:int=gl.prefs.get_preference_value("minPrimerOverlapLength")
    for sequenceRecord in Model.modelInstance.sequenceRecordList:
        sequenceRecord:MySeqRecord

        if sequenceRecord.primer:
            complementedPrimerSeq:Seq=sequenceRecord.seq.complement()
            for strandRegularRecord in Model.modelInstance.sequenceRecordList:
                strandRegularRecord: MySeqRecord
                if not strandRegularRecord.primer:      
                    if not strandRegularRecord.fiveTo3: # <----
                        print("Testing 5to3",strandRegularRecord.seq) 
                        # overlaps, largestInStrand, largestInPrimer =findPrimerOverlaps(targetDnaRecord=strandRegularRecord, primerRecordSeq=complementedPrimerSeq, minOverlapLength=minOverlapLength)
                    else:  
                        print("Testing 3to5 ",strandRegularRecord.seq)                 
                        overlaps, largestOverlapInStrand, largestOverlapInPrimer =findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq, minOverlapLength=minOverlapLength)
                    if largestOverlapInStrand:                            
                        print("Overlaps:", overlaps)
                        print("Largest Overlaps in longString:", largestOverlapInStrand)
                        print("Largest Overlaps in shortString:", largestOverlapInPrimer)                        
                        print("In:", strandRegularRecord)  
        # newSequenceWidth,yFinal=drawStrand(canvas, sequenceRecord, sequenceIndex, canvasHorizontalMargin,yPos)
		# if newSequenceWidth>sequenceWidth:
		# 	sequenceWidth=newSequenceWidth
		# yPos=yFinal
		# sequenceIndex+=1


def elongate( uiApp:UiApp):
  messagebox.showinfo("TBD", f"elongate") 
