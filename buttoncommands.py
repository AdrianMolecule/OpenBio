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
from primers import PrimerUtils

def addPrimerHandler(canvas:Canvas)->Seq:
    seqRecList, filePath= loadFile( )
    if len(seqRecList)>1:
        messagebox.showerror("Too many sequences", f" A primer should contain ony one sequence and this one contains {len(seqRecList)}") 
        return None
    newRecord: MySeqRecord = seqRecList[0]
    if not newRecord.annotations.get("molecule_type")=="ss-DNA":
        messagebox.showerror("Not a primer candidate", f" A primer should be single stranded but this record does not have molecule_type =ss-DNA") 
        return None
    # add a feature spanning the full length of the primer
    mandatoryFeatureText="None"
    if newRecord.description !="":
        mandatoryFeatureText=newRecord.description# this is what is shown
    else:        
        mandatoryFeatureText="AddedFeature"
    mandatoryFeature:SeqFeature=SeqFeature(SimpleLocation(0, len(newRecord.seq), strand=None), type="primer", id="a primer", qualifiers={"label": [mandatoryFeatureText]  })
    newRecord.features.append(mandatoryFeature)   
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
    myRecord.isPrimer=True
    myRecord.singleStranded=True
    myRecord.fiveTo3=True
    myRecord.hybridizedTo=None
    Model.modelInstance.sequenceRecordList.append(myRecord)
    # Model.modelInstance.dumpModel("in main")
    drawCanvas(canvas)
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def denaturate( canvas:Canvas):
    for sequenceRecord in Model.modelInstance.sequenceRecordList:
        sequenceRecord:MySeqRecord
        if sequenceRecord.hybridizedTo:
            sequenceRecord.hybridizedTo=False
    drawCanvas(canvas)

def anealPrimers(canvas:Canvas ):
    found:bool=False
    minOverlapLength:int=gl.prefs.get_preference_value("minPrimerOverlapLength")
    for p, sequenceRecord in enumerate(Model.modelInstance.sequenceRecordList):
        sequenceRecord:MySeqRecord
        if sequenceRecord.isPrimer and not sequenceRecord.hybridizedTo :
            complementedPrimerSeq:Seq=sequenceRecord.seq.complement()
            complementedReversedPrimerSeq:Seq=sequenceRecord.seq.reverse_complement()
            for s, strandRegularRecord in enumerate(Model.modelInstance.sequenceRecordList):
                strandRegularRecord: MySeqRecord
                if not strandRegularRecord.isPrimer and not strandRegularRecord.hybridizedTo:
                    if strandRegularRecord.fiveTo3: # <----
                        # print("Testing 5to3",strandRegularRecord.seq) 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedReversedPrimerSeq, minOverlapLength=minOverlapLength)
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecord.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecord.xStartOffset=largestOverlapsInStrand[0][0]
                            Model.modelInstance.sequenceRecordList[p].hybridizedTo=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedTo=Model.modelInstance.sequenceRecordList[p]
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=False
                            primRec.seq=seqToString(primRec.seq)[::-1]
                            Model.modelInstance.sequenceRecordList.insert(s+1,primRec)
                    else:  
                        # print("Testing 3to5",strandRegularRecord.seq)                 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq, minOverlapLength=minOverlapLength)                      
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecord.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecord.xStartOffset=largestOverlapsInStrand[0][0]
                            Model.modelInstance.sequenceRecordList[p].hybridizedTo=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedTo=Model.modelInstance.sequenceRecordList[p]
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=True
                            Model.modelInstance.sequenceRecordList.insert(s,primRec)    
    if found:
        gl.prefs.set_preference_value("shrink", False)
        drawCanvas(canvas)                                            
    else:
        messagebox.showinfo("No Anealing", f"The specific primer does not anneal to any present sequence") 




def elongate( ):
  messagebox.showinfo("TBD", f"elongate") 
