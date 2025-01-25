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
from myseqrecord import MySeqRecord
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location

import gl
from model import Model
from primers import PrimerUtils
from util import *


def addPrimerHandler(canvas:Canvas)->Seq:
    seqRecList, filePath= loadFile()
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
    minLen=gl.prefs.getPreferenceValue("minPrimerOverlapLength")
    maxLen=gl.prefs.getPreferenceValue("maxPrimerLength")
    if leng< minLen:
            messagebox.showerror("Invalid Primer", f"primer length: {len} is smaller than minimal primar length {minLen}") 
            return
    if leng> maxLen:
            messagebox.showerror("Invalid Primer", f"primer length: {len} is larger than maximum primar length {maxLen}") 
            return
    myRecord.isPrimer=True
    myRecord.singleStranded=True
    myRecord.fiveTo3=True
    myRecord.hybridizedToStrand=None
    Model.modelInstance.sequenceRecordList.append(myRecord)
    # Model.modelInstance.dumpModel("in main")
    drawCanvas(canvas)
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def denaturate( canvas:Canvas):
    for sequenceRecord in Model.modelInstance.sequenceRecordList:
        sequenceRecord:MySeqRecord
        if  not sequenceRecord.isPrimer and sequenceRecord.hybridizedToStrand:
            sequenceRecord.hybridizedToStrand=False
    drawCanvas(canvas)

def anealPrimers(canvas:Canvas ):
    found:bool=False
    minOverlapLength:int=gl.prefs.getPreferenceValue("minPrimerOverlapLength")
    for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
        sequenceRecordPrimer:MySeqRecord
        if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
            complementedPrimerSeq:Seq=sequenceRecordPrimer.seq.complement()
            complementedReversedPrimerSeq:Seq=sequenceRecordPrimer.seq.reverse_complement()
            for s, strandRegularRecord in enumerate(Model.modelInstance.sequenceRecordList):
                strandRegularRecord: MySeqRecord
                if not strandRegularRecord.isPrimer and (not strandRegularRecord.hybridizedToPrimer or(strandRegularRecord.hybridizedToPrimer and sequenceRecordPrimer.uniqueId!=strandRegularRecord.hybridizedToPrimer.uniqueId)): #and not strandRegularRecord.hybridizedToPrimer: # adrian avoid adding the same primer twice on same strand
                    if strandRegularRecord.fiveTo3: # <----
                        # print("Testing 5to3",strandRegularRecord.seq) 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedReversedPrimerSeq, minOverlapLength=minOverlapLength)
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecordPrimer.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecordPrimer.xStartOffsetAsLetters=largestOverlapsInStrand[0][0]
                            # change feature location
                            Model.modelInstance.sequenceRecordList[p].hybridizedToStrand=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedToPrimer=Model.modelInstance.sequenceRecordList[p]  
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=False
                            primRec.seq=Seq(seqToString(primRec.seq)[::-1])
                            Model.modelInstance.sequenceRecordList.insert(s+1,primRec)
                    else:  
                        # print("Testing 3to5",strandRegularRecord.seq)                 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq, minOverlapLength=minOverlapLength)                      
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecordPrimer.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecordPrimer.xStartOffsetAsLetters=largestOverlapsInStrand[0][0]
                            Model.modelInstance.sequenceRecordList[p].hybridizedToStrand=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedToPrimer=Model.modelInstance.sequenceRecordList[p]
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=True
                            Model.modelInstance.sequenceRecordList.insert(s,primRec)    
    if found:
        # gl.prefs.set_preference_value("shrink", False)
        drawCanvas(canvas)                                            
    else:
        names:str=""
        for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
            if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
                 names+=(", "+sequenceRecordPrimer.description)
    
        messagebox.showinfo("No Anealing", f"Primers {names} do not anneal to any present sequence") 


def toggleShrink( canvas:Canvas):
    p=gl.prefs.getPreferenceValue(preference_name="shrink")
    if p:
          gl.prefs.setPreferenceValue("shrink", False)
    else:
          gl.prefs.setPreferenceValue("shrink", True)
    drawCanvas(canvas)           
              

def elongate( canvas:Canvas):
    p=gl.prefs.getPreferenceValue(preference_name="shrink")
    if p:
          gl.prefs.setPreferenceValue("shrink", False)
    else:
          gl.prefs.setPreferenceValue("shrink", True)
    drawCanvas(canvas)           
              
def clickOnSeqRecord( event: tk.Event, canvas:Canvas, mySeqRecord:MySeqRecord) -> None:
    # Get the coordinates of the click
    x, y = event.x, event.y
    button:EnhancedButton=event.widget
    # Find the bounding box of the text
    bbox = canvas.bbox("current")
    item = canvas.find_closest(x, y)  # Find the closest item to the mouse click
    if item and canvas.type(item)=="text" : # Get the type of the item    
        text = canvas.itemcget("current", "text")
        # If the click is within the bounding box, calculate the letter clicked
        if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
            # Find the character index in the text that was clicked
            char_index = int((x - bbox[0]) / (bbox[2] - bbox[0]) * len(text))
            clicked_char = text[char_index]
            print(f"Action 1 triggered: {text}, Clicked character: '{clicked_char}'")
    else:
            print(f"Action 1 was not over a character")
    for i,r in enumerate(Model.modelInstance.sequenceRecordList):
        r:MySeqRecord
        if r.uniqueId ==mySeqRecord.uniqueId:
            Model.modelInstance.sequenceRecordList.pop(i)
            # print(f"the clicked sequence is found{r.uniqueId}")
            drawCanvas(canvas)
            break

