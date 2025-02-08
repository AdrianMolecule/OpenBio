from copy import deepcopy
import math
from multiprocessing import log_to_stderr
from tkinter import Canvas, Button
import tkinter as tk
from tkinter import filedialog, messagebox, Menu
#
import sys
from pathlib import Path
from typing import no_type_check  
from PIL import ImageGrab
#
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import time
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location
# mine
import gl
from myseqrecord import MySeqRecord
from model import Model
from enhancedbutton import EnhancedButton
from primers import PrimerUtils
from preferences import Preferences

def drawCanvas( )->int:
	yPos:int = 0  # y is 0 at top and increases downwards	
	# Clear any previous drawings
	gl.canvas.delete("all")
	sequenceWidth=0
	sequenceIndex=0
	yFinal=0
	for i, sequenceRecord in enumerate(Model.modelInstance.sequenceRecordList):
		if not sequenceRecord.isPrimer:
			newSequenceWidth,yFinal=drawStrand( sequenceRecord,yPos)
		else:# primer
				newSequenceWidth,yFinal=drawPrimer( sequenceRecord,yPos)
				lastY=yPos
				if i>0 and Model.modelInstance.sequenceRecordList[i-1].isPrimer and  Model.modelInstance.sequenceRecordList[i-1].hybridizedToStrand and sequenceRecord.hybridizedToStrand and Model.modelInstance.sequenceRecordList[i-1].hybridizedToStrand.uniqueId==sequenceRecord.hybridizedToStrand.uniqueId:
					yFinal=yPos# do not advance if 2 primers are on the same strand
		if newSequenceWidth>sequenceWidth:
			sequenceWidth=newSequenceWidth
		yPos=yFinal
		sequenceIndex+=1
	drawRulerAndMask(yFinal)
	# Adjust the scrollable region based on the length of the string
	gl.canvas.config(scrollregion=(0, 0, sequenceWidth, yFinal+40))  # Set the scrollable area
	return 2*gl.canvasLeftPadding+sequenceWidth

#draw features from original or cached features
def drawFeatures( mySequenceRecord: MySeqRecord, yStart: int, font: tuple[str, int], bgColor):
	for feature in mySequenceRecord.features:
		if feature.type=="XXXX":
			continue
		if ((feature.location.strand == 1 and mySequenceRecord.fiveTo3) or (feature.location.strand == -1 and not mySequenceRecord.fiveTo3) or feature.location.strand==None):
			loc: Location = feature.location
			if gl.shrink:
				x=gl.canvasLeftPadding + gl.baseRectangleSymbolXPixelSize * (loc.start+mySequenceRecord.xStartOffsetAsLetters+gl.maskSkipped[loc.start+mySequenceRecord.xStartOffsetAsLetters])
			else:
				x=gl.canvasLeftPadding + gl.baseRectangleSymbolXPixelSize * (loc.start+mySequenceRecord.xStartOffsetAsLetters)
			if getLabel(feature):
				text=getLabel(feature)[0]
			else:	
				text=feature.type
			drawTextInRectangle(text, gl.canvas, 	x, yStart, 
					   gl.baseRectangleSymbolXPixelSize, gl.baseRectangleSymbolYPixelSize, 
					bgColor, loc.end - loc.start, font)

def drawPrimer(myRecordPrimer:MySeqRecord, yStart:int)->int:	
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	verticalSequenceSpacing:int=gl.verticalSequenceSpacing+5 # blank space between 2 sequences Adrian
	upsideDownLetter:bool=not myRecordPrimer.fiveTo3 and gl.rotated
	dnaSequenceStr=MySeqRecord.seqToString(myRecordPrimer.seq)	
	featureYStart, sequenceYStart, bandYEnd = calculateYs( myRecordPrimer,  yStart)		
	i=0
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]		
		if gl.shrink and myRecordPrimer.hybridizedToStrand: # non attached primers do not have a clear x position as they float in the liquid
			xLett=gl.canvasLeftPadding + (myRecordPrimer.xStartOffsetAsLetters+i+gl.maskSkipped[myRecordPrimer.xStartOffsetAsLetters+i])*gl.baseRectangleSymbolXPixelSize
		else:
			xLett=gl.canvasLeftPadding + (myRecordPrimer.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize
		color="white" if not gl.coloredBases else gl.prefs.getPreferenceValue(letter)
		drawBase(letter, xLett, sequenceYStart,  color=color, font=font, notAnealedLocation=myRecordPrimer.notAnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter,fiveTo3=myRecordPrimer.fiveTo3)
	drawFeatures( myRecordPrimer, featureYStart,  font, 'pink')
	if gl.shrink:
		maxX=gl.canvasLeftPadding + (len(dnaSequenceStr)-(gl.maskSkipped[i+myRecordPrimer.xStartOffsetAsLetters])*gl.baseRectangleSymbolXPixelSize)
	else:	
		maxX=gl.canvasLeftPadding + (len(dnaSequenceStr)-(gl.maskSkipped[i])*gl.baseRectangleSymbolXPixelSize)
		
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton(myRecordPrimer.description[:8], 1, yStart,myRecordPrimer,  labelHeightPx=bandYEnd-yStart)	
	if  gl.hydrogen : 
		bandYEnd+=2*gl.hydrogenLinesHalfLength# add a gap to draw the hydrogen bond lines		
	return  maxX,bandYEnd


def drawStrand(mySequenceRecord:MySeqRecord, yStart:int)->int:
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	# gl.prefs.dump()
	featureYStart, sequenceYStart, bandYEnd = calculateYs( mySequenceRecord,  yStart)
	seq:Seq=mySequenceRecord.seq
	dnaSequenceStr=MySeqRecord.seqToString(seq)
	spamCount=0
	upsideDownLetter:bool=not mySequenceRecord.fiveTo3 and gl.rotated
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]			
		if gl.shrink:
			xLett=gl.canvasLeftPadding + (mySequenceRecord.xStartOffsetAsLetters+i+gl.maskSkipped[mySequenceRecord.xStartOffsetAsLetters+i])*gl.baseRectangleSymbolXPixelSize	
			excitingLetter: bool= gl.mask[i+mySequenceRecord.xStartOffsetAsLetters]
			if excitingLetter==True:
				spamCount=0
				color =gl.prefs.getPreferenceValue(letter)
			else:
				color="grey"
				spamCount+=1
			if spamCount<=3:# keep 3 letters
				drawBase(letter,  xLett, sequenceYStart, color=color, font=font,notAnealedLocation=mySequenceRecord.notAnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter, fiveTo3=mySequenceRecord.fiveTo3)
		else: # NO SHRINK
			xLett=gl.canvasLeftPadding + (mySequenceRecord.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize	
			if not gl.coloredBases:
				color="white"
			else:
				color =gl.prefs.getPreferenceValue(letter)
			drawBase(letter,  xLett, sequenceYStart, color=color, font=font, notAnealedLocation=mySequenceRecord.notAnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter, fiveTo3=mySequenceRecord.fiveTo3)

	# Replace the placeholder with the call to the new method
	drawFeatures(mySequenceRecord, featureYStart,  font,  "white")
	
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton( mySequenceRecord.description[:8], 1, yStart,mySequenceRecord, labelHeightPx=bandYEnd-yStart)	
	if mySequenceRecord.fiveTo3:
		if  gl.hydrogen : 
			bandYEnd+=2*gl.hydrogenLinesHalfLength# add a gap to draw the hydrogen bond lines	
	return  xLett+gl.baseRectangleSymbolXPixelSize,bandYEnd

#calculate the yop relative Ys for features and primers and the final y
def calculateYs(mySequenceRecord, yStart): # adrian for some reason somebody added 5 to verticalSequenceSpacing
	if mySequenceRecord.fiveTo3:
		featureYStart=yStart+gl.verticalSequenceSpacing
		sequenceYStart=yStart+gl.verticalSequenceSpacing+gl.baseRectangleSymbolYPixelSize
		if mySequenceRecord.hybridizedToStrand or mySequenceRecord.hybridizedToPrimer:
			bandEnd=yStart+gl.verticalSequenceSpacing+2*gl.baseRectangleSymbolYPixelSize
		else:
			bandEnd=yStart+gl.verticalSequenceSpacing+2*gl.baseRectangleSymbolYPixelSize+gl.verticalSequenceSpacing
	else: #  3 to 5
		if mySequenceRecord.hybridizedToStrand or mySequenceRecord.hybridizedToPrimer:
			sequenceYStart=yStart
			featureYStart=yStart+gl.baseRectangleSymbolYPixelSize
			bandEnd=yStart+gl.verticalSequenceSpacing+2*gl.baseRectangleSymbolYPixelSize
		else:
			sequenceYStart=yStart+gl.verticalSequenceSpacing
			featureYStart=yStart+gl.verticalSequenceSpacing+gl.baseRectangleSymbolYPixelSize			
			bandEnd=yStart+gl.verticalSequenceSpacing+2*gl.baseRectangleSymbolYPixelSize+gl.verticalSequenceSpacing
	return featureYStart,sequenceYStart,bandEnd

def drawTextInRectangle(tex:str,canvas:Canvas, xLeft, yTop, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, charsLength,font,rotated=None):
	canvas.create_rectangle( xLeft, yTop, xLeft + charsLength*baseRectangleSymbolXPixelSize ,yTop + baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if rotated:
		# The default anchor for text in Tkinter's Canvas.create_text() method is "center". T
		canvas.create_text(xLeft+ charsLength*baseRectangleSymbolXPixelSize/2, yTop+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(xLeft+ charsLength*baseRectangleSymbolXPixelSize/2, yTop+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black")	
	# upside_down_text = tk.Text(canvas, width=gl.fontSize, height=1 , wrap=tk.WORD, bd=1, relief="solid", highlightbackground="red",font, padx=0, pady=-3) # height is the number of lines
	# upside_down_text.insert(tk.END, tex)  # Insert the 
	# textpushDown=baseRectangleSymbolYPixelSize+10
	# canvas.create_window(x+length*baseRectangleSymbolXPixelSize/2, textpushDown+y, width=length*baseRectangleSymbolXPixelSize+1, height=baseRectangleSymbolYPixelSize+3,  window=upside_down_text)
	
# Define functions to draw each DNA base x,y relative from upper 	left corner
def drawBase(base:str, x, y, color, font, notAnealedLocation:tuple[int,int]=None, letterIndex:int=0,upsideDownLetter=None, fiveTo3=None):
	gl.canvas.create_rectangle( x, y, x + gl.baseRectangleSymbolXPixelSize ,y + gl.baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if upsideDownLetter:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)			
	else:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")
	if gl.hydrogen and (notAnealedLocation==None or (notAnealedLocation[0]<=letterIndex<notAnealedLocation[1])):
		drawHydrogenBondLines( base, x+1,y + gl.baseRectangleSymbolYPixelSize+1, fiveTo3, color)

	gl.canvas.tag_bind(textId, "<Button-1>", lambda event: clickOnSeqRecord(event, gl.canvas, None))

def drawHydrogenBondLines( base:str, x: int, topY: int, fiveTo3, color)->int:#y is always top
	lineHalfLength = gl.hydrogenLinesHalfLength-1  # Length of each vertical line
	lineSpacing = 2  # Spacing between the lines
	if fiveTo3:
		for i in range(3 if base=="C" or base=="G" else 2):
			gl.canvas.create_line(x+3, topY, x+3, topY + lineHalfLength, fill="black") #color
			x += lineSpacing	
	else:
		for i in range(3 if base=="C" or base=="G" else 2):
			gl.canvas.create_line(x+3, topY-gl.baseRectangleSymbolYPixelSize-2, x+3, topY-gl.baseRectangleSymbolYPixelSize -2- lineHalfLength-1, fill="black")
			x += lineSpacing			
	return lineHalfLength

# # Define functions to draw each DNA base x,y relative from upper left corner
# def drawTextInRectangleWithoutWidgets(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, col, font, length,rotated=None):
# 	canvas.create_rectangle( x, y, x + length*baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=col, outline="black" )
# 	if rotated:
# 		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=base, font=font, fill="black", angle=180)	
# 	else:
# 		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=base, font=font, fill="black")	

# # Define functions to draw each DNA base x,y relative from upper left corner
# def drawBaseWithoutWidgets(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, col, font, coloredBases, rotated=None):
# 	if coloredBases==False:
# 		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill="white", outline="black" )
# 	else:
# 		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=col, outline="black" )
# 	if rotated:
# 		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)	
# 	else:
# 		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")	

def buildMask():# cell is True is visible
	# for rec in Model.modelInstance.sequenceRecordList:
	size: int= max(((len(rec.seq)+rec.xStartOffsetAsLetters) for rec in  Model.modelInstance.sequenceRecordList), default=0)
	
	gl.mask=[0] * size
	gl.maskSkipped=[0] * size
	# skipped:int=0
	for rec in Model.modelInstance.sequenceRecordList:
		rec:MySeqRecord
		if rec.isPrimer:
			if  rec.hybridizedToStrand:
				for cell in range(rec.xStartOffsetAsLetters,rec.xStartOffsetAsLetters+len(rec.seq)):
					gl.mask[cell]=1
		else:#strand
			for feature in rec.features:
				if rec.toIgnore(feature): # skip features that span the full sequence
					continue
				for cell in range(feature.location.start+rec.xStartOffsetAsLetters,feature.location.end+rec.xStartOffsetAsLetters):
					gl.mask[cell]=1
	updateMaskSkipped()

def updateMaskSkipped():# todo call it only if skipped and the current one is None
	visibleSpamCount = 0
	skip = 0
	for i, _ in enumerate(gl.mask):
		excitingLetter: bool = gl.mask[i]
		if excitingLetter:
			visibleSpamCount = 0
		else:
			visibleSpamCount += 1
			if visibleSpamCount > 3:  # to keep 3 letters
				skip += 1
		gl.maskSkipped[i] = -skip  # gray but visible

def printCanvas():
	# filePath=str(Path(__file__).resolve().parent)+"/junk"
	# canvas.postscript(file=filePath, colormode='color')
 	# # Get the canvas's total scrollable area (scrollregion)
	scrollregion = gl.canvas.bbox("all")  # Returns (x1, y1, x2, y2)    
	if scrollregion:
		x1, y1, x2, y2 = scrollregion
		
		# Get the position of the canvas on the screen
		canvas_x = gl.canvas.winfo_rootx()
		canvas_y = gl.canvas.winfo_rooty()

		# Calculate the full area of the canvas to capture
		capture_x1 = canvas_x + x1
		capture_y1 = canvas_y + y1
		capture_x2 = canvas_x + x2
		capture_y2 = canvas_y + y2

		# Capture the entire scrollable area of the canvas
		img = ImageGrab.grab(bbox=(capture_x1, capture_y1, capture_x2, capture_y2), all_screens=True)
		
		# Ask the user for a file name and save the image
		file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
		if file_path:
			img.save(file_path)
			print(f"Canvas content saved to {file_path}")
	else:
		print("Canvas has no scrollable content.")

def saveModel():	
	currentDir=str(Path(__file__).resolve().parent)
	format:list=gl.prefs.getPreferenceValue(preference_name="format")
	formatPrompt=format.split(",")[0]
	formatExtension=format.split(",")[2]
	filePath = filedialog.asksaveasfilename(defaultextension="."+formatExtension,  
			filetypes=[(formatPrompt, "*."+format), (formatPrompt, "*."+formatExtension), ("All files", "*.*")], title="Save in "+formatPrompt+" file format")
	if filePath:
		if not filePath.endswith("."+formatExtension):
			filePath += "."+formatExtension
		formatName=gl.prefs.getPreferenceValue(preference_name="format").split(",")[1]
		if formatName !="genbank" and formatName !="embl":
			messagebox.showerror("Unrecognized format", f"Format names can be either genbank or embl. Please change preferences") 
			return 				
		record=Model.modelInstance.sequenceRecordList[0]
		try:
			with open(filePath, 'w') as file:
				SeqIO.write(record, file, formatName)
		except Exception as e:
			messagebox.showerror("Error", f"An error occurred while writing the file: {e}")	
		messagebox.showinfo("Success", f"Sequence exported as {filePath}")
	else:
		messagebox.showwarning("Warning", "No file name was provided!")		
		
def updateModel(seqRecList, filePath=None):
	if Model.modelInstance!=None:
		for newRecord in seqRecList:
			Model.modelInstance.sequenceRecordList.append(newRecord)
	else:
		newModel=Model(filePath,seqRecList)
		Model.modelInstance=newModel	

# the ID line is parsed in  C:\a\diy\pythonProjects\DNAPrinting\.venv\Lib\site-packages\Bio\GenBank\Scanner.py EmblScanner._feed_first_line and the parsing in line 788
def loadSequencesFile(filePath=None)->tuple[list[MySeqRecord], str]:	
	currentDir=str(Path(__file__).resolve().parent)
	format:list=gl.prefs.getPreferenceValue(preference_name="format")
	formatPrompt=format.split(",")[0]
	formatName=format.split(",")[1]
	formatExtension=format.split(",")[2]
	if not filePath:
		if not gl.debug:
			filePath = filedialog.askopenfilename(title="Open "+formatPrompt,filetypes=[(formatPrompt, "*."+formatExtension), ("All Files", "*.*")])  
		else:
			filePath = filedialog.askopenfilename(title="Open "+formatPrompt, initialdir=currentDir+"/samples/", filetypes=[(formatPrompt, "*."+formatExtension), ("All Files", "*.*")])  	
	seqRecList:list[MySeqRecord]=None			
	if filePath:
		try:
			seqRecList=loadAndSeparateSequences(filePath, formatName)
		except Exception as e:
			messagebox.showerror("Error", f"An error occurred while reading the file: {e}")
	else:
			messagebox.showwarning("No file", "Please select a file")
	if seqRecList is None or len(seqRecList)==0:
		messagebox.showerror("No Sequences", f" Please select a file that has at least one sequence") 
		# raise FileNotFoundError("no model loaded")
	return seqRecList, filePath	


def loadAndSeparateSequences(filePath:str, formatName:str)->list[MySeqRecord]:
	sequenceRecordIterator=SeqIO.parse(filePath, format=formatName)
	sequenceRecordList:list[MySeqRecord]=[]
	for mySeqRecord in sequenceRecordIterator:
		mySeqRecord:MySeqRecord
		if mySeqRecord.features==None:
			mySeqRecord.features=list()
		singleStranded=True if mySeqRecord.annotations.get("molecule_type")=="ss-DNA" else False
		if singleStranded:
			myRecord=MySeqRecord(mySeqRecord,True, True, primer=False)
			myRecord.removeFullSpanningFeatures()
			myRecord.singleStranded=True
			sequenceRecordList.append(myRecord)	
		else:# create 2 ss strands
			myRecord53=MySeqRecord(mySeqRecord,False, True, primer=False)
			myRecord53.singleStranded=True
			myRecord53.removeFullSpanningFeatures()
			sequenceRecordList.append(myRecord53)
			threeTo5Record=deepcopy(myRecord53)	
			threeTo5Record.seq=threeTo5Record.seq.complement()
			myRecord35=MySeqRecord(threeTo5Record,False, False, primer=False)
			myRecord35.removeFullSpanningFeatures()
			sequenceRecordList.append(myRecord35)		
			myRecord53.hybridizedToStrand=myRecord35
			myRecord35.hybridizedToStrand=myRecord53
	return 	sequenceRecordList 


def canvasZoom(zoomin):# 1 for  zoom In or bigger
	if zoomin!=True and gl.prefs.getPreferenceValue("fontSize")>3: # ZOOM OUT no to negative font sizes
		gl.prefs.setFontSize(gl.fontSize-1)    
		refresh()
	else:   	
		gl.prefs.setFontSize(gl.fontSize+1 )          
	refresh()

def printRed(message:str):
	print('\033[1;31m' + message + '\033[0m') 

# def checkTranslation( inSeq, outSeq):
# 	"""Checks the translation of the engineered sequence against the wild-type sequence"""
# 	myInSeq=MySeqRecord(Seq(inSeq))
# 	myOutSeq=MySeqRecord(Seq(outSeq))
# 	if myInSeq.translate().seq==myOutSeq.translate().seq:
# 		successFlag=True
# 	else:
# 		successFlag=False
# 	return successFlag

def getLabel (feature:SeqFeature):
	return feature.qualifiers.get("label")


def drawRulerAndMask( yStart)->int:
	font:tuple[str,int]=(gl.fontName,6)
	x=gl.canvasLeftPadding
	y=yStart	+gl.baseRectangleSymbolYPixelSize
	if gl.ruler:
		maxLen=0
		for i, sequenceRecord in enumerate(Model.modelInstance.sequenceRecordList):
				sequenceRecord:MySeqRecord
				if len(sequenceRecord.seq)+sequenceRecord.xStartOffsetAsLetters>maxLen:
					maxLen=len(sequenceRecord.seq)+sequenceRecord.xStartOffsetAsLetters
		for p in range(0 if gl.debug else 1,maxLen if gl.debug else maxLen-1):
			gl.canvas.create_line(x+gl.baseRectangleSymbolYPixelSize/2, y,x+gl.baseRectangleSymbolYPixelSize/2,y-(gl.baseRectangleSymbolYPixelSize), fill="grey", width=1)
			# gl.canvas.create_line(x+gl.baseRectangleSymbolYPixelSize/2, y,x+gl.baseRectangleSymbolYPixelSize/2,0, fill="grey", width=1)
			drawTextInRectangle( str(p),gl.canvas, x, y,
						gl.baseRectangleSymbolXPixelSize , gl.baseRectangleSymbolYPixelSize , 
							"white",1,font)	
			x+=	gl.baseRectangleSymbolXPixelSize			
		y+=	4
	x=gl.canvasLeftPadding	
	if gl.debug:
		for i, bit in enumerate(gl.mask):
			# canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , 
			# 				fill="Yellow" if gl.mask[i] else "grey", outline="black" )
			drawTextInRectangle( str((-gl.maskSkipped[i])),gl.canvas, x, y+gl.baseRectangleSymbolYPixelSize,
						gl.baseRectangleSymbolXPixelSize , gl.baseRectangleSymbolYPixelSize , 
							"Yellow" if gl.mask[i] else "lightgrey",1,font)	
			x+=	gl.baseRectangleSymbolXPixelSize
	return yStart+4*gl.baseRectangleSymbolYPixelSize		


def refresh():
	Preferences.updateGlobalCache()	
	buildMask()
	gl.canvasLeft.delete("all")
	drawCanvas()     

def addPrimer(filePath=None)->Seq:
	seqRecList, filePath= loadSequencesFile(filePath)
	if not seqRecList:
		return
	if len(seqRecList)>1:
		messagebox.showerror("Too many sequences", f" A primer should contain ony one sequence and this one contains {len(seqRecList)}") 
		return None
	newPrimerRecord:MySeqRecord=seqRecList[0]
	newPrimerRecord.isPrimer=True
	# newPrimerRecord.5to3 is unknown at this point because the primer is not anealed
	if not newPrimerRecord.annotations.get("molecule_type")=="ss-DNA":
		messagebox.showerror("Not a primer candidate", f" A primer should be single stranded but this record does not have molecule_type =ss-DNA") 
		return None
	if not newPrimerRecord.features or (len(newPrimerRecord.features)==1 and newPrimerRecord.toIgnore(newPrimerRecord.features[0])):
		# add a feature spanning the full length of the primer
		mandatoryFeatureText="None"
		if newPrimerRecord.description !="":
			mandatoryFeatureText=newPrimerRecord.description# this is what is shown
		else:        
			mandatoryFeatureText="AddedFeature"
		newPrimerRecord.addFeature(0, len(newPrimerRecord.seq), strand=None, type="primer", id="a primer", label=mandatoryFeatureText)	
	myRecord:MySeqRecord=MySeqRecord(newPrimerRecord, True,fiveTo3=True,primer=True)
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
	refresh() 
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def denaturate( ):
	for sequenceRecord in Model.modelInstance.sequenceRecordList:
		sequenceRecord:MySeqRecord
		if  sequenceRecord.hybridizedToStrand or sequenceRecord.hybridizedToPrimer:
			sequenceRecord.hybridizedToStrand=False
			sequenceRecord.hybridizedToPrimer=False
			sequenceRecord.singleStranded=True
			sequenceRecord.notAnealedLocation=None
	refresh() 

def anealPrimers(anealEvenIfWeCannotElongate:bool=False): # aneal only if the anealed primer can be elongated. This option is not good for determining loops as loops need to be shown even if there is no elongation possible
	found:set= set()
	added:bool=False
	for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
		sequenceRecordPrimer:MySeqRecord
		if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
			complementedPrimerSeq:Seq=sequenceRecordPrimer.seq.complement()
			complementedReversedPrimerSeq:Seq=sequenceRecordPrimer.seq.reverse_complement()
			for s, strandRegularRecord in enumerate(Model.modelInstance.sequenceRecordList):
				strandRegularRecord: MySeqRecord
				# check is not blocked by another strand or the same primer
				if not strandRegularRecord.hybridizedToStrand and not strandRegularRecord.isPrimer and (not strandRegularRecord.hybridizedToPrimer or
					(strandRegularRecord.hybridizedToPrimer and sequenceRecordPrimer.uniqueId!=strandRegularRecord.hybridizedToPrimer.uniqueId)): #and not strandRegularRecord.hybridizedToPrimer: # adrian avoid adding the same primer twice on same strand
					if strandRegularRecord.fiveTo3: # <----
						# print("Testing 5to3",strandRegularRecord.seq) 
						overlaps, largestOverlapsInStrand, largestOverlapsInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedReversedPrimerSeq)
						# if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
						# 	messagebox.showinfo("Problem", f"primer {sequenceRecordPrimer.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
						# 	return
						if largestOverlapsInStrand and len (largestOverlapsInStrand)>=1:
							if not added:
								if not anealEvenIfWeCannotElongate :
									can,where, perfectMatch= canElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer), fiveTo3Strand=strandRegularRecord.fiveTo3) # we add because the tail of the 5to3 primer coincides with the tail of the strand overlap region
								else:
									where, perfectMatch=whereToElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer))
								if anealEvenIfWeCannotElongate or can:
									if not perfectMatch:
										sequenceRecordPrimer.setNotAnealedLocation((where[0], where[1]+1))
									sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapsInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
									if not perfectMatch:
										delta=sequenceRecordPrimer.xStartOffsetAsLetters-strandRegularRecord.xStartOffsetAsLetters
										strandRegularRecord.setNotAnealedLocation((where[0]+delta, where[1]+1+delta))											
									# change feature location
									sequenceRecordPrimer.hybridizedToStrand=strandRegularRecord
									strandRegularRecord.hybridizedToPrimer=sequenceRecordPrimer 
									primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
									primRec.fiveTo3=False
									primRec.seq=Seq(MySeqRecord.seqToString(primRec.seq)[::-1])# reverses the string
									Model.modelInstance.sequenceRecordList.insert(s+1,primRec)
									added=True
							else:
								messagebox.showwarning("Warning",f"multiple anealing sites\n{found}\nAnealing to the first target that can be elongated")														
							found.add(( tuple(largestOverlapsInStrand), tuple(largestOverlapsInPrimer)))
					else:  # strand is 3 to 5
						# print("Testing 3to5",strandRegularRecord.seq)                 
						overlaps, largestOverlapsInStrand, largestOverlapsInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq)                      
						if largestOverlapsInStrand and len (largestOverlapsInStrand)>=1:
							if not added:
								if not anealEvenIfWeCannotElongate :
									can, where,perfectMatch= canElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer), fiveTo3Strand=strandRegularRecord.fiveTo3) 
								else:
									where, perfectMatch=whereToElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer))
								if anealEvenIfWeCannotElongate or can:
									if not perfectMatch:
										sequenceRecordPrimer.setNotAnealedLocation((where[0], where[1]+1))
									sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapsInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
									if not perfectMatch:
										delta=sequenceRecordPrimer.xStartOffsetAsLetters-strandRegularRecord.xStartOffsetAsLetters
										strandRegularRecord.setNotAnealedLocation((where[0]+delta, where[1]+1+delta))											
									sequenceRecordPrimer.hybridizedToStrand=strandRegularRecord
									strandRegularRecord.hybridizedToPrimer=sequenceRecordPrimer
									primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
									primRec.fiveTo3=True
									Model.modelInstance.sequenceRecordList.insert(s,primRec)    
									added=True
							else:
								messagebox.showwarning("Warning",f"multiple anealing sites\n{found}\n Anealing to first target that can elongate")
							found.add(( tuple(largestOverlapsInStrand), tuple(largestOverlapsInPrimer)))
	if added:
		refresh()                                           
	else:
		# collect names of checked primers
		names:str=""
		for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
			if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
				names+=(", "+sequenceRecordPrimer.description)
		if not found:		
			messagebox.showinfo("No Anealing", f"Primers {names} do not anneal to any present sequence") 
		else:# found match but not added
			messagebox.showinfo("No Anealing", f"Primers {names} CAN anneal to {len(found)} sequences but they were not placed on the corresponding strands since they cannot elongate") 


def workflow1():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F1CF2.gb")  
	denaturate()
	anealPrimers()
	elongate()	
	denaturate()
	# time.sleep(.5) 
	deleteSequence(uniqueId=0)
	deleteSequence(1)
	leftAlignSequence(3)
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(Model.modelInstance.sequenceRecordList[0].seq, Model.modelInstance.sequenceRecordList[0].seq.reverse_complement())
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]
	myRec.features.clear()	
	myRec.addFeature(largestOverlapsInStrand[0][0], largestOverlapsInStrand[0][1]+1,strand=None, type="misc_feature",id=None,label=f"Str_Overlap_Beg{largestOverlapsInStrand[0][0]}-{ largestOverlapsInStrand[0][1]}")
	myRec.addFeature(largestOverlapsInStrand[1][0], largestOverlapsInStrand[1][1]+1,strand=None, type="misc_feature", id=None,label=f"Str_Overlap_End{largestOverlapsInStrand[1][0]}-{largestOverlapsInStrand[1][1]}")

	# print(f"largestOverlapsInPrimer {largestOverlapsInPrimer[0]}  {largestOverlapsInPrimer[1]}")
	refresh()

def workflow():
	None

def hairpins():
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(Model.modelInstance.sequenceRecordList[0].seq, Model.modelInstance.sequenceRecordList[0].seq.reverse_complement())
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]	
	halfUpperLoopSize=math.ceil((largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]-1)/2)
	remainder=(largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]+1)%2
	newUpperRec, newLowerRec=myRec.splitRecord(halfUpperLoopSize+largestOverlapsInStrand[0][1]+1,remainder )
	# i=findSequenceIndexInModel(myRec.uniqueId)
	# Model.modelInstance.sequenceRecordList[i]=newUpperRec
	Model.modelInstance.sequenceRecordList.append(newUpperRec)
	Model.modelInstance.sequenceRecordList.append(newLowerRec)
	# delete the original record from the model but leave it in the upper strand that loops
	newUpperRec.preLoop=myRec
	deleteSequence(myRec.uniqueId)
	anealPrimers(anealEvenIfWeCannotElongate=True)
	refresh()

def loopAneal():
	anealPrimers(anealEvenIfWeCannotElongate=True)
	refresh()
		
def workflowAnealForLoopPrep():
	messagebox("Not used")
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/righttofold.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)	

	###############
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(Model.modelInstance.sequenceRecordList[0].seq, Model.modelInstance.sequenceRecordList[0].seq.reverse_complement())
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]	
	halfUpperLoopSize=math.ceil((largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]-1)/2)
	remainder=(largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]+1)%2
	newUpperRec, newLowerRec=myRec.splitRecord(halfUpperLoopSize+largestOverlapsInStrand[0][1]+1,remainder )
	# i=findSequenceIndexInModel(myRec.uniqueId)
	# Model.modelInstance.sequenceRecordList[i]=newUpperRec
	Model.modelInstance.sequenceRecordList.append(newUpperRec)
	Model.modelInstance.sequenceRecordList.append(newLowerRec)
	refresh()


def loopPrep():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/righttofold.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)		
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(Model.modelInstance.sequenceRecordList[0].seq, Model.modelInstance.sequenceRecordList[0].seq.reverse_complement())
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]	
	halfUpperLoopSize=math.floor((largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]-1)/2)
	remainder=(largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]+1)%2
	newUpperRec, newLowerRec=myRec.splitRecord(halfUpperLoopSize+largestOverlapsInStrand[0][1]+1,remainder )
	# i=findSequenceIndexInModel(myRec.uniqueId)
	# Model.modelInstance.sequenceRecordList[i]=newUpperRec
	Model.modelInstance.sequenceRecordList.append(newUpperRec)
	Model.modelInstance.sequenceRecordList.append(newLowerRec)

	# print(f"largestOverlapsInPrimer {largestOverlapsInPrimer[0]}  {largestOverlapsInPrimer[1]}")
	refresh()


#only if primer.start ==0 when strand is 5 to 3 or when primer.end==primerLen for 3 to 5 strands
def canElongate(largestOverlapsInPrimer, primerLen, fiveTo3Strand):# canElongate, which primer overlap , isPerfectOverlap
	#adrian todo this looks until it finds one
	for i in range(len(largestOverlapsInPrimer)):
		if primerLen==largestOverlapsInPrimer[i][1]-largestOverlapsInPrimer[i][0]+1:
			return True, None, True
		if  fiveTo3Strand:#works well
			if largestOverlapsInPrimer[i][0]==0:
				return True, largestOverlapsInPrimer[i], False
		else: #3 to 5 strand
			if largestOverlapsInPrimer[i][1]+1==primerLen:
				return True, largestOverlapsInPrimer[i], False
	return False, None, False


def whereToElongate(largestOverlapsInPrimer, primerLen):#  which primer overlap , isPerfectOverlap
	#adrian todo this looks until it finds one
	for i in range(len(largestOverlapsInPrimer)):
		if primerLen==largestOverlapsInPrimer[i][1]-largestOverlapsInPrimer[i][0]+1:#perfect overlap
			return  None, True
	else:
			return largestOverlapsInPrimer[i], False


def toggleShrink( ):
	p=gl.shrink
	if p:
		gl.prefs.setPreferenceValue("shrink", value=False)
	else:
		gl.prefs.setPreferenceValue("shrink", True)
	refresh()            
			  
def xxx():
	#http://rna.tbi.univie.ac.at//cgi-bin/RNAWebSuite/RNAfold.cgi
	None        

def elongate():
	found:bool=False
	for i, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
		sequenceRecordPrimer:MySeqRecord
		if sequenceRecordPrimer.isPrimer and sequenceRecordPrimer.hybridizedToStrand:    
			if sequenceRecordPrimer.fiveTo3:     
				subsequence: Seq = Seq(sequenceRecordPrimer.seq+ sequenceRecordPrimer.hybridizedToStrand.seq[len(sequenceRecordPrimer)
					+sequenceRecordPrimer.xStartOffsetAsLetters -sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters:].complement())
			else:
				#subsequence: Seq = sequenceRecordPrimer.hybridizedToStrand.seq[:sequenceRecordPrimer.xStartOffsetAsLetters+len(sequenceRecordPrimer.seq)].complement()
				subsequence: Seq = Seq(sequenceRecordPrimer.hybridizedToStrand.seq[:sequenceRecordPrimer.xStartOffsetAsLetters-
								sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters].complement()+sequenceRecordPrimer.seq)				
			newSeqRec=SeqRecord(subsequence, id=sequenceRecordPrimer.id, name=sequenceRecordPrimer.name, annotations=sequenceRecordPrimer.annotations,
								 description=f" {sequenceRecordPrimer.description}")
			newMySequenceRec = MySeqRecord(newSeqRec, singleStranded=None,fiveTo3=sequenceRecordPrimer.fiveTo3,primer=False)			
			featureLabel=f"seed primer "+sequenceRecordPrimer.description
			if sequenceRecordPrimer.fiveTo3:     
				newMySequenceRec.xStartOffsetAsLetters=sequenceRecordPrimer.xStartOffsetAsLetters  
				oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(0, len(sequenceRecordPrimer.seq), strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
			else:
				oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(sequenceRecordPrimer.xStartOffsetAsLetters, sequenceRecordPrimer.xStartOffsetAsLetters+len(sequenceRecordPrimer.seq), strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
				newMySequenceRec.xStartOffsetAsLetters=sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters     
			newMySequenceRec.hybridizedToPrimer=False
			newMySequenceRec.hybridizedToStrand=sequenceRecordPrimer.hybridizedToStrand	
			newMySequenceRec.uniqueId=sequenceRecordPrimer.uniqueId    
			newMySequenceRec.setNotAnealedLocation(sequenceRecordPrimer.notAnealedLocation)    
			newMySequenceRec.features.insert(0,oldPrimerFeature)
			Model.modelInstance.sequenceRecordList.pop(i)# remove the primer to replace it with elongated sequence
			Model.modelInstance.sequenceRecordList.insert(i, newMySequenceRec)
			newMySequenceRec.singleStranded=False
			found=True
	if not found:
		messagebox.showerror("Not found", "No primer ready to elongate") 
	refresh()   

def deleteSequence(uniqueId:int):
	for i,r in enumerate(Model.modelInstance.sequenceRecordList):
		r:MySeqRecord
		if r.uniqueId ==uniqueId:
			Model.modelInstance.sequenceRecordList.pop(i)#deletion happens here
			if r.hybridizedToPrimer:
				r.hybridizedToPrimer.hybridizedToStrand=None
			if r.hybridizedToStrand:
				r.hybridizedToStrand.hybridizedToPrimer=None                
			# print(f"the clicked sequence is found{r.uniqueId}")
			refresh()
			break

def findSequenceIndexInModel(uniqueId:int):
	for i,r in enumerate(Model.modelInstance.sequenceRecordList):
		r:MySeqRecord
		if r.uniqueId ==uniqueId:
			return i

def leftAlignSequence(uniqueId:int):
	for i,r in enumerate(Model.modelInstance.sequenceRecordList):
		r:MySeqRecord
		if r.uniqueId ==uniqueId:
			r.xStartOffsetAsLetters=0
			refresh()
			break
			  
def clickOnSeqRecordToDelete( event: tk.Event, mySeqRecord:MySeqRecord) -> None:
	# Get the coordinates of the click
	x, y = event.x, event.y
	button:EnhancedButton=event.widget
	# Find the bounding box of the text
	bbox: tuple[int, int, int, int] = gl.canvasLeft.bbox("current")
	item: tuple[int, ...] = gl.canvasLeft.find_closest(x, y)  # Find the closest item to the mouse click
	if item and gl.canvasLeft.type(item)=="text" : # Get the type of the item    
		text = gl.canvasLeft.itemcget("current", "text")
		# If the click is within the bounding box, calculate the letter clicked
		if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
			# Find the character index in the text that was clicked
			char_index = int((x - bbox[0]) / (bbox[2] - bbox[0]) * len(text))
			clicked_char = text[char_index]
			# print(f"Action 1 triggered: {text}, Clicked character: '{clicked_char}'")
	else:
			None
			# print(f"Action 1 was not over a character")
	deleteSequence(mySeqRecord.uniqueId)


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
			screenIndex=(x-gl.canvasLeftPadding)/gl.baseRectangleSymbolXPixelSize
			print(f" Clicked character: '{clicked_char}' at index {screenIndex} where should be {mySeqRecord._seq[int(screenIndex)]}")
	else:
			print(f"Action 1 was not over a character")
	# for i,r in enumerate(Model.modelInstance.sequenceRecordList):
	#     r:MySeqRecord
	#     if r.uniqueId ==mySeqRecord.uniqueId:
	#         Model.modelInstance.sequenceRecordList.pop(i)
	#         # print(f"the clicked sequence is found{r.uniqueId}")#         drawCanvas(canvas)
	#         break

# def clickOnSeqRecordToDisplayInfo( event: tk.Event, canvas:Canvas, mySeqRecord:MySeqRecord) -> None:
#     # Get the coordinates of the click
#     x, y = event.x, event.y
#     button:EnhancedButton=event.widget
#     # Find the bounding box of the text
#     bbox = canvas.bbox("current")
#     item = canvas.find_closest(x, y)  # Find the closest item to the mouse click
#     if item and canvas.type(item)=="text" : # Get the type of the item    
#         text = canvas.itemcget("current", "text")
#         # If the click is within the bounding box, calculate the letter clicked
#         if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
#             # Find the character index in the text that was clicked
#             char_index = int((x - bbox[0]) / (bbox[2] - bbox[0]) * len(text))
#             clicked_char = text[char_index]
#             print(f"Action 1 triggered: {text}, Clicked character: '{clicked_char}'")
#     else:
#             print(f"Action 1 was not over a character")
#     for i,r in enumerate(Model.modelInstance.sequenceRecordList):
#         r:MySeqRecord
#         if r.uniqueId ==mySeqRecord.uniqueId:
#             Model.modelInstance.sequenceRecordList.pop(i)
#             # print(f"the clicked sequence is found{r.uniqueId}")
#             drawCanvas(canvas)
#             break
	# import traceback
    # tb = traceback.extract_tb(e.__traceback__)
    # for filename, line, func, text in tb:
    #     print(f"Error in {func} at {filename}:{line} - {text}")