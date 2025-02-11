from copy import deepcopy
import math
from multiprocessing import log_to_stderr
from tkinter import Canvas, Button
import tkinter as tk
from tkinter import filedialog, messagebox, Menu
#
import sys
from pathlib import Path
from turtle import left
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
	gl.canvas.config(scrollregion=(0, 0, sequenceWidth+50, yFinal+40))  # Set the scrollable area
	return 2*gl.canvasLeftPadding+sequenceWidth

#draw features from original or cached features
def drawFeatures( myRec: MySeqRecord, yStart: int, font: tuple[str, int], bgColor):
	for feature in myRec.features:
		# if ((feature.location.strand == 1 and myRec.fiveTo3) or (feature.location.strand == -1 and not myRec.fiveTo3) or feature.location.strand==None):
		loc: Location = feature.location
		if gl.shrink and not (myRec.isPrimer and not myRec.hybridizedToStrand):# don't shrink unAnnealed primers
			x=gl.canvasLeftPadding + gl.baseRectangleSymbolXPixelSize * (loc.start+myRec.xStartOffsetAsLetters+gl.maskSkipped[loc.start+myRec.xStartOffsetAsLetters])
		else:
			x=gl.canvasLeftPadding + gl.baseRectangleSymbolXPixelSize * (loc.start+myRec.xStartOffsetAsLetters)
		if getLabel(feature):
			text=getLabel(feature)[0]
		else:	
			text=feature.type
		drawTextInRectangle(text, gl.canvas, 	x, yStart, 
					gl.baseRectangleSymbolXPixelSize, gl.baseRectangleSymbolYPixelSize, 
				bgColor, loc.end - loc.start, font)

def drawPrimer(myRec:MySeqRecord, yStart:int)->int:	
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	verticalSequenceSpacing:int=gl.verticalSequenceSpacing+5 # blank space between 2 sequences Adrian
	upsideDownLetter:bool=not myRec.fiveTo3 and gl.rotated
	dnaSequenceStr=MySeqRecord.seqToString(myRec.seq)	
	featureYStart, sequenceYStart, bandYEnd = calculateYs( myRec,  yStart)		
	i=0
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]		
		if gl.shrink and myRec.hybridizedToStrand: # non attached primers do not have a clear x position as they float in the liquid
			xLett=gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters+i+gl.maskSkipped[myRec.xStartOffsetAsLetters+i])*gl.baseRectangleSymbolXPixelSize
		else:
			xLett=gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize
		color="white" if not gl.coloredBases else gl.prefs.getPreferenceValue(letter)
		drawBase(letter, xLett, sequenceYStart,  color=color, font=font, notAnnealedLocation=myRec.notAnnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter,fiveTo3=myRec.fiveTo3)
	drawFeatures( myRec, featureYStart,  font, 'pink')
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton(myRec.description[:5], 1, yStart,myRec,  labelHeightPx=bandYEnd-yStart)	
	if myRec.loopInfo:# draw the loop
		myRec.visualModel=(gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters)*gl.baseRectangleSymbolXPixelSize, sequenceYStart, xLett+gl.baseRectangleSymbolXPixelSize, bandYEnd)
		drawLoop(myRec,bandYEnd)
	if  gl.hydrogen : 
		bandYEnd+=2*gl.hydrogenLinesHalfLength# add a gap to draw the hydrogen bond lines		
	return  xLett+gl.baseRectangleSymbolXPixelSize,bandYEnd

def drawStrand(myRec:MySeqRecord, yStart:int)->int:
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	# gl.prefs.dump()
	featureYStart, sequenceYStart, bandYEnd = calculateYs( myRec,  yStart)
	dnaSequenceStr=MySeqRecord.seqToString(myRec.seq)
	spamCount=0
	upsideDownLetter:bool=not myRec.fiveTo3 and gl.rotated
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]			
		if gl.shrink:
			xLett=gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters+i+gl.maskSkipped[myRec.xStartOffsetAsLetters+i])*gl.baseRectangleSymbolXPixelSize	
			excitingLetter: bool= gl.mask[i+myRec.xStartOffsetAsLetters]
			if excitingLetter==True:
				spamCount=0
				color =gl.prefs.getPreferenceValue(letter)
			else:
				color="grey"
				spamCount+=1
			if spamCount<=3:# keep 3 letters
				drawBase(letter,  xLett, sequenceYStart, color=color, font=font,notAnnealedLocation=myRec.notAnnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter, fiveTo3=myRec.fiveTo3)
		else: # NO SHRINK
			xLett=gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize	
			if not gl.coloredBases:
				color="white"
			else:
				color =gl.prefs.getPreferenceValue(letter)
			drawBase(letter,  xLett, sequenceYStart, color=color, font=font, notAnnealedLocation=myRec.notAnnealedLocation, letterIndex=i, upsideDownLetter=upsideDownLetter, fiveTo3=myRec.fiveTo3)
	# Replace the placeholder with the call to the new method
	drawFeatures(myRec, featureYStart,  font,  "white")		
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton( myRec.description[:5], 1, yStart,myRec, labelHeightPx=bandYEnd-yStart)	
	if myRec.loopInfo:# draw the loop
		myRec.visualModel=(gl.canvasLeftPadding + (myRec.xStartOffsetAsLetters)*gl.baseRectangleSymbolXPixelSize, sequenceYStart, xLett+gl.baseRectangleSymbolXPixelSize, bandYEnd)
		drawLoop(myRec,bandYEnd)
	if myRec.fiveTo3:
		if  gl.hydrogen : 
			bandYEnd+=2*gl.hydrogenLinesHalfLength# add a gap to draw the hydrogen bond lines	
	return  xLett+gl.baseRectangleSymbolXPixelSize,bandYEnd

def drawLoop(myRec:MySeqRecord,bandYEnd):
	color="black"
	# gl.canvas.create_rectangle( myRec.visualModel[0], myRec.visualModel[1], myRec.visualModel[2] ,myRec.visualModel[3] ,  outline="red", width=1 )
	lineL=gl.baseRectangleSymbolXPixelSize-2# might need to check with canvasLeftBorder
	# originalUnsplitRecord, coonectedRecord,isLeft, isTop
	if myRec.loopInfo[2]:#left one
		if not myRec.loopInfo[3]:#bottom 
			#left bottom
			# horiz
			gl.canvas.create_line(myRec.visualModel[0], myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
									myRec.loopInfo[1].visualModel[0]-lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2 , fill=color)    
			# vert
			gl.canvas.create_line(myRec.loopInfo[1].visualModel[0]-lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
								  myRec.loopInfo[1].visualModel[0]-lineL, myRec.visualModel[1]-(gl.baseRectangleSymbolYPixelSize+gl.hydrogenLinesHalfLength), fill=color) 
		else:#left top
			gl.canvas.create_line(myRec.visualModel[0], myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
								myRec.visualModel[0]	-lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2 , fill=color)    
	else:# right one
		if not myRec.loopInfo[3]:#bottom right one
			gl.canvas.create_line(myRec.visualModel[2], myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
								myRec.loopInfo[1].visualModel[2]+lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2 , fill=color) 	
			gl.canvas.create_line(myRec.loopInfo[1].visualModel[2]+lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
								  myRec.loopInfo[1].visualModel[2]+lineL, myRec.visualModel[1]-(gl.baseRectangleSymbolYPixelSize+gl.hydrogenLinesHalfLength), fill=color) 			
		else: #right top
			gl.canvas.create_line(myRec.visualModel[2], myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2,
								myRec.visualModel[2]+lineL, myRec.visualModel[1]+gl.baseRectangleSymbolYPixelSize/2 , fill=color)		

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
def drawBase(base:str, x, y, color, font, notAnnealedLocation:tuple[int,int]=None, letterIndex:int=0,upsideDownLetter=None, fiveTo3=None):
	gl.canvas.create_rectangle( x, y, x + gl.baseRectangleSymbolXPixelSize ,y + gl.baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if upsideDownLetter:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)			
	else:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")
	if gl.hydrogen and (notAnnealedLocation==None or (letterIndex<notAnnealedLocation[0] or letterIndex>notAnnealedLocation[1])):
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
			return None,None
	else:
			messagebox.showwarning("No file", "Please select a file")
	if seqRecList is None or len(seqRecList)==0:
		messagebox.showerror("No Sequences", f" Please select a file that has at least one sequence") 
		return None,None
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
			myRecord=MySeqRecord(mySeqRecord,singleStranded=True, fiveTo3=True, primer=False)
			myRecord.removeFullSpanningFeatures()
			myRecord.removeFeaturesNotOnThisStrand()
			sequenceRecordList.append(myRecord)	
		else:# create 2 ss strands
			myRecord53=MySeqRecord(mySeqRecord,True, True, primer=False)
			myRecord53.singleStranded=True
			myRecord53.removeFullSpanningFeatures()
			sequenceRecordList.append(myRecord53)
			threeTo5Record=deepcopy(myRecord53)	
			myRecord53.removeFeaturesNotOnThisStrand()
			threeTo5Record.seq=threeTo5Record.seq.complement()
			myRecord35=MySeqRecord(threeTo5Record,True, False, primer=False)
			myRecord35.removeFeaturesNotOnThisStrand()
			sequenceRecordList.append(myRecord35)		
			myRecord53.hybridizedToStrand=myRecord35
			myRecord53.hybridizedToPrimer=None
			myRecord35.hybridizedToStrand=myRecord53
			myRecord35.hybridizedToPrimer=None
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
	# newPrimerRecord.5to3 is unknown at this point because the primer is not annealed
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
	if not Model.modelInstance:
		updateModel(seqRecList, filePath=None)
	else:
		Model.modelInstance.sequenceRecordList.append(myRecord)
	# Model.modelInstance.dumpModel("in main")
	refresh() 
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def annealPrimers(annealForLoops:bool=False): # anneal only if the annealed primer can be elongated. This option is not good for determining loops as loops need to be shown even if there is no elongation possible
	found:list= list()
	added:bool=False
	for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
		sequenceRecordPrimer:MySeqRecord
		if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
			complementedPrimerSeq:Seq=sequenceRecordPrimer.seq.complement()
			complementedReversedPrimerSeq:Seq=sequenceRecordPrimer.seq.reverse_complement()
			for s, strandRegularRecord in enumerate(Model.modelInstance.sequenceRecordList):#[::-1]
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
								if not annealForLoops :
									can,where, perfectMatch= canElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer), fiveTo3Strand=strandRegularRecord.fiveTo3) # we add because the tail of the 5to3 primer coincides with the tail of the strand overlap region
								else:
									where, perfectMatch=whereToElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer))
								if annealForLoops or can:
									if not perfectMatch:
										notAnnealedLoc=(0, where[0]-1) if where[0]>0 else (where[1]+1,len(sequenceRecordPrimer)-1)	#potential bad code for partial annealing hydrogen link drawing TODO 
										sequenceRecordPrimer.setNotAnnealedLocation(notAnnealedLoc)
									sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapsInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
									if not perfectMatch:
										delta=sequenceRecordPrimer.xStartOffsetAsLetters-strandRegularRecord.xStartOffsetAsLetters
										strandRegularRecord.setNotAnnealedLocation((notAnnealedLoc[0]+delta, notAnnealedLoc[1]+delta))																				
									# change feature location
									sequenceRecordPrimer.hybridizedToStrand=strandRegularRecord
									sequenceRecordPrimer.hybridizedToPrimer=None
									strandRegularRecord.hybridizedToPrimer=sequenceRecordPrimer 
									strandRegularRecord.hybridizedToStrand=None 
									primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)	# todo understantd if popping messes up the list loop
									primRec.fiveTo3=False
									primRec.seq=Seq(MySeqRecord.seqToString(primRec.seq)[::-1])# reverses the string
									Model.modelInstance.sequenceRecordList.insert(s+1,primRec)
									found.append((largestOverlapsInStrand, largestOverlapsInPrimer))
									added=True
							else: # already added
								if annealForLoops:
									None
								else:
									found.append((largestOverlapsInStrand, largestOverlapsInPrimer))
									messagebox.showwarning("Warning",f"multiple annealing sites\n{found}\nAnnealing to the first target that can be elongated")														
					else:  # strand is 3 to 5
						# print("Testing 3to5",strandRegularRecord.seq)                 
						overlaps, largestOverlapsInStrand, largestOverlapsInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq)                      
						if largestOverlapsInStrand and len (largestOverlapsInStrand)>=1:
							if not added:
								if not annealForLoops :
									can, where,perfectMatch= canElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer), fiveTo3Strand=strandRegularRecord.fiveTo3) 
								else:
									where, perfectMatch=whereToElongate(largestOverlapsInPrimer,len(sequenceRecordPrimer))
								if annealForLoops or can:
									if not perfectMatch:
										notAnnealedLoc=(0, where[0]-1) if where[0]>0 else (where[1]+1,len(sequenceRecordPrimer)-1)
										# sequenceRecordPrimer.setNotAnnealedLocation((where[0], where[1]+1))
										sequenceRecordPrimer.setNotAnnealedLocation(notAnnealedLoc)
									sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapsInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
									if not perfectMatch:
										delta=sequenceRecordPrimer.xStartOffsetAsLetters-strandRegularRecord.xStartOffsetAsLetters
										strandRegularRecord.setNotAnnealedLocation((notAnnealedLoc[0]+delta, notAnnealedLoc[1]+delta))	
									sequenceRecordPrimer.hybridizedToStrand=strandRegularRecord
									sequenceRecordPrimer.hybridizedToPrimer=None
									strandRegularRecord.hybridizedToPrimer=sequenceRecordPrimer
									strandRegularRecord.hybridizedToStrand=None
									primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)	# todo understantd if popping messes up the list loop
									primRec.fiveTo3=True
									Model.modelInstance.sequenceRecordList.insert(s,primRec)
									found.append((largestOverlapsInStrand, largestOverlapsInPrimer))    
									added=True
							else: # already added
								found.append((largestOverlapsInStrand, largestOverlapsInPrimer))
								messagebox.showwarning("Warning",f"multiple annealing sites\n{found}\n Annealead to first target that could elongate")
	if added:
		refresh()                                           
	else:
		# collect names of checked primers
		names:str=""
		for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
			if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
				names+=(", "+sequenceRecordPrimer.description)
		if not found:		
			messagebox.showinfo("No Annealing", f"Primers {names} do not annneal to any present sequence") 
		else:# found match but not added
			messagebox.showinfo("No Annealing", f"Primers {names} CAN annneal to {len(found)} sequences but they were not placed on the corresponding strands since they cannot elongate") 

def workflow():
	None

def hairpins(anneal=True):# anneal False is just for debug purpose only
	if not Model.modelInstance.sequenceRecordList or  len(Model.modelInstance.sequenceRecordList)!=1:
		if len(Model.modelInstance.sequenceRecordList)!=2 or Model.modelInstance.sequenceRecordList[0].loopInfo==None or Model.modelInstance.sequenceRecordList[1].loopInfo==None or Model.modelInstance.sequenceRecordList[0].loopInfo[1].uniqueId !=Model.modelInstance.sequenceRecordList[1].uniqueId:
			messagebox.showerror("More than one sequence or 2 unrelated sequences", f"Hairpin analisys is done on only one unannealed sequence and we have { len(Model.modelInstance.sequenceRecordList)} sequences") 
			return
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]	
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(myRec.seq, myRec.seq.reverse_complement())
	if not largestOverlapsInStrand:
		messagebox.showerror("No annealing", f"Cannot find hairpins") 
		return
	splitPointIndex=(largestOverlapsInStrand[1][0]-largestOverlapsInStrand[0][1]-1)/2+largestOverlapsInStrand[0][1]+1
	if myRec.isLeft(splitPointIndex):
		splitPointIndex=math.floor(splitPointIndex)
	else:
		splitPointIndex=math.ceil(splitPointIndex)
	newUpperSideRec, newLowerSideRec=myRec.splitRecord(splitPointIndex=splitPointIndex )
	newUpperSideRec: MySeqRecord
	newLowerSideRec: MySeqRecord
	# i=findSequenceIndexInModel(myRec.uniqueId)	
	Model.modelInstance.sequenceRecordList.append(newUpperSideRec)
		# Model.modelInstance.sequenceRecordList[i]=newUpperRec
	Model.modelInstance.sequenceRecordList.append(newLowerSideRec)
	newUpperSideRec.loopInfo=(myRec,newLowerSideRec, myRec.isLeft(splitPointIndex),True, splitPointIndex)
	newLowerSideRec.loopInfo=(myRec,newUpperSideRec, myRec.isLeft(splitPointIndex), False, splitPointIndex)
	deleteSequenceFromModel(myRec.uniqueId)
	if( anneal):
		annealPrimers(annealForLoops=True)
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

def elongate():#todo elongate when going left should change the noAnealed region index by the number of letters elongated
	found:bool=False
	for i, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
		sequenceRecordPrimer:MySeqRecord
		if sequenceRecordPrimer.isPrimer and sequenceRecordPrimer.hybridizedToStrand:    
			if sequenceRecordPrimer.fiveTo3:     
				subsequence: Seq = Seq(sequenceRecordPrimer.seq+ sequenceRecordPrimer.hybridizedToStrand.seq[len(sequenceRecordPrimer)
					+sequenceRecordPrimer.xStartOffsetAsLetters -sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters:].complement())
			else:# 3 to 5 strand
				subsequence: Seq = Seq(sequenceRecordPrimer.hybridizedToStrand.seq[:sequenceRecordPrimer.xStartOffsetAsLetters-
								sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters].complement()+sequenceRecordPrimer.seq)				
			newSeqRec=SeqRecord(subsequence, id=sequenceRecordPrimer.id, name=sequenceRecordPrimer.name, annotations=sequenceRecordPrimer.annotations,
								 description=f" {sequenceRecordPrimer.description}")
			newElongatedRecord = MySeqRecord(newSeqRec, singleStranded=None,fiveTo3=sequenceRecordPrimer.fiveTo3,primer=False)			
			featureLabel=f"seed primer "+sequenceRecordPrimer.description
			if sequenceRecordPrimer.fiveTo3:     
				newElongatedRecord.xStartOffsetAsLetters=sequenceRecordPrimer.xStartOffsetAsLetters  
				oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(0, len(sequenceRecordPrimer.seq), strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
				if sequenceRecordPrimer.notAnnealedLocation:
					newElongatedRecord.setNotAnnealedLocation(sequenceRecordPrimer.notAnnealedLocation)    
			else:# primer is 3 to 5
				oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(sequenceRecordPrimer.xStartOffsetAsLetters-
														  sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters,
														   sequenceRecordPrimer.xStartOffsetAsLetters-
														   sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters+len(sequenceRecordPrimer.seq),
															 strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
				newElongatedRecord.xStartOffsetAsLetters=sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters  
				if sequenceRecordPrimer.notAnnealedLocation:
					shift= len(newElongatedRecord.seq)-len(sequenceRecordPrimer)
					newElongatedRecord.setNotAnnealedLocation((sequenceRecordPrimer.notAnnealedLocation[0]+shift,sequenceRecordPrimer.notAnnealedLocation[1]+shift))    
				# shift the Not annealed locations  
			#change the old sequence to point to the new one where 
			sequenceRecordPrimer.hybridizedToStrand.hybridizedToPrimer=None		
			sequenceRecordPrimer.hybridizedToStrand.hybridizedToStrand=newElongatedRecord		
			newElongatedRecord.hybridizedToStrand=sequenceRecordPrimer.hybridizedToStrand	
			newElongatedRecord.hybridizedToPrimer=None	
			newElongatedRecord.uniqueId=sequenceRecordPrimer.uniqueId    
			newElongatedRecord.features.insert(0,oldPrimerFeature)
			Model.modelInstance.sequenceRecordList.pop(i)# remove the primer to replace it with elongated sequence 	# to do understantd if popping messes up the list loop
			Model.modelInstance.sequenceRecordList.insert(i, newElongatedRecord)
			newElongatedRecord.singleStranded=False
			found=True
	if not found:
		messagebox.showerror("Not found", "No primer ready to elongate") 
	refresh()   

def deleteSequenceFromModel(uniqueId:int):# only one so no indexing errors are created
	for i,myRec in enumerate(Model.modelInstance.sequenceRecordList):
		myRec:MySeqRecord
		if myRec.uniqueId ==uniqueId:
			if myRec.hybridizedToPrimer:
				myRec.hybridizedToPrimer.hybridizedToStrand=None
				myRec.hybridizedToPrimer.notAnnealedLocation=None
				myRec.hybridizedToPrimer.singleStranded=True
			if myRec.hybridizedToStrand:
				myRec.hybridizedToStrand.hybridizedToPrimer=None                
				myRec.hybridizedToStrand.notAnnealedLocation=None  
				myRec.hybridizedToStrand.singleStranded=True              
			Model.modelInstance.sequenceRecordList.pop(i)#deletion happens here
			# print(f"the clicked sequence is found{r.uniqueId}")
			refresh()
			break	

def deleteFirstSequenceFromModel():# only one so no indexing errors are created
	for i,myRec in enumerate(Model.modelInstance.sequenceRecordList):
		myRec:MySeqRecord
		if i ==0:
			if myRec.hybridizedToPrimer:
				myRec.hybridizedToPrimer.hybridizedToStrand=None
				myRec.hybridizedToPrimer.notAnnealedLocation=None
				myRec.hybridizedToPrimer.singleStranded=True
			if myRec.hybridizedToStrand:
				myRec.hybridizedToStrand.hybridizedToPrimer=None                
				myRec.hybridizedToStrand.notAnnealedLocation=None  
				myRec.hybridizedToStrand.singleStranded=True              
			Model.modelInstance.sequenceRecordList.pop(i)#deletion happens here
			refresh()
			break	

def denaturate( ):
	for sequenceRecord in Model.modelInstance.sequenceRecordList:
		sequenceRecord:MySeqRecord
		sequenceRecord.hybridizedToStrand=None
		sequenceRecord.hybridizedToPrimer=None
		sequenceRecord.notAnnealedLocation=None
		sequenceRecord.singleStranded=True
	refresh() 

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
	deleteSequenceFromModel(mySeqRecord.uniqueId)


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


def workflow1():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F1CF2.gb")  
	denaturate()
	annealPrimers()
	elongate()	
	denaturate()
	# time.sleep(.5) 
	deleteSequenceFromModel(uniqueId=0)
	deleteSequenceFromModel(1)
	leftAlignSequence(3)
	overlaps, largestOverlapsInStrand, largestOverlapsInPrimer=PrimerUtils.findPrimerOverlaps(Model.modelInstance.sequenceRecordList[0].seq, Model.modelInstance.sequenceRecordList[0].seq.reverse_complement())
	myRec: MySeqRecord=Model.modelInstance.sequenceRecordList[0]
	myRec.features.clear()	
	myRec.addFeature(largestOverlapsInStrand[0][0], largestOverlapsInStrand[0][1]+1,strand=None, type="misc_feature",id=None,label=f"Str_Overlap_Beg{largestOverlapsInStrand[0][0]}-{ largestOverlapsInStrand[0][1]}")
	myRec.addFeature(largestOverlapsInStrand[1][0], largestOverlapsInStrand[1][1]+1,strand=None, type="misc_feature", id=None,label=f"Str_Overlap_End{largestOverlapsInStrand[1][0]}-{largestOverlapsInStrand[1][1]}")

	# print(f"largestOverlapsInPrimer {largestOverlapsInPrimer[0]}  {largestOverlapsInPrimer[1]}")
	refresh()

def addDirectMenus(menuBar):
	menuBar.add_command(label="testLoadPorkDenaturate", command=testLoadPorkDenaturate)
	menuBar.add_command(label="testFullNoHairStartWithF1cF2", command=testFullNoHairStartWithF1cF2)
	menuBar.add_command(label="testFullNoHairStartWithB1cB2", command=testFullNoHairStartWithB1cB2)
	# menuBar.add_command(label="testB1B2partial", command=testB1B2partial)
	menuBar.add_command(label="pCRsample", command=pCRsample)
	menuBar.add_command(label="testLeftLoopSplit", command=testLeftLoopSplit)
	menuBar.add_command(label="testRightLoopSplit", command=testRightLoopSplit)
	menuBar.add_command(label="testLoopAnneal", command=testLoopAnneal)
	# menuBar.add_command(label="testFeatureLabelBug", command=testFeatureLabelBug)
	# menuBar.add_command(label="testF1F2all", command=testF1F2all)

def testFeatureLabelBug():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	denaturate()
	deleteFirstSequenceFromModel()
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F3.gb")  
	annealPrimers()
	elongate()	
	denaturate()
	deleteSequenceFromModel(1)
	leftAlignSequence(3)
	# #
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/B3.gb")  
	# annealPrimers()
	# elongate()	


def testFullNoHairStartWithF1cF2():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	deleteFirstSequenceFromModel()
	denaturate()
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F1CF2.gb")  
	annealPrimers()
	elongate()	
	denaturate()
	deleteSequenceFromModel(1)
	leftAlignSequence(3)
	# #
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/B1cB2.gb")  
	annealPrimers()
	elongate()	
	denaturate()
	deleteFirstSequenceFromModel()
	# time.sleep(.5) 

def testFullNoHairStartWithB1cB2():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	denaturate()
	deleteSequenceFromModel(uniqueId=1)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/B1CB2.gb")  
	annealPrimers()
	elongate()
	denaturate()

	deleteFirstSequenceFromModel()
	# leftAlignSequence(3)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F1CF2.gb")  
	annealPrimers()
	elongate()	
	denaturate()
	deleteSequenceFromModel(uniqueId=3)	
	# leftAlignSequence(3)	
	
def pCRsample():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	denaturate()
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F3.gb")  
	annealPrimers()
	elongate()	
	denaturate()	
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/B3.gb") 
	deleteFirstSequenceFromModel()
	deleteSequenceFromModel(uniqueId=1)
	annealPrimers()
	elongate()	
	denaturate()	
	deleteSequenceFromModel(uniqueId=3)


def testB1B2all():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/B1CB2.gb")  
	denaturate()
	annealPrimers()
	elongate()	
	denaturate()
	# time.sleep(.5) 
	deleteSequenceFromModel(uniqueId=0)
	deleteSequenceFromModel(1)
	leftAlignSequence(3)
	
def testF1F2all():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/porkcomplete.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	addPrimer(filePath=str(Path(__file__).resolve().parent)+"/samples/F1CF2.gb")  
	denaturate()
	annealPrimers()
	elongate()	
	denaturate()
	# time.sleep(.5) 
	deleteSequenceFromModel(uniqueId=0)
	deleteSequenceFromModel(1)
	leftAlignSequence(3)

def testLoadPorkDenaturate():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0		
	defaultTestFileValue=gl.prefs.getPreferenceValue("defaultTestFileValue")
	filePath=str(Path(__file__).resolve().parent)+defaultTestFileValue
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	updateModel(seqRecList, filePath=filePath)
	# app.root.title("OpenBio "+Model.modelInstance.loadedFileName)
	denaturate()


def testLoopAnneal():
	annealPrimers(annealForLoops=True)
	refresh()

def testLeftLoopSplit():
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/leftSplitTestCase.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	assert seqRecList and len(seqRecList)==1
	updateModel(seqRecList, filePath=filePath)	
	gl.prefs.setPreferenceValue("minPrimerOverlapLength",4                                                                                                                  )
	Preferences.updateGlobalCache()
	hairpins(False)

def testRightLoopSplit():# cccttttcgc tcaaaa
	Model.modelInstance=None
	MySeqRecord.uniqueId=0
	filePath=str(Path(__file__).resolve().parent)+"/samples/rightSplitTestCase.gb"
	seqRecList, filePath=loadSequencesFile(filePath=filePath)
	assert seqRecList and len(seqRecList)==1
	updateModel(seqRecList, filePath=filePath)	
	gl.prefs.setPreferenceValue("minPrimerOverlapLength",4)
	Preferences.updateGlobalCache()
	hairpins(False)
