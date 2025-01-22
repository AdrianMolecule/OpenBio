"""Utility methods"""
from copy import deepcopy
from tkinter import Canvas, Button
import tkinter as tk
from tkinter import filedialog, messagebox, Menu
#
import sys
from pathlib import Path
from typing import no_type_check  
from PIL import ImageGrab
import bisect
#
from Bio.SeqRecord import SeqRecord
from wrappers import MySeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location
# from upsideDownText import UpsideDownText
# mine
import gl
from model import Model
from enhancedbutton import EnhancedButton
#######
# some 'constants' that can be changed if push comes to shove !
# some initial values that might remain constant
#canvasVerticalMargin = 30
horizontalPixelsMargin=2 # head room between the base letter and it's holding box

def drawCanvas(canvas:Canvas )->int:
	yPos:int = 0  # y is 0 at top and increases downwards	
	# Clear any previous drawings
	canvas.delete("all")
	sequenceWidth=0
	sequenceIndex=0
	for sequenceRecord in Model.modelInstance.sequenceRecordList:
		if not sequenceRecord.isPrimer:
			newSequenceWidth,yFinal=drawStrand(canvas, sequenceRecord, gl.canvasHorizontalMargin,yPos)
		else:# primer
				newSequenceWidth,yFinal=drawPrimer(canvas, sequenceRecord, gl.canvasHorizontalMargin,yPos)
		if newSequenceWidth>sequenceWidth:
			sequenceWidth=newSequenceWidth
		yPos=yFinal
		sequenceIndex+=1
	#
	# Adjust the scrollable region based on the length of the string
	canvas.config(scrollregion=(0, 0, sequenceWidth,yFinal+20))  # Set the scrollable area
	return 2*gl.canvasHorizontalMargin+sequenceWidth

# the ID line is parsed in  C:\a\diy\pythonProjects\DNAPrinting\.venv\Lib\site-packages\Bio\GenBank\Scanner.py EmblScanner._feed_first_line and the parsing in line 788
def loadFile(default=False)->tuple[list[MySeqRecord],str]:	
	if default:
		filePath=str(Path(__file__).resolve().parent)+gl.prefs.get_preference_value("defaultTestFileValue")
	else:
		filePath = filedialog.askopenfilename(title="Open EMBL File", filetypes=[("EMBL Files", "*.embl"), ("All Files", "*.*")])    
	secRecList:list[MySeqRecord]		=None
	if filePath:
		try:
			secRecList=loadAndSeparateEmblSequences(filePath)
		except Exception as e:
			messagebox.showerror("Error", f"An error occurred while reading the file: {e}")
	else:
			messagebox.showwarning("No file", "Please select a file")
	return secRecList, filePath

def loadModel(default:False, append=False):	
	seqRecList, filePath= loadFile(default )
	if append and Model.modelInstance!=None:
		for newRecord in seqRecList:
			Model.modelInstance.sequenceRecordList.append(newRecord)
	else:
		newModel=Model(filePath,seqRecList)
		Model.modelInstance=newModel
	# Model.modelInstance.dumpModel("in main")
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def findBoringRegions(longString, locations, extraLocation=None)->list[tuple[int,int]]:
	# Initialize an empty list to store non-overlapping regions
	boringRegions:list[tuple[int,int]] = []
	# Assume that the string starts from index 0
	currentEnd = 0
	stringLength: int = len(longString)
	for start, end in locations:
		# If there is a gap between the current end and the start of the next substring, add the gap as a region
		if start > currentEnd:
			boringRegions.append((currentEnd, start))
		# Update the currentEnd to be the maximum of the currentEnd and the current substring's end
		currentEnd: int = max(currentEnd, end)
	# If there's any remaining region after the last substring, add it
	if currentEnd < stringLength:
		boringRegions.append((currentEnd, stringLength))    
	return boringRegions

def isExcitingLetterIndex( index, boringRegions, potentialAnealedPrimer=None): 
	exciting= not any(start <= index < end for start, end in boringRegions)		
	if not exciting and potentialAnealedPrimer:
		exciting=potentialAnealedPrimer.features[0].location.start <= index < potentialAnealedPrimer.features[0].location.end
	return exciting		

# def drawCanvasCircle(canvas:Canvas):
#         canvas.create_oval(100, 150, 200, 250, outline="blue", width=2)

def canvasZoom(zoomin):# 1 for  zoom In or bigger
    if zoomin!=True and gl.prefs.get_preference_value("fontSize")>3: # ZOOM OUT no to negative font sizes
            gl.prefs.setFontSize(gl.prefs.get_preference_value("fontSize")-1)    
            drawCanvas()
    else:   
        gl.prefs.setFontSize(gl.prefs.get_preference_value("fontSize")+1 )          
        drawCanvas()

def calculateBaseRectangleSymbolXPixelSize(fontSize):
	return fontSize+horizontalPixelsMargin

def calculateBaseRectangleSymbolYPixelSize(fontSize):
	return fontSize+horizontalPixelsMargin

def loadAndSeparateEmblSequences(emblName:str)->list[MySeqRecord]:
	#fIn=open(embl,'r')
	sequenceRecordIterator=SeqIO.parse(emblName, "embl")
	sequenceRecordList:list[MySeqRecord]=[]
	for record in sequenceRecordIterator:
		if record.features==None:
			record.features=list()
		singleStranded=True if record.annotations.get("molecule_type")=="ss-DNA" else False
		if singleStranded:
			myRecord=MySeqRecord(record,True, True, primer=False)
			sequenceRecordList.append(myRecord)	
		else:
			myRecord53=MySeqRecord(record,False, True, primer=False)
			sequenceRecordList.append(myRecord53)
			threeTo5Record=deepcopy(myRecord53)
			threeTo5Record.seq=threeTo5Record.seq.complement()
			myRecord35=MySeqRecord(threeTo5Record,False, False, primer=False)
			sequenceRecordList.append(myRecord35)		
			myRecord53.hybridizedTo=myRecord35
			myRecord35.hybridizedTo=myRecord53
	return 	sequenceRecordList 

def loadFastas():
	fName="short.fa" #sys.argv[1]
	fIn=open(fName,'r')
	sequenceRecordIterator=SeqIO.parse(fName, "fasta")
	sequenceRecordList=[]
	for record in sequenceRecordIterator:
		print("Fasta record name:%s length:%i Sequence: %s" % (record.id, len(record), record.seq))
		sequenceRecordList.append(record)
	return 	sequenceRecordList

def seqToString(seq:Seq)->str:
	return seq._data.decode('ASCII') # from bytes to stringead and strip any extra whitespace or newline characters

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

def shiftLeftShrinkedFeaturesLocations(i:int, record:MySeqRecord):
	shrinkedFeatures=record.shrinkedFeatures
	for f in range(len(shrinkedFeatures)):		
		if record.features[f].location.start >i:
			if  isinstance(shrinkedFeatures[f].location,CompoundLocation):
				parts = list()
				for loc in shrinkedFeatures[f].location.parts:
					loc=loc-1
					parts.append(loc)
				shrinkedFeatures[f].location=CompoundLocation(parts)
			else:
				shrinkedFeatures[f].location=shrinkedFeatures[f].location-1
	if record.hybridizedTo:
		if record.hybridizedTo.features[0].location.start >i:
			record.hybridizedTo.features[0].location=record.hybridizedTo.features[0].location-1

#draw features from original or cached features
def drawFeatures(canvas: Canvas, mySequenceRecord: MySeqRecord, xStart: int, yStart: int, baseRectangleSymbolXPixelSize: int, baseRectangleSymbolYPixelSize: int, verticalSequenceSpacing: int, font: tuple[str, int], shrink: bool):
	if mySequenceRecord.isPrimer:
		source: list[SeqFeature] = mySequenceRecord.features
	else:# regular sequences
		if shrink:
			source = mySequenceRecord.shrinkedFeatures
		else:
			source: list[SeqFeature] = mySequenceRecord.features
	for feature in source:
		if ((feature.location.strand == 1 and mySequenceRecord.fiveTo3) or (feature.location.strand == -1 and not mySequenceRecord.fiveTo3) or feature.location.strand==None):
			loc: Location = feature.location
			drawTextInRectangle(feature.qualifiers.get("label")[0], canvas, xStart + baseRectangleSymbolXPixelSize * loc.start, yStart, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, 'white', loc.end - loc.start, font)

def drawPrimer(canvas:Canvas,mySequenceRecord:MySeqRecord,xStart:int, yStart:int)->int:
	fontSize:int=gl.prefs.get_preference_value(preference_name="fontSize")
	shrink:int=gl.prefs.get_preference_value(preference_name="shrink")	
	font:tuple[str,int]=(gl.prefs.get_preference_value(preference_name="fontName"),fontSize)
	coloredBases:bool=gl.prefs.get_preference_value(preference_name="coloredBases")
	rotated:bool=gl.prefs.get_preference_value(preference_name="rotated")
	baseRectangleSymbolXPixelSize:int=calculateBaseRectangleSymbolXPixelSize(fontSize) #in pixels
	baseRectangleSymbolYPixelSize:int=calculateBaseRectangleSymbolYPixelSize(fontSize) #in pixels	
	verticalSequenceSpacing:int=gl.prefs.get_preference_value(preference_name="verticalSequenceSpacing")+5 # blank space between 2 sequences
	xStart=xStart	+mySequenceRecord.xStartOffset*	baseRectangleSymbolXPixelSize
	x: int=xStart
	seq:Seq=mySequenceRecord.seq
	upsideDownLetter:bool=not mySequenceRecord.fiveTo3 and rotated
	dnaSequenceStr=seqToString(seq)	
	featureYStart, sequenceYStart, bandYEnd = calculateYs(canvas, mySequenceRecord, xStart, yStart, baseRectangleSymbolYPixelSize, verticalSequenceSpacing)		
	i=0
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]
		color="white" if not coloredBases else gl.prefs.get_preference_value(letter)
		drawBase(letter, canvas, x, sequenceYStart, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize,
		color=color, font=font, upsideDownLetter=upsideDownLetter)
		x += baseRectangleSymbolXPixelSize # Move to the next position
	drawFeatures(canvas, mySequenceRecord, xStart, featureYStart, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, shrink)
	def enhancedButtonAction(event: tk.Event) -> None:
		eb.handle_click(event)

	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton(canvas, mySequenceRecord.description[:2], 0, yStart, enhancedButtonAction, labelHeightPx=bandYEnd-yStart)	

	return  x,bandYEnd

def drawStrand(canvas:Canvas,mySequenceRecord:MySeqRecord,xStart:int, yStart:int)->int:
	fontSize:int=gl.prefs.get_preference_value(preference_name="fontSize")
	shrink:int=gl.prefs.get_preference_value(preference_name="shrink")	
	font:tuple[str,int]=(gl.prefs.get_preference_value(preference_name="fontName"),fontSize)
	# gl.prefs.dump()
	coloredBases:bool=gl.prefs.get_preference_value(preference_name="coloredBases")
	rotated:bool=gl.prefs.get_preference_value(preference_name="rotated")
	baseRectangleSymbolXPixelSize:int=calculateBaseRectangleSymbolXPixelSize(fontSize) #in pixels
	baseRectangleSymbolYPixelSize:int=calculateBaseRectangleSymbolYPixelSize(fontSize) #in pixels	
	verticalSequenceSpacing:int=gl.prefs.get_preference_value(preference_name="verticalSequenceSpacing")+5 # blank space between 2 sequences
	featureYStart, sequenceYStart, bandYEnd = calculateYs(canvas, mySequenceRecord, xStart, yStart, baseRectangleSymbolYPixelSize, verticalSequenceSpacing)
	x: int=xStart		
	seq:Seq=mySequenceRecord.seq
	if shrink:
		mySequenceRecord.shrinkedFeatures= deepcopy(mySequenceRecord.features)
		if mySequenceRecord.hybridizedTo:
			mySequenceRecord.gostFeatures= deepcopy(mySequenceRecord.features)
	dnaSequenceStr=seqToString(seq)
	if mySequenceRecord.features:
		nonOverlappingRegions:list[tuple[int,int]]=findBoringRegions(dnaSequenceStr, [(feature.location.start, feature.location.end) for feature in mySequenceRecord.features])
	spamCount=1
	upsideDownLetter:bool=not mySequenceRecord.fiveTo3 and rotated
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]
		excitingLetter: bool=isExcitingLetterIndex(i, nonOverlappingRegions, mySequenceRecord.hybridizedTo)
		if excitingLetter==True:
			spamCount=0
		if not shrink or spamCount<=3:
			color=None
			if not coloredBases:
				if not shrink:
					color="white"
				elif not excitingLetter:
					color="grey"
			else:
				if shrink and not excitingLetter:
					color="grey"
				else:
					color =gl.prefs.get_preference_value(letter)

			drawBase(letter, canvas, x, sequenceYStart, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize,
			color=color, font=font, upsideDownLetter=upsideDownLetter)
			spamCount+=1
			x += baseRectangleSymbolXPixelSize # Move to the next position
		elif shrink:
			# skip a letter at index i so decrease by 1 all the following locations
			shiftLeftShrinkedFeaturesLocations(i, record=mySequenceRecord )

	# Replace the placeholder with the call to the new method
	drawFeatures(canvas, mySequenceRecord, xStart, featureYStart, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, shrink)
	
	def enhancedButtonAction(event: tk.Event) -> None:
		eb.handle_click(event)

	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton(canvas, mySequenceRecord.description[:2], 0, yStart, enhancedButtonAction, labelHeightPx=bandYEnd-yStart)	
	return  x,bandYEnd

#calculate the yop relative Ys for features and primers and the final y
def calculateYs(canvas, mySequenceRecord, xStart, yStart, baseRectangleSymbolYPixelSize, verticalSequenceSpacing):
	if mySequenceRecord.fiveTo3:
		featureYStart=yStart+verticalSequenceSpacing
		sequenceYStart=yStart+verticalSequenceSpacing+baseRectangleSymbolYPixelSize
		if mySequenceRecord.hybridizedTo:
			bandEnd=yStart+verticalSequenceSpacing+2*baseRectangleSymbolYPixelSize
		else:
			bandEnd=yStart+verticalSequenceSpacing+2*baseRectangleSymbolYPixelSize+verticalSequenceSpacing
	else: #  3 to 5
		if mySequenceRecord.hybridizedTo:
			sequenceYStart=yStart
			featureYStart=yStart+baseRectangleSymbolYPixelSize
			bandEnd=yStart+verticalSequenceSpacing+2*baseRectangleSymbolYPixelSize
		else:
			sequenceYStart=yStart+verticalSequenceSpacing
			featureYStart=yStart+verticalSequenceSpacing+baseRectangleSymbolYPixelSize			
			bandEnd=yStart+verticalSequenceSpacing+2*baseRectangleSymbolYPixelSize+verticalSequenceSpacing
	return featureYStart,sequenceYStart,bandEnd

def drawTextInRectangle(tex:str,canvas:Canvas, xLeft, yTop, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, length,font,rotated=None):
	canvas.create_rectangle( xLeft, yTop, xLeft + length*baseRectangleSymbolXPixelSize ,yTop + baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if rotated:
		# The default anchor for text in Tkinter's Canvas.create_text() method is "center". T
		canvas.create_text(xLeft+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, yTop+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(xLeft+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, yTop+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black")	
	# upside_down_text = tk.Text(canvas, width=gl.fontSize, height=1 , wrap=tk.WORD, bd=1, relief="solid", highlightbackground="red",font, padx=0, pady=-3) # height is the number of lines
	# upside_down_text.insert(tk.END, tex)  # Insert the 
	# textpushDown=baseRectangleSymbolYPixelSize+10
	# canvas.create_window(x+length*baseRectangleSymbolXPixelSize/2, textpushDown+y, width=length*baseRectangleSymbolXPixelSize+1, height=baseRectangleSymbolYPixelSize+3,  window=upside_down_text)
	
# Define functions to draw each DNA base x,y relative from upper 	left corner
def drawBase(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, font, upsideDownLetter=None):
	canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	# 
	if upsideDownLetter:
		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")

# Define functions to draw each DNA base x,y relative from upper left corner
def drawTextInRectangleWithoutWidgets(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, col, font, length,rotated=None):
	canvas.create_rectangle( x, y, x + length*baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=col, outline="black" )
	if rotated:
		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=base, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=base, font=font, fill="black")	

# Define functions to draw each DNA base x,y relative from upper left corner
def drawBaseWithoutWidgets(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, col, font, coloredBases, rotated=None):
	if coloredBases==False:
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill="white", outline="black" )
	else:
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=col, outline="black" )
	if rotated:
		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(x+baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")	

def printCanvas(canvas):
	# filePath=str(Path(__file__).resolve().parent)+"/junk"
	# canvas.postscript(file=filePath, colormode='color')
 	# # Get the canvas's total scrollable area (scrollregion)
    scrollregion = canvas.bbox("all")  # Returns (x1, y1, x2, y2)    
    if scrollregion:
        x1, y1, x2, y2 = scrollregion
        
        # Get the position of the canvas on the screen
        canvas_x = canvas.winfo_rootx()
        canvas_y = canvas.winfo_rooty()

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