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
# from upsideDownText import UpsideDownText
# mine
import gl
from model import Model
#######
# some 'constants' that can be changed if push comes to shove !
# some initial values that might remain constant
canvasHorizontalMargin = 0 # blank space from the left edge of canvas
#canvasVerticalMargin = 30
horizontalPixelsMargin=2 # total room between the base letter and it's holding box

def drawCanvas( canvas:Canvas)->int:
	verticalSequenceSpacing:int=gl.prefs.get_preference_value(preference_name="verticalSequenceSpacing")+5 # blank space between 2 sequences
	fontSize:int=gl.prefs.get_preference_value(preference_name="fontSize")
	font:tuple[str,int]=(gl.prefs.get_preference_value(preference_name="fontName"),fontSize)
	# gl.prefs.dump()
	coloredBases:bool=gl.prefs.get_preference_value(preference_name="coloredBases")
	rotated:bool=gl.prefs.get_preference_value(preference_name="rotated")
	baseRectangleSymbolXPixelSize:int=calculateBaseRectangleSymbolXPixelSize(fontSize) #in pixels
	baseRectangleSymbolYPixelSize:int=calculateBaseRectangleSymbolYPixelSize(fontSize) #in pixels
	yPos:int = verticalSequenceSpacing  # y is 0 at top and increases downwards
	
	# Clear any previous drawings
	canvas.delete("all")
	sequenceWidth=0
	sequenceIndex=0
	for sequenceRecord in Model.modelInstance.sequenceRecordList:
		newSequenceWidth,yFinal=drawSequenceRecord(canvas, sequenceRecord, sequenceIndex, canvasHorizontalMargin,yPos,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize,verticalSequenceSpacing, font,coloredBases, rotated)
		if newSequenceWidth>sequenceWidth:
			sequenceWidth=newSequenceWidth
		yPos=yFinal
		sequenceIndex+=1
	#
	# Adjust the scrollable region based on the length of the string
	canvas.config(scrollregion=(0, 0, sequenceWidth,100))  # Set the scrollable area
	return 2*canvasHorizontalMargin+sequenceWidth

# the ID line is parsed in  C:\a\diy\pythonProjects\DNAPrinting\.venv\Lib\site-packages\Bio\GenBank\Scanner.py EmblScanner._feed_first_line and the parsing in line 788
def loadFile(default=False)->tuple[list[MySeqRecord],str]:	
	if default:
		filePath=str(Path(__file__).resolve().parent)+gl.prefs.get_preference_value("defaultTestFileValue")
	else:
		filePath = filedialog.askopenfilename(title="Open EMBL File", filetypes=[("EMBL Files", "*.embl"), ("All Files", "*.*")])    
	secRecList:list[MySeqRecord]		=None
	if filePath:
		try:
			secRecList=loadEmblSequences(filePath)
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

	Model.modelInstance.dumpModel("in main")
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def findNonOverlappingRegions(longString, locations)->list[tuple[int,int]]:
    # Initialize an empty list to store non-overlapping regions
    nonOverlappingRegions:list[tuple[int,int]] = []
    # Assume that the string starts from index 0
    currentEnd = 0
    stringLength: int = len(longString)
    for start, end in locations:
        # If there is a gap between the current end and the start of the next substring, add the gap as a region
        if start > currentEnd:
            nonOverlappingRegions.append((currentEnd, start))
        # Update the currentEnd to be the maximum of the currentEnd and the current substring's end
        currentEnd: int = max(currentEnd, end)
    # If there's any remaining region after the last substring, add it
    if currentEnd < stringLength:
        nonOverlappingRegions.append((currentEnd, stringLength))    
    return nonOverlappingRegions

def isIndexOverlapping( index, nonOverlappingRegions): 
	return not any(start <= index < end for start, end in nonOverlappingRegions)		

def drawCanvasCircle(canvas:Canvas):
        canvas.create_oval(100, 150, 200, 250, outline="blue", width=2)

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

def loadEmblSequences(emblName:str)->list[MySeqRecord]:
	#fIn=open(embl,'r')
	sequenceRecordIterator=SeqIO.parse(emblName, "embl")
	sequenceRecordList:list[MySeqRecord]=[]
	for record in sequenceRecordIterator:
		myRecord=MySeqRecord(record)
		myRecord.singleStranded=True if record.annotations.get("molecule_type")=="ss-DNA" else False
		sequenceRecordList.append(myRecord)
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

def checkTranslation( inSeq, outSeq):
	"""Checks the translation of the engineered sequence against the wild-type sequence"""
	myInSeq=MySeqRecord(Seq(inSeq))
	myOutSeq=MySeqRecord(Seq(outSeq))
	if myInSeq.translate().seq==myOutSeq.translate().seq:
		successFlag=True
	else:
		successFlag=False
	return successFlag

def drawSequenceRecord(canvas:Canvas,sequenceRecord:MySeqRecord,sequenceIndex, xStart:int, yStart:int,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, coloredBases, rot):
	shri=gl.prefs.get_preference_value(preference_name="shrink")
	x=drawStrand( canvas, sequenceRecord=sequenceRecord,sequenceIndex=sequenceIndex, xStart=canvasHorizontalMargin,yStart=yStart,baseRectangleSymbolXPixelSize=baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize=baseRectangleSymbolYPixelSize,verticalSequenceSpacing=verticalSequenceSpacing, font=font, coloredBases=coloredBases, shrink=shri)	
	yFin=yStart+verticalSequenceSpacing+baseRectangleSymbolYPixelSize
	if not sequenceRecord.singleStranded:
		x=drawStrand( canvas, sequenceRecord, sequenceIndex,canvasHorizontalMargin,yStart+2*verticalSequenceSpacing,baseRectangleSymbolXPixelSize,
			   baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, coloredBases,
			    shri,complemented=True, rotated=True if rot else False)
		yFin=yFin+verticalSequenceSpacing+baseRectangleSymbolYPixelSize
	return  x,yFin

def adjustCachedFeatureLocations(i:int, record:MySeqRecord):
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
			

def drawStrand(canvas:Canvas,sequenceRecord:MySeqRecord,sequenceIndex:int,xStart:int, yStart:int,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize, verticalSequenceSpacing,  font, coloredBases,
			    shrink=None, complemented=False, rotated=None)->int:
	x=xStart		
	if not isinstance(sequenceRecord,SeqRecord):
		print ("stop here")
	seq:Seq=sequenceRecord.seq
	if shrink and not complemented: # build the cache twice
		sequenceRecord.shrinkedFeatures= deepcopy(sequenceRecord.features)

	if complemented:
		dnaSequenceStr=seqToString(seq.complement())
	else:
		dnaSequenceStr=seqToString(seq)
	nonOverlappingRegions:list[tuple[int,int]]=findNonOverlappingRegions(dnaSequenceStr, [(feature.location.start, feature.location.end) for feature in sequenceRecord.features])
	spamCount=1
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]
		overlapping: bool=isIndexOverlapping(i, nonOverlappingRegions)
		if overlapping==True:
			spamCount=0
		if not shrink or spamCount<=3:
			drawBase(letter, canvas, x, verticalSequenceSpacing+yStart, baseRectangleSymbolXPixelSize,
				baseRectangleSymbolYPixelSize,gl.prefs.get_preference_value(letter), font, coloredBases if not shrink  else True if overlapping else False, overlapping, rotated)
			spamCount+=1
			x += baseRectangleSymbolXPixelSize # Move to the next position
		elif shrink and not complemented:
			# skip a letter at intex i
			adjustCachedFeatureLocations(i, record=sequenceRecord )
	#draw features from original or cached features
	if shrink:
		source=sequenceRecord.shrinkedFeatures
	else:
		source=sequenceRecord.features

	for feature in source:
		if((feature.location.strand==1 and rotated is None) or (feature.location.strand==-1 and rotated==True)):
			# print("Feature in drawSequenceRecord:","id:",feature.id, feature.qualifiers.get("label")[0] , "location:",feature.location,"strand",feature.location.strand, "type", feature.type,)
			loc:Location=feature.location
			topY:int=verticalSequenceSpacing+yStart
			#labelId=canvas.create_rectangle( xStart+baseRectangleSymbolXPixelSize*loc.start, topY+baseRectangleSymbolYPixelSize, xStart+baseRectangleSymbolXPixelSize*loc.end,topY+ 2*baseRectangleSymbolYPixelSize , fill="white", outline="black" )
			# Get the width of the second string
			# bbox = canvas.bbox(labelId)  # bbox returns (x1, y1, x2, y2)
			# width = bbox[2] - bbox[0]  # The width is the difference between x2 and x1
			# print("Location:", loc.start, loc.end, loc.strand)
			drawTextInRectangle(feature.qualifiers.get("label")[0],canvas, xStart+baseRectangleSymbolXPixelSize*loc.start, topY+baseRectangleSymbolYPixelSize, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, 'white',loc.end-loc.start, font, )
	return  x

def drawTextInRectangle(tex:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, length,font,rotated=None):
	canvas.create_rectangle( x, y, x + length*baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if rotated:
		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black", angle=180)	
	else:
		canvas.create_text(x+baseRectangleSymbolXPixelSize+ length*baseRectangleSymbolXPixelSize/2, y+baseRectangleSymbolYPixelSize/2 , text=tex, font=font, fill="black")	
	# upside_down_text = tk.Text(canvas, width=gl.fontSize, height=1 , wrap=tk.WORD, bd=1, relief="solid", highlightbackground="red",font, padx=0, pady=-3) # height is the number of lines
	# upside_down_text.insert(tk.END, tex)  # Insert the 
	# textpushDown=baseRectangleSymbolYPixelSize+10
	# canvas.create_window(x+length*baseRectangleSymbolXPixelSize/2, textpushDown+y, width=length*baseRectangleSymbolXPixelSize+1, height=baseRectangleSymbolYPixelSize+3,  window=upside_down_text)
	
# Define functions to draw each DNA base x,y relative from upper left corner
def drawBase(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, font, colored, isIndexOverlapping=None, rotated=None):
	if not colored:
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=("white" if isIndexOverlapping else "grey"), outline="black" )
	else:			 #colored bases
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	# 
	if rotated:
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