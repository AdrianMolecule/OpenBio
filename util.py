"""Utility methods"""
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
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location
# from upsideDownText import UpsideDownText
# mine
import gl
import model
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
	baseRectangleSymbolXPixelSize:int=calculateBaseRectangleSymbolXPixelSize(fontSize) #in pixels
	baseRectangleSymbolYPixelSize:int=calculateBaseRectangleSymbolYPixelSize(fontSize) #in pixels
	yPos:int = verticalSequenceSpacing  # y is 0 at top and increases downwards
	# Clear any previous drawings
	canvas.delete("all")
	sequenceWidth=0
	for sequenceRecord in model.sequenceRecordList:
		newSequenceWidth,yFinal=drawSequenceRecord(canvas, sequenceRecord, canvasHorizontalMargin,yPos,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize,verticalSequenceSpacing, font,coloredBases)#, True)
		if newSequenceWidth>sequenceWidth:
			sequenceWidth=newSequenceWidth
		yPos=yFinal
	#
	# Adjust the scrollable region based on the length of the string
	canvas.config(scrollregion=(0, 0, sequenceWidth,100))  # Set the scrollable area
	return 2*canvasHorizontalMargin+sequenceWidth

# the ID line is parsed in  C:\a\diy\pythonProjects\DNAPrinting\.venv\Lib\site-packages\Bio\GenBank\Scanner.py EmblScanner._feed_first_line and the parsing in line 788
def loadFile()->list[SeqRecord]:	
	defaultTestFileValue:int=gl.prefs.get_preference_value(preference_name="defaultTestFileValue") 
	if defaultTestFileValue is None or defaultTestFileValue=="":	
		filePath = filedialog.askopenfilename(title="Open EMBL File", filetypes=[("EMBL Files", "*.embl"), ("All Files", "*.*")])    
	else:
		filePath=str(Path(__file__).resolve().parent)+gl.prefs.get_preference_value("defaultTestFileValue")
	if filePath:
		try:
			secRecList:list[SeqRecord]=loadEmblSequences(filePath)
			print(secRecList,secRecList.__class__)                
			model.loadedFileName=filePath
			return secRecList
		except Exception as e:
			messagebox.showerror("Error", f"An error occurred while reading the file: {e}")
	else:
			messagebox.showwarning("No file", "Please select a file")

def loadModel():	
	seqRecList:list[SeqRecord]= loadFile()
	print ("New load of model", seqRecList )
	model.sequenceRecordList=seqRecList
	model.dumpModel("in main")
	model.appendSequenceRecord(newSequenceRecord=SeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

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
	return any(start <= index < end for start, end in nonOverlappingRegions)		

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

def loadEmblSequences(emblName:str)->list[SeqRecord]:
	#fIn=open(embl,'r')
	sequenceRecordIterator=SeqIO.parse(emblName, "embl")
	sequenceRecordList=[]
	for record in sequenceRecordIterator:
		sequenceRecordList.append(record)
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
	myInSeq=SeqRecord(Seq(inSeq))
	myOutSeq=SeqRecord(Seq(outSeq))
	if myInSeq.translate().seq==myOutSeq.translate().seq:
		successFlag=True
	else:
		successFlag=False
	return successFlag

def drawSequenceRecord(canvas:Canvas,sequenceRecord:SeqRecord,xStart:int, yStart:int,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, coloredBases, rot=None):
	shri=gl.prefs.get_preference_value(preference_name="shrink")
	"""Draws a sequence record on the canvas"""
	x=drawStrand( canvas, sequenceRecord, canvasHorizontalMargin,yStart,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize,verticalSequenceSpacing, font, coloredBases, shri)	
	yFin=yStart+verticalSequenceSpacing+baseRectangleSymbolYPixelSize
	if sequenceRecord.annotations.get("molecule_type")!="ss-DNA":
		x=drawStrand( canvas, sequenceRecord, canvasHorizontalMargin,yStart+2*verticalSequenceSpacing,baseRectangleSymbolXPixelSize,
			   baseRectangleSymbolYPixelSize, verticalSequenceSpacing, font, coloredBases,
			    shri,complemented=True, rotated=True if rot else False)
		yFin=yFin+verticalSequenceSpacing+baseRectangleSymbolYPixelSize
	return  x,yFin

def drawStrand(canvas:Canvas,sequenceRecord:SeqRecord,xStart:int, yStart:int,baseRectangleSymbolXPixelSize,baseRectangleSymbolYPixelSize, verticalSequenceSpacing,  font, coloredBases,
			    shrink=None, complemented=False, rotated=None):
	x=xStart
	seq:Seq=sequenceRecord._seq
	if complemented:
		dnaSequenceStr=seqToString(seq.complement())
	else:
		dnaSequenceStr=seqToString(seq)
	nonOverlappingRegions:list[tuple[int,int]]=findNonOverlappingRegions(dnaSequenceStr, [(feature.location.start, feature.location.end) for feature in sequenceRecord.features])
	spamCount=0
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]
		overlapping: bool=isIndexOverlapping(i, nonOverlappingRegions)
		if overlapping==True:
			spamCount=0
		if True:#not shrink or spamCount<4:
			drawBase(letter, canvas, x, verticalSequenceSpacing+yStart, baseRectangleSymbolXPixelSize,
				baseRectangleSymbolYPixelSize,gl.prefs.get_preference_value(letter), font, coloredBases, overlapping, rotated)
			spamCount+=1
		x += baseRectangleSymbolXPixelSize # Move to the next position
	#draw features
	for feature in sequenceRecord.features:
		if((feature.location.strand==1 and rotated is None) or (feature.location.strand==-1 and rotated==True)):
			print("Feature in drawSequenceRecord:","id:",feature.id, feature.qualifiers.get("label")[0] , "location:",feature.location,"strand",feature.location.strand, "type", feature.type,)
			loc:Location=feature.location
			topY:int=verticalSequenceSpacing+yStart
			#labelId=canvas.create_rectangle( xStart+baseRectangleSymbolXPixelSize*loc.start, topY+baseRectangleSymbolYPixelSize, xStart+baseRectangleSymbolXPixelSize*loc.end,topY+ 2*baseRectangleSymbolYPixelSize , fill="white", outline="black" )
			# Get the width of the second string
			# bbox = canvas.bbox(labelId)  # bbox returns (x1, y1, x2, y2)
			# width = bbox[2] - bbox[0]  # The width is the difference between x2 and x1
			print("Location:", loc.start, loc.end, loc.strand)
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
def drawBase(base:str,canvas:Canvas, x, y, baseRectangleSymbolXPixelSize, baseRectangleSymbolYPixelSize, color, font, coloredBases, isIndexOverlapping=None, rotated=None):
	if coloredBases==False:
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=("grey" if isIndexOverlapping else "white"), outline="black" )
	else:			 #colored bases
		canvas.create_rectangle( x, y, x + baseRectangleSymbolXPixelSize ,y + baseRectangleSymbolYPixelSize , fill=("grey" if isIndexOverlapping else color), outline="black" )
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