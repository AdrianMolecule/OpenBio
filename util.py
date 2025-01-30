from copy import deepcopy
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
from myseqrecord import MySeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition, UnknownPosition, Location
# mine
import gl
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
	drawMask(yFinal)
	# Adjust the scrollable region based on the length of the string
	gl.canvas.config(scrollregion=(0, 0, sequenceWidth, yFinal+40))  # Set the scrollable area
	return 2*gl.canvasLeftPadding+sequenceWidth

#draw features from original or cached features
def drawFeatures( mySequenceRecord: MySeqRecord, yStart: int, font: tuple[str, int], bgColor):
	for feature in mySequenceRecord.features:
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

def drawPrimer(mySequenceRecordPrimer:MySeqRecord, yStart:int)->int:	
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	verticalSequenceSpacing:int=gl.verticalSequenceSpacing+5 # blank space between 2 sequences Adrian
	upsideDownLetter:bool=not mySequenceRecordPrimer.fiveTo3 and gl.rotated
	dnaSequenceStr=seqToString(mySequenceRecordPrimer.seq)	
	featureYStart, sequenceYStart, bandYEnd = calculateYs( mySequenceRecordPrimer,  yStart)		
	i=0
	for i in range(len(dnaSequenceStr)):
		letter: str=dnaSequenceStr[i]		
		if gl.shrink and mySequenceRecordPrimer.hybridizedToStrand: # non attached primers do not have a clear x position as they float in the liquid
			xLett=gl.canvasLeftPadding + (mySequenceRecordPrimer.xStartOffsetAsLetters+i+gl.maskSkipped[mySequenceRecordPrimer.xStartOffsetAsLetters+i])*gl.baseRectangleSymbolXPixelSize
		else:
			xLett=gl.canvasLeftPadding + (mySequenceRecordPrimer.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize
		color="white" if not gl.coloredBases else gl.prefs.getPreferenceValue(letter)
		drawBase(letter, xLett, sequenceYStart,  color=color, font=font, upsideDownLetter=upsideDownLetter,fiveTo3=mySequenceRecordPrimer.fiveTo3)
	drawFeatures( mySequenceRecordPrimer, featureYStart,  font, 'pink')
	if gl.shrink:
		maxX=gl.canvasLeftPadding + (len(dnaSequenceStr)-(gl.maskSkipped[i+mySequenceRecordPrimer.xStartOffsetAsLetters])*gl.baseRectangleSymbolXPixelSize)
	else:	
		maxX=gl.canvasLeftPadding + (len(dnaSequenceStr)-(gl.maskSkipped[i])*gl.baseRectangleSymbolXPixelSize)
		
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton(mySequenceRecordPrimer.description[:8], 1, yStart,mySequenceRecordPrimer,  labelHeightPx=bandYEnd-yStart)	
	return  maxX,bandYEnd


def drawStrand(mySequenceRecord:MySeqRecord, yStart:int)->int:
	font:tuple[str,int]=(gl.fontName,gl.fontSize)
	# gl.prefs.dump()
	featureYStart, sequenceYStart, bandYEnd = calculateYs( mySequenceRecord,  yStart)
	seq:Seq=mySequenceRecord.seq
	dnaSequenceStr=seqToString(seq)
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
				drawBase(letter,  xLett, sequenceYStart, color=color, font=font, upsideDownLetter=upsideDownLetter, fiveTo3=mySequenceRecord.fiveTo3)
		else: # NO SHRINK
			xLett=gl.canvasLeftPadding + (mySequenceRecord.xStartOffsetAsLetters+i)*gl.baseRectangleSymbolXPixelSize	
			if not gl.coloredBases:
				color="white"
			else:
				color =gl.prefs.getPreferenceValue(letter)
			drawBase(letter,  xLett, sequenceYStart, color=color, font=font, upsideDownLetter=upsideDownLetter, fiveTo3=mySequenceRecord.fiveTo3)

	# Replace the placeholder with the call to the new method
	drawFeatures(mySequenceRecord, featureYStart,  font,  "white")
	
	# Create labels using the createLabel method
	eb:EnhancedButton=EnhancedButton( mySequenceRecord.description[:8], 1, yStart,mySequenceRecord, labelHeightPx=bandYEnd-yStart)	
	if mySequenceRecord.fiveTo3:
		if  gl.hydrogen : 
			bandYEnd+=gl.hydrogenLinesLength# add a gap to draw the hydrogen bond lines	
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
def drawBase(base:str, x, y, color, font, upsideDownLetter=None, fiveTo3=None):
	gl.canvas.create_rectangle( x, y, x + gl.baseRectangleSymbolXPixelSize ,y + gl.baseRectangleSymbolYPixelSize , fill=color, outline="black" )
	if upsideDownLetter:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black", angle=180)			
	else:
		textId=gl.canvas.create_text(x+gl.baseRectangleSymbolXPixelSize/2, y+gl.baseRectangleSymbolYPixelSize/2+1 , text=base, font=font, fill="black")
	if fiveTo3 and gl.hydrogen:
		drawHydrogenBondLines( base, x+1,y + gl.baseRectangleSymbolYPixelSize+1)

	gl.canvas.tag_bind(textId, "<Button-1>", lambda event: clickOnSeqRecord(event, gl.canvas, None))

def drawHydrogenBondLines( base:str, x: int, topY: int)->int:#y is always top
	lineLength = gl.hydrogenLinesLength-1  # Length of each vertical line
	lineSpacing = 2  # Spacing between the lines
	for i in range(3 if base=="C" or base=="G" else 2):
		gl.canvas.create_line(x+3, topY, x+3, topY + lineLength, fill="white")
		x += lineSpacing	
	return lineLength

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


# the ID line is parsed in  C:\a\diy\pythonProjects\DNAPrinting\.venv\Lib\site-packages\Bio\GenBank\Scanner.py EmblScanner._feed_first_line and the parsing in line 788
def loadFile(default=False)->tuple[list[MySeqRecord],str]:	
	currentDir=str(Path(__file__).resolve().parent)
	if default:
		filePath=currentDir+"/samples/"+gl.prefs.getPreferenceValue("defaultTestFileValue")
	else:
		if not gl.debug:
			filePath = filedialog.askopenfilename(title="Open EMBL File",filetypes=[("EMBL Files", "*.embl"), ("All Files", "*.*")])  
		else:
			filePath = filedialog.askopenfilename(title="Open EMBL File",  initialdir=currentDir+"/samples/",filetypes=[("EMBL Files", "*.embl"),("All Files", "*.*")])    
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
	if seqRecList is None or len(seqRecList)==0:
		messagebox.showerror("No Sequences", f" Please select a file that has at least one sequence") 
		return None	
	if append and Model.modelInstance!=None:
		for newRecord in seqRecList:
			Model.modelInstance.sequenceRecordList.append(newRecord)
	else:
		newModel=Model(filePath,seqRecList)
		Model.modelInstance=newModel
	# Model.modelInstance.dumpModel("in main")
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

# def drawCanvasCircle(canvas:Canvas):
#         canvas.create_oval(100, 150, 200, 250, outline="blue", width=2)

def canvasZoom(zoomin):# 1 for  zoom In or bigger
	if zoomin!=True and gl.prefs.getPreferenceValue("fontSize")>3: # ZOOM OUT no to negative font sizes
		gl.prefs.setFontSize(gl.fontSize-1)    
		refresh()
	else:   	
		gl.prefs.setFontSize(gl.fontSize+1 )          
	refresh()

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
			myRecord.singleStranded=True
			sequenceRecordList.append(myRecord)	
		else:
			myRecord53=MySeqRecord(record,False, True, primer=False)
			sequenceRecordList.append(myRecord53)
			threeTo5Record=deepcopy(myRecord53)
			threeTo5Record.seq=threeTo5Record.seq.complement()
			myRecord35=MySeqRecord(threeTo5Record,False, False, primer=False)
			sequenceRecordList.append(myRecord35)		
			myRecord53.hybridizedToStrand=myRecord35
			myRecord35.hybridizedToStrand=myRecord53
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

def getLabel (feature:SeqFeature):
	return feature.qualifiers.get("label")


def drawMask( yStart)->int:
	font:tuple[str,int]=(gl.fontName,6)
	x=gl.canvasLeftPadding
	y=yStart	+gl.baseRectangleSymbolYPixelSize
	if gl.ruler:
		maxLen=0
		for i, sequenceRecord in enumerate(Model.modelInstance.sequenceRecordList):
				sequenceRecord:MySeqRecord
				if len(sequenceRecord.seq)+sequenceRecord.xStartOffsetAsLetters>maxLen:
					maxLen=len(sequenceRecord.seq)+sequenceRecord.xStartOffsetAsLetters
		for p in range(0,maxLen):
			gl.canvas.create_line(x+gl.baseRectangleSymbolYPixelSize/2, y,x+gl.baseRectangleSymbolYPixelSize/2,y-(gl.baseRectangleSymbolYPixelSize), fill="black", width=1)
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
    refresh() 
	# Model.modelInstance.appendSequenceRecord(newSequenceRecord=MySeqRecord(seq=Seq(data="GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))

def denaturate( canvas:Canvas):
    for sequenceRecord in Model.modelInstance.sequenceRecordList:
        sequenceRecord:MySeqRecord
        if  sequenceRecord.hybridizedToStrand or sequenceRecord.hybridizedToPrimer:
            sequenceRecord.hybridizedToStrand=False
            sequenceRecord.hybridizedToPrimer=False
            sequenceRecord.singleStranded=True
    refresh() 

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
                # check is not blocked by another strand or the same primer
                if not strandRegularRecord.hybridizedToStrand and not strandRegularRecord.isPrimer and (not strandRegularRecord.hybridizedToPrimer or
                    (strandRegularRecord.hybridizedToPrimer and sequenceRecordPrimer.uniqueId!=strandRegularRecord.hybridizedToPrimer.uniqueId)): #and not strandRegularRecord.hybridizedToPrimer: # adrian avoid adding the same primer twice on same strand
                    if strandRegularRecord.fiveTo3: # <----
                        # print("Testing 5to3",strandRegularRecord.seq) 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedReversedPrimerSeq, minOverlapLength=minOverlapLength)
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecordPrimer.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
                            # change feature location
                            Model.modelInstance.sequenceRecordList[p].hybridizedToStrand=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedToPrimer=Model.modelInstance.sequenceRecordList[p]  
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=False
                            primRec.seq=Seq(seqToString(primRec.seq)[::-1])# reverses the string
                            Model.modelInstance.sequenceRecordList.insert(s+1,primRec)
                    else:  
                        # print("Testing 3to5",strandRegularRecord.seq)                 
                        overlaps, largestOverlapsInStrand, largestOverlapInPrimer =PrimerUtils.findPrimerOverlaps(targetDnaRecordSequence=strandRegularRecord.seq, primerRecordSequence=complementedPrimerSeq, minOverlapLength=minOverlapLength)                      
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)>1:
                            messagebox.showinfo("Problem", f"primer {sequenceRecordPrimer.seq} can bind {len (largestOverlapsInStrand)} times to strand {strandRegularRecord.seq} ") 
                            return
                        if largestOverlapsInStrand and len (largestOverlapsInStrand)==1:
                            found=True
                            sequenceRecordPrimer.xStartOffsetAsLetters=(largestOverlapsInStrand[0][0]-largestOverlapInPrimer[0][0])+strandRegularRecord.xStartOffsetAsLetters
                            Model.modelInstance.sequenceRecordList[p].hybridizedToStrand=Model.modelInstance.sequenceRecordList[s]
                            Model.modelInstance.sequenceRecordList[s].hybridizedToPrimer=Model.modelInstance.sequenceRecordList[p]
                            primRec:MySeqRecord=Model.modelInstance.sequenceRecordList.pop(p)
                            primRec.fiveTo3=True
                            Model.modelInstance.sequenceRecordList.insert(s,primRec)    
    if found:
        # gl.prefs.set_preference_value("shrink", False)
        refresh()                                           
    else:
        names:str=""
        for p, sequenceRecordPrimer in enumerate(Model.modelInstance.sequenceRecordList):
            if sequenceRecordPrimer.isPrimer and not sequenceRecordPrimer.hybridizedToStrand :
                 names+=(", "+sequenceRecordPrimer.description)
    
        messagebox.showinfo("No Anealing", f"Primers {names} do not anneal to any present sequence") 


def toggleShrink( canvas:Canvas):
    p=gl.shrink
    if p:
          gl.prefs.setPreferenceValue("shrink", False)
    else:
          gl.prefs.setPreferenceValue("shrink", True)
    refresh()            
              

def elongate( canvas:Canvas):
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
                                                                                   sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters].complement()+
                sequenceRecordPrimer.seq)
                
            seqRec=SeqRecord(subsequence, id=f"elongated {sequenceRecordPrimer.id}", name=f"from elongated primer {sequenceRecordPrimer.description}",
                                 description=f" {sequenceRecordPrimer.description}")
            newRec = MySeqRecord(seqRec, singleStranded=None,fiveTo3=sequenceRecordPrimer.fiveTo3,primer=False)
            
            featureLabel=f"seed primer "+sequenceRecordPrimer.description
            if sequenceRecordPrimer.fiveTo3:     
                newRec.xStartOffsetAsLetters=sequenceRecordPrimer.xStartOffsetAsLetters  
                oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(0, len(sequenceRecordPrimer.seq), strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
            else:
                oldPrimerFeature:SeqFeature=SeqFeature(SimpleLocation(sequenceRecordPrimer.xStartOffsetAsLetters, sequenceRecordPrimer.xStartOffsetAsLetters+len(sequenceRecordPrimer.seq), strand=None), type="old_sequence", id="elongated primer", qualifiers={"label": [featureLabel]})
                newRec.xStartOffsetAsLetters=sequenceRecordPrimer.hybridizedToStrand.xStartOffsetAsLetters     
            newRec.hybridizedToPrimer=False
            newRec.hybridizedToStrand=sequenceRecordPrimer.hybridizedToStrand
            newRec.uniqueId=sequenceRecordPrimer.uniqueId    
            newRec.features.insert(0,oldPrimerFeature)
            Model.modelInstance.sequenceRecordList.pop(i)
            Model.modelInstance.sequenceRecordList.insert(i, newRec)
            newRec.singleStranded=False
            found=True
    if not found:
        messagebox.showerror("Not found", "No primer ready to elongate") 
    refresh()            
              
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
    for i,r in enumerate(Model.modelInstance.sequenceRecordList):
        r:MySeqRecord
        if r.uniqueId ==mySeqRecord.uniqueId:
            Model.modelInstance.sequenceRecordList.pop(i)#deletion happens here
            gl.canvasLeft.delete("all")
            if r.hybridizedToPrimer:
                r.hybridizedToPrimer.hybridizedToStrand=None
            if r.hybridizedToStrand:
                r.hybridizedToStrand.hybridizedToPrimer=None                

            # print(f"the clicked sequence is found{r.uniqueId}")
            refresh()
            break

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

