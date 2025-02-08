from  tkinter import Canvas
from enhancedbutton import EnhancedButton
from preferences import *

prefs:Preferences=None
canvasLeft:Canvas=None
canvas:Canvas=None
#
canvasLeftPadding=None # blank space from the left edge of canvas
canvas=None # blank space from the left edge of canvas
mask:list[int]=None
maskSkipped:list[int]=None
canvasHeight=800
debug=True
# cache below
fontName=None
fontSize=None
canvasLeftPadding=None
horizontalPixelsMargin=None
baseRectangleSymbolXPixelSize=None
baseRectangleSymbolYPixelSize=None
shrink:int=None
coloredBases:bool=None
rotated:bool=None
verticalSequenceSpacing:int=None
hydrogen:bool=None
hydrogenLinesHalfLength:int=None
leftButtonsWidth=None
minPrimerOverlapLength=None
a=None
c=None
g=None
t=None