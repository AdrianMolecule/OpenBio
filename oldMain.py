import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import Menu
from tkinter import Canvas, Label, Scrollbar, Frame, RIGHT, Y, BOTTOM, X, HORIZONTAL, LEFT, BOTH
#
from pathlib import Path 
#  
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas 
#
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
# mine
from util import *
import model
import gl 

#
def loadFileCommandHandler()->Seq:
    model.sequenceRecordList=loadFile()
    # Build the UI components
    canvasDraw()

def loadFile()->list[SeqRecord]:
    """Prompts the user to load a .embl file, reads its contents, and draws the string on the canvas."""
    defaultTestFileValue:int=gl.prefs.get_preference_value(preference_name="defaultTestFileValue") 
    if defaultTestFileValue is None or defaultTestFileValue=="":
        filePath: str = filedialog.askopenfilename(title="Open EMBL File", filetypes=[("EMBL Files", "*.embl"), ("All Files", "*.*")])    
    else:
        filePath=str(Path(__file__).resolve().parent)+DEFAULT_FILE
    if filePath:
        try:
            secRecList:list[SeqRecord]=loadEmblSequences(filePath)
            print(secRecList,secRecList.__class__)                
            return secRecList
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred while reading the file: {e}")
    else:
            messagebox.showwarning("No file", "Please select a file")


def exitApplication():
    """Exits the application."""
    root.quit()

def createLowerButtons(buttonFrame:Frame):
    buttonLoad = tk.Button(buttonFrame, text="Load File", command=loadFileCommandHandler)
    buttonLoad.grid(row=0, column=0, padx=10, pady=10)
    # Create additional buttons for drawing shapes
    buttonCircle = tk.Button(buttonFrame, text="Draw Circle", command=lambda: canvasDrawCircle())
    buttonCircle.grid(row=0, column=1, padx=10, pady=10)
    # Create zoom In button
    buttonZi = tk.Button(buttonFrame, text="Zoom In", command=lambda: canvasZoom(True))
    buttonZi.grid(row=0, column=2, padx=10, pady=10)
    # Create zoom out button
    buttonZo = tk.Button(buttonFrame, text="Zoom Out", command=lambda: canvasZoom(False))
    buttonZo.grid(row=0, column=3, padx=10, pady=10)
    # Create print button
    buttonPrint = tk.Button(buttonFrame, text="Print", command=lambda: printCanvas(canvas))
    buttonPrint.grid(row=0, column=4, padx=10, pady=10)

def createRightButtons(buttonFrame:Frame):
    # Create a button to load a file
    buttonLoad = tk.Button(buttonFrame, text="Load File", command=loadFileCommandHandler)
    buttonLoad.grid(row=0, column=0, padx=10, pady=10)
    # Create additional buttons for drawing shapes
    buttonCircle = tk.Button(buttonFrame, text="Draw Circle", command=lambda: canvasDrawCircle())
    buttonCircle.grid(row=0, column=1, padx=10, pady=10)
    # Create zoom In button
    buttonZi = tk.Button(buttonFrame, text="Zoom In", command=lambda: canvasZoom(True))
    buttonZi.grid(row=0, column=2, padx=10, pady=10)
    # Create zoom out button
    buttonZo = tk.Button(buttonFrame, text="Zoom Out", command=lambda: canvasZoom(False))
    buttonZo.grid(row=0, column=3, padx=10, pady=10)
    # Create print button
    buttonPrint = tk.Button(buttonFrame, text="Print", command=lambda: printCanvas(canvas))
    buttonPrint.grid(row=0, column=4, padx=10, pady=10)

def createMenu():
# Create menus
    menuBar = Menu(root)
    # Left menu: File Load
    fileMenu = Menu(menuBar, tearoff=0)
    fileMenu.add_command(label="Load File", command=loadFile)
    menuBar.add_cascade(label="File", menu=fileMenu)
    # Right menu: Exit
    exitMenu = Menu(menuBar, tearoff=0)
    exitMenu.add_command(label="Exit", command=exitApplication)
    menuBar.add_cascade(label="Exit", menu=exitMenu)
    # Configuring the root window to use the menu
    root.config(menu=menuBar)         

def createCanvas(CanvasFrame:Frame):
    global canvas
    # Create a canvas to draw on
    canvas = Canvas(CanvasFrame, bg="white", width=800, height=600, scrollregion=(0, 0, 800, 600))
    # Create a vertical scrollbar for the canvas
    vScroll = Scrollbar(CanvasFrame, orient='vertical')
    vScroll.pack(side=RIGHT, fill=Y)
    vScroll.config(command=canvas.yview)
    # Create a horizontal scrollbar for the canvas
    hScroll = Scrollbar(CanvasFrame, orient=HORIZONTAL)
    hScroll.pack(side=BOTTOM, fill=X)
    hScroll.config(command=canvas.xview)
    # Configure the canvas to use the scrollbars
    canvas.config(yscrollcommand=vScroll.set, xscrollcommand=hScroll.set)
    canvas.pack(side=LEFT, expand=True, fill=BOTH)
    return canvas

def buildUI():
    """Builds the UI components including canvas, buttons, and menus."""
    canvasFrame = tk.Frame(root)
    canvasFrame.pack(side=tk.BOTTOM, fill=tk.X)
    createCanvas(canvasFrame)
    bottomFrame = tk.Frame(root)
    bottomFrame.pack(side=tk.BOTTOM, fill=tk.X)
    createLowerButtons(bottomFrame)
    rightFrame = tk.Frame(root)
    rightFrame.pack(side=tk.RIGHT, fill=tk.Y)
    createRightButtons(rightFrame)
    createMenu()
    
def loadModel():
    seqRecList:list[SeqRecord]=loadFile()
    print ("Initial load of model", seqRecList )
    model.sequenceRecordList=seqRecList
    model.dumpModel("in main")
    model.appendSequenceRecord(SeqRecord(Seq("GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))    

def main():
    global root
    root = tk.Tk()
    """Main method to set up and run the Tkinter application."""
    root.title("Resizable Canvas with Buttons")
    screen_width = root.winfo_screenwidth()-20
    screen_height = int(root.winfo_screenheight()/3)
    root.geometry(f"{screen_width}x{screen_height}+0+0")
    # Maximize the window on startup
    #root.state('zoomed')  # This maximizes the window on most systems
    loadModel()
    # Build the UI components
    buildUI()
    canvasDraw(canvas)
    # Start the Tkinter event loop
    root.mainloop()
    return


# Run the main method to start the application
if __name__ == "__main__":
    main()
