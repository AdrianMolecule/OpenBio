import tkinter as tk
from util import *

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Tkinter Layout Example")

        # Create the canvas with a horizontal scrollbar
        self.createCanvasWithScrollbar()
        self.createButtonBars()

    def createCanvasWithScrollbar(self):
        # Create a frame to hold the canvas and the scrollbar
        canvasFrame = tk.Frame(self.root)
        canvasFrame.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")

        # Create a canvas
        self.canvas = tk.Canvas(canvasFrame, bg="lightgray", height=400)
        
        # Create a horizontal scrollbar
        self.hScrollbar = tk.Scrollbar(canvasFrame, orient="horizontal", command=self.canvas.xview)
        self.canvas.config(xscrollcommand=self.hScrollbar.set)

        # Pack the canvas and scrollbar
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.hScrollbar.grid(row=1, column=0, sticky="ew")

        # Add a 500-character long string to the canvas
        long_string = "This is a long string " * 25  # Creates a string of 500 characters
        self.canvas.create_text(10, 20, anchor="nw", text=long_string, font=("Arial", 12))

        # Add some content to the canvas for testing the scrollbar
        self.canvas.create_rectangle(0, 0, 1200, 400, fill="blue")

        # Configure grid to allow resizing
        canvasFrame.grid_columnconfigure(0, weight=1)
        canvasFrame.grid_rowconfigure(0, weight=1)

    def createButtonBars(self):
        # Create horizontal button bar at the bottom
        self.bottomButtonBar = tk.Frame(self.root)
        self.bottomButtonBar.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

        buttons = ["Button 1", "Button 2", "Button 3"]
        for i, button in enumerate(buttons):
            tk.Button(self.bottomButtonBar, text=button).grid(row=0, column=i, padx=5)

        # Create vertical button bar on the right
        self.rightButtonBar = tk.Frame(self.root)
        self.rightButtonBar.grid(row=0, column=1, padx=10, pady=10, sticky="ns")

        verticalButtons = ["Button A", "Button B", "Button C"]
        for i, button in enumerate(verticalButtons):
            tk.Button(self.rightButtonBar, text=button).grid(row=i, column=0, pady=5)

        # Configure grid to allow resizing
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)


def loadModel():
    seqRecList:list[SeqRecord]=loadFile()
    print ("Initial load of model", seqRecList )
    model.sequenceRecordList=seqRecList
    model.dumpModel("in main")
    model.appendSequenceRecord(SeqRecord(Seq("GATATAT"),id="AdrianShortSeq", name="AdrianSecondSeqName"))    

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    loadModel()
    canvasDraw(app.canvas)
    root.mainloop()
