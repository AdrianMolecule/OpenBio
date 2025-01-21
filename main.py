import tkinter as tk
from util import *
from buttoncommands import *
from preferences import Preferences

class UiApp:

    def __init__(self, root):
        self.canvas=None
        self.root = root
        self.root.title("OpenBio")
        screen_width = root.winfo_screenwidth()-20
        screen_height = int(root.winfo_screenheight()/3)
        self.root.geometry(f"{screen_width}x{screen_height}+0+0")
        # Create the canvas with a horizontal scrollbar
        self.createCanvasWithScrollbar()
        gl.prefs=Preferences(self.root, self.refresh)
        self.createButtonBars()
        self.createMenu()

    def createCanvasWithScrollbar(self):
        # Create a frame to hold the canvas and the scrollbar
        canvasFrame = tk.Frame(self.root)
        canvasFrame.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
        # Create a canvas
        self.canvas = tk.Canvas(canvasFrame, bg="lightgray", height=450)
        # Create a horizontal scrollbar
        self.hScrollbar = tk.Scrollbar(canvasFrame, orient="horizontal", command=self.canvas.xview)
        self.canvas.config(xscrollcommand=self.hScrollbar.set)
        # Create a vertical scrollbar        
        self.vScrollbar = tk.Scrollbar(canvasFrame, orient="vertical", command=self.canvas.yview)
        self.canvas.config(xscrollcommand=self.hScrollbar.set)
        self.canvas.config(yscrollcommand=self.vScrollbar.set)
        # Pack the canvas and scrollbar
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.hScrollbar.grid(row=1, column=0, sticky="ew")
        self.vScrollbar.grid(row=0, column=1, sticky="ns")
        # Configure grid to allow resizing
        canvasFrame.grid_columnconfigure(0, weight=1)
        canvasFrame.grid_rowconfigure(0, weight=1)

    def createButtonBars(self):
        # Create horizontal button bar at the bottom
        self.bottomButtonBar = tk.Frame(self.root)
        column=0
        self.bottomButtonBar.grid(row=1, column=column, padx=5, pady=5, sticky="ew")
        column+=1
        loadFileButton=tk.Button(self.bottomButtonBar, text="Load Sequences", command=lambda:self.loadSequencesHandler())
        loadFileButton.grid(row=0, column=column, padx=5)
        column+=1
        addFileButton=tk.Button(self.bottomButtonBar, text="Add New Sequences", command=lambda:self.addSequencesHandler())
        addFileButton.grid(row=0, column=column, padx=5)
        column+=1        
        addPrimerButton=tk.Button(self.bottomButtonBar, text="Add Primer", command=lambda:addPrimerHandler(self.canvas))
        addPrimerButton.grid(row=0, column=column, padx=5)
        column+=1        
        buttonZi = tk.Button(self.bottomButtonBar, text="Zoom In", command=lambda: self.canvasZoom(True))
        buttonZi.grid(row=0, column=column, padx=10, pady=10)
        column+=1
        # Create zoom out button
        buttonZo = tk.Button(self.bottomButtonBar, text="Zoom Out", command=lambda: self.canvasZoom(False))
        buttonZo.grid(row=0, column=column, padx=10, pady=10)
        column+=1        
        # Create print button
        buttonPrint = tk.Button(self.bottomButtonBar, text="Print", command=lambda: self.printCanvas(self.canvas))
        buttonPrint.grid(row=0, column=column, padx=10, pady=10)
        column+=1   
        buttonRefresh = tk.Button(self.bottomButtonBar, text="Refresh", command=self.refresh)
        column+=1   
        buttonRefresh.grid(row=0, column=0, pady=5)
        #
        # Create Steps button bar
        def incrementRowOrCol(row,col,vertical):
            if vertical:
                return row+1,col
            else:
                return row,col+1
        
        self.stepsButtonBar = tk.Frame(self.root)
        vertSteps=gl.prefs.get_preference_value("verticalSteps")
        row=0 if vertSteps else 2
        col=1

        self.stepsButtonBar.grid(row=row, column=col if vertSteps else "0", padx=5, pady=5, sticky=("ns" if vertSteps else "ew")) 
        
        buttonDenaturate = tk.Button(self.stepsButtonBar, text="Denaturate", command=lambda:denaturate(self.canvas))
        buttonDenaturate.grid(row=row, column=col, pady=5)      

        row,col=incrementRowOrCol(row,col,vertSteps)
        buttonAneal = tk.Button(self.stepsButtonBar, text="Aneal", command=lambda:anealPrimers(self.canvas))
        buttonAneal.grid(row=row, column=col, pady=5)   

        row,col=incrementRowOrCol(row, col, vertSteps)
        buttonElongate = tk.Button(self.stepsButtonBar, text="Elongate", command=lambda:elongate())
        buttonElongate.grid(row=row, column=col, pady=5)   
        row,col=incrementRowOrCol(row, col, vertSteps)
        
        #
        # Configure grid to allow resizing
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(0, weight=1)

    def exitApplication(self):
        self.root.quit()      

    def createMenu(self):
        menuBar = Menu(root)
        # Left menu: File Load
        fileMenu = Menu(menuBar, tearoff=0)
        fileMenu.add_command(label="Load File", command=lambda:self.loadSequencesHandler())
        menuBar.add_cascade(label="File", menu=fileMenu)
        # Right menu: Exit
        exitMenu = Menu(menuBar, tearoff=0)
        exitMenu.add_command(label="Exit", command=self.exitApplication)
        menuBar.add_cascade(label="Exit", menu=exitMenu)
        # menu: Preferences
        preferencesMenu = Menu(menuBar, tearoff=0)
        preferencesMenu.add_command(label="Preferences", command=self.updatePreferences)
        menuBar.add_cascade(label="Preferences", menu=preferencesMenu)
        # Configuring the root window to use the menu
        root.config(menu=menuBar)   

    def updatePreferences(self):
        gl.prefs.open_preferences_popup()       

    def loadSequencesHandler(self)->Seq:
        loadModel(False, append=False)
        self.root.title("OpenBio "+Model.modelInstance.loadedFileName)
        drawCanvas(self.canvas) 

    def addSequencesHandler(self)->Seq:
        loadModel(False, append=True)
        self.root.title("OpenBio "+Model.modelInstance.loadedFileName)
        drawCanvas(self.canvas) 

    def refresh(self):
        drawCanvas(self.canvas) 

    def canvasDrawCircle(self):
        self.canvas.create_oval(100, 150, 200, 250, outline="blue", width=2)

    def canvasZoom(self,zoomin):# 1 for  zoom In or bigger
        name="fontSize"
        oldSize=gl.prefs.get_preference_value(name)
        if zoomin!=True and oldSize>3: # ZOOM OUT no to negative font sizes
                gl.prefs.get_preference_value(name)
                gl.prefs.set_preference_value(name, oldSize-1)
                drawCanvas(self.canvas)
        else:   
            gl.prefs.set_preference_value(name,oldSize+1 )           
            drawCanvas(self.canvas)


if __name__ == "__main__":
    root = tk.Tk()
    app = UiApp(root)
    loadModel(default=True)
    app.root.title("OpenBio "+Model.modelInstance.loadedFileName)
    drawCanvas(app.canvas)
    root.mainloop()
