import tkinter as tk
from util import *
from gl import *
from buttoncommands import *
from preferences import Preferences

class UiApp:

    def __init__(self, root:tk):
        self.canvas=None
        self.root:tk = root
        self.root.title("OpenBio")
        screen_width = root.winfo_screenwidth()-20
        screen_height = int(root.winfo_screenheight()/2)
        self.root.geometry(f"{screen_width}x{screen_height}+0+0")
        gl.prefs=Preferences(self.root, self.refresh)
        self.createLeftCanvasForEnhancedBottons()
        self.createCanvasWithScrollbar()  
        # self.root.grid_columnconfigure(0, weight=0)    
        self.root.grid_columnconfigure(1, weight=1)    # make the main canvas take all the horizontal space
        Preferences.updateGlobalCache()
        self.createButtonBars()
        self.createMenu()
        # Configure grid to allow resizing  
        self.root.grid_rowconfigure(0, weight=1)
        self.root.grid_columnconfigure(1, weight=1)

    def createLeftCanvasForEnhancedBottons(self):
        self.canvasLeft = tk.Canvas(root, width=gl.prefs.getPreferenceValue("leftButtonsWidth"), height=gl.canvasHeight, bg='lightblue')   
        self.canvasLeft.grid(row=0, column=0, padx=0, pady=0, sticky="nw")
        gl.canvasLeft=self.canvasLeft

    def createCanvasWithScrollbar(self):
        # Create a frame to hold the canvas and the scrollbar
        canvasFrame = tk.Frame(self.root)
        canvasFrame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")
        # Create a canvas
        self.canvas = tk.Canvas(canvasFrame, bg="lightgray", height=gl.canvasHeight)
        
        # Create a horizontal scrollbar
        self.hScrollbar = tk.Scrollbar(canvasFrame, orient="horizontal", command=self.canvas.xview)
        self.canvas.config(xscrollcommand=self.hScrollbar.set)
        # Create a vertical scrollbar        
        self.vScrollbar = tk.Scrollbar(canvasFrame, orient="vertical", command=self.canvas.yview)
        self.canvas.config(yscrollcommand=self.vScrollbar.set)
        # Pack the canvas and scrollbar
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.hScrollbar.grid(row=1, column=0, sticky="ew")
        self.vScrollbar.grid(row=0, column=1, sticky="ns")
        # Configure grid to allow resizing
        canvasFrame.grid_columnconfigure(0, weight=1)
        canvasFrame.grid_rowconfigure(0, weight=2)
        gl.canvas=self.canvas

    def createButtonBars(self):
        # Create horizontal button bar at the bottom
        self.bottomButtonBar = tk.Frame(self.root)
        self.bottomButtonBar.grid(row=1, column=1, padx=5, pady=0, sticky="ew")
        column=0
        #load sequences
        loadFileButton=tk.Button(self.bottomButtonBar, text="Load Sequences", command=lambda:self.loadSequencesHandler())
        loadFileButton.grid(row=0, column=column, padx=3)
        column+=1
        #add new sequences
        addFileButton=tk.Button(self.bottomButtonBar, text="Add New Sequences", command=lambda:self.addSequencesHandler())
        addFileButton.grid(row=0, column=column, padx=3)
        column+=1     
        # add primer   
        addPrimerButton=tk.Button(self.bottomButtonBar, text="Add Primer", command=lambda:addPrimerHandler(self.canvas))
        addPrimerButton.grid(row=0, column=column, padx=3)
        column+=1  
        # Create zoom in button      
        buttonZi = tk.Button(self.bottomButtonBar, text="Zoom In", command=lambda: self.canvasZoom(True))
        buttonZi.grid(row=0, column=column, padx=0)
        column+=1
        # Create zoom out button
        buttonZo = tk.Button(self.bottomButtonBar, text="Zoom Out", command=lambda: self.canvasZoom(False))
        buttonZo.grid(row=0, column=column, padx=3)
        # print
        column+=1           
        buttonPrint = tk.Button(self.bottomButtonBar, text="Print", command=lambda: self.printCanvas(self.canvas))
        buttonPrint.grid(row=0, column=column, padx=3)
        column+=1
        #shrink
        buttonShrink = tk.Button(self.bottomButtonBar, text="Toggle Shrink", command=lambda:toggleShrink(self.canvas))
        buttonShrink.grid(row=0, column=column, padx=3)   
        #refresh
        column+=1        
        buttonRefresh = tk.Button(self.bottomButtonBar, text="Refresh", command=self.refresh)
        buttonRefresh.grid(row=0, column=column, padx=3)   
        #debug
        column+=1 
        def toggleDebug():
            gl.debug = not gl.debug  
            self.refresh()  
        buttonDebug = tk.Button(self.bottomButtonBar, text="Toggle Debug", command=toggleDebug)
        buttonDebug.grid(row=0, column=column, padx=20)   
        #
        column+=1        
        buttonXXX = tk.Button(self.bottomButtonBar, text="Button", command=None)
        buttonXXX.grid(row=0, column=column, padx=20)   
        # next row with buttons
        self.stepsButtonBar = tk.Frame(self.root)
        self.stepsButtonBar.grid(row=2, column=1, padx=5, pady=3, sticky=( "ew")) 
        # denaturate
        column=0
        buttonDenaturate = tk.Button(self.stepsButtonBar, text="Denaturate", command=lambda:denaturate(self.canvas))
        buttonDenaturate.grid(row=0, column=column)      
        column+=1   
        # aneal
        buttonAneal = tk.Button(self.stepsButtonBar, text="Aneal", command=lambda:anealPrimers(self.canvas))
        buttonAneal.grid(row=0, column=column)   
        column+=1   
        #Elongate
        buttonElongate = tk.Button(self.stepsButtonBar, text="Elongate", command=lambda:elongate(self.canvas))
        buttonElongate.grid(row=0, column=column)   


    def exitApplication(self):
        self.root.quit()      

    def createMenu(self):
        menuBar = Menu(root)
        # Left menu: File Load
        fileMenu = Menu(menuBar, tearoff=0)
        fileMenu.add_command(label="Load File", command=lambda:self.loadSequencesHandler())
        menuBar.add_cascade(label="File", menu=fileMenu)
        # menu: Preferences
        preferencesMenu = Menu(menuBar, tearoff=0)
        preferencesMenu.add_command(label="Preferences", command=self.updatePreferences)
        menuBar.add_cascade(label="Preferences", menu=preferencesMenu)
        # Right menu: Exit
        exitMenu = Menu(menuBar, tearoff=0)
        exitMenu.add_command(label="Exit", command=self.exitApplication)
        menuBar.add_cascade(label="Exit", menu=exitMenu)
        # Configuring the root window to use the menu
        root.config(menu=menuBar)   

    def updatePreferences(self):
        gl.prefs.openPreferencesPopup()       

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
        oldSize=gl.prefs.getPreferenceValue(name)
        if zoomin!=True and oldSize>3: # ZOOM OUT no to negative font sizes
                gl.prefs.getPreferenceValue(name)
                gl.prefs.setPreferenceValue(name, oldSize-1)
                drawCanvas(self.canvas)
        else:   
            gl.prefs.setPreferenceValue(name,oldSize+1 )           
            drawCanvas(self.canvas)


if __name__ == "__main__":
    root = tk.Tk()
    app = UiApp(root)
    loadModel(default=True)
    app.root.title("OpenBio "+Model.modelInstance.loadedFileName)
    drawCanvas(app.canvas)
    root.mainloop()
