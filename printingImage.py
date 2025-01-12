import tkinter as tk
from tkinter import filedialog
from PIL import ImageGrab

def printCanvas(canvas):
    # Get the canvas's total scrollable area (scrollregion)
    print("in printCanvas")
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
        #img = ImageGrab.grab(bbox=(capture_x1, capture_y1, capture_x2, capture_y2))
        img=ImageGrab.grab(all_screens=True)
        # Ask the user for a file name and save the image
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
        if file_path:
            img.save(file_path)
            print(f"Canvas content saved to {file_path}")
    else:
        print("Canvas has no scrollable content.")

# Example Tkinter window with a scrollable canvas
def buildUI():
    window = tk.Tk()
    window.title("Tkinter Print Scrollable Canvas Example")
    window.geometry("600x400")  # Set window size

    # Create a canvas frame 
    canvasFrame = tk.Frame(window)
    canvasFrame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    # Create a canvas inside a frame to hold the scrollbars
    global canvas  # Making canvas a global variable to use it across different functions
    canvas = tk.Canvas(canvasFrame, bg="lightgray", yscrollincrement=10, xscrollincrement=10)

    # Create vertical and horizontal scrollbars for the canvas
    verticalScrollbar = tk.Scrollbar(canvasFrame, orient="vertical", command=canvas.yview)
    horizontalScrollbar = tk.Scrollbar(canvasFrame, orient="horizontal", command=canvas.xview)
    canvas.config(yscrollcommand=verticalScrollbar.set, xscrollcommand=horizontalScrollbar.set)

    # Place the scrollbars around the canvas
    canvas.grid(row=0, column=0, sticky="nsew")
    verticalScrollbar.grid(row=0, column=1, sticky="ns")
    horizontalScrollbar.grid(row=1, column=0, sticky="ew")

    # Configure the grid to expand the canvas as the window resizes
    canvasFrame.grid_rowconfigure(0, weight=1)
    canvasFrame.grid_columnconfigure(0, weight=1)

    # Draw some content on the canvas
    canvas.create_rectangle(50, 50, 200, 150, fill="yellow", outline="black")
    canvas.create_oval(900, 50, 450, 150, fill="green", outline="black")
    canvas.create_text(250, 400, text="Canvas Content", font=("Arial", 16))
    
    # Add extra content to make it scrollable
    for i in range(100, 800, 50):
        canvas.create_line(i, 0, i, 800, fill="gray", dash=(4, 4))
    
    # Update the scroll region after drawing content
    canvas.config(scrollregion=canvas.bbox("all"))

    # Add a button to trigger the print canvas method
    print_button = tk.Button(window, text="Print Canvas", command=lambda: printCanvas(canvas))
    print_button.pack(pady=20)

    window.mainloop()

# Run the Tkinter UI
buildUI()
