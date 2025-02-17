import tkinter as tk
import math

# Set up the main window
root = tk.Tk()
root.title("Logarithmic Scale Rectangles")

# Create a canvas widget
canvas = tk.Canvas(root, width=600, height=500)
canvas.pack()

# Heights to plot (in logarithmic scale)
heights = [100, 200, 300, 400, 500, 650]

logBase = 10# Define the base for the logarithmic scale (e.g., base 10)
maxCanvasHeight = 1000
maxLogHeight = math.log(max(heights), logBase)# Find the maximum logarithmic height for scaling
scaleFactor = maxCanvasHeight / maxLogHeight# Define a scaling factor to fit the maximum log-transformed height to maxCanvasHeight
max_log_height = math.log(max(heights), logBase) * scaleFactor# Find the maximum logarithmic height to normalize the scaling
# Draw each rectangle along the X-axis, transformed by log scale on Y-axis
for i, height in enumerate(heights):
    # Logarithmic transformation for the height
    logHeight = math.log(height, logBase) * scaleFactor
    x1=20
    y1 = maxCanvasHeight - logHeight  # Y position, subtracted to flip the log height (higher values go down)
    canvas.create_line(x1, y1, x1+100, y1)

# Run the tkinter event loop
root.mainloop()
