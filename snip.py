import tkinter as tk

# Create the main window
root = tk.Tk()
root.title("Two Canvases Example")

# Create the first canvas (200 pixels wide)
canvas1 = tk.Canvas(root, width=50, height=100, bg='lightblue')
canvas1.grid(row=0, column=0, padx=10, pady=10)

# Add a label to the first canvas
label1 = tk.Label(canvas1, text="This is Canvas 1", bg='lightblue')
label1.pack(pady=20)  # Centering the label vertically

# Create the second canvas that fills the remaining space
canvas2 = tk.Canvas(root, bg='lightgreen')
canvas2.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)

# Add a label to the second canvas
label2 = tk.Label(canvas2, text="This is Canvas 2", bg='lightgreen')
label2.pack(pady=20)  # Centering the label vertically

# Configure the grid to allow the second column to expand
root.grid_columnconfigure(1, weight=1)

# Start the Tkinter event loop
root.mainloop()
