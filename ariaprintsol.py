import tkinter as tk
from tkinter import Canvas, Scrollbar, Frame, RIGHT, Y, BOTTOM, TOP, LEFT, X, BOTH, NW, SE, NONE

class ScrollableCanvas(Frame):
    def __init__(self, parent, width=500, height=500):
        super().__init__(parent)

        # Create a canvas with a white background
        self.canvas = Canvas(self, width=width, height=height, bg="white")

        # Create horizontal and vertical scrollbars
        self.hbar = Scrollbar(self, orient='horizontal', command=self.canvas.xview)
        self.vbar = Scrollbar(self, orient='vertical', command=self.canvas.yview)

        # Configure the canvas to use the scrollbars
        self.canvas.config(xscrollcommand=self.hbar.set, yscrollcommand=self.vbar.set)

        # Pack the canvas and scrollbars
        self.canvas.pack(side=LEFT, fill=BOTH, expand=True)
        self.hbar.pack(side=BOTTOM, fill=X)
        self.vbar.pack(side=RIGHT, fill=Y)

        # Bind the canvas to the scrollbars
        self.canvas.bind("<Configure>", self._configure_canvas)

    def _configure_canvas(self, event):
        # Resize the canvas scroll region to fit the content
        self.canvas.config(scrollregion=self.canvas.bbox("all"))

    def draw_circle(self, x, y, radius):
        # Draw a blue circle on the canvas
        self.canvas.create_oval(x - radius, y - radius, x + radius, y + radius, fill="blue")

    def print_canvas(self, filename):
        # Create a new canvas with the same dimensions as the original
        print_canvas = Canvas(width=self.canvas.winfo_width(), height=self.canvas.winfo_height(), bg="white")

        # Draw all the items from the original canvas onto the new canvas
        for item in self.canvas.find_all():
            coords = self.canvas.coords(item)
            print_canvas.create_oval(coords, fill=self.canvas.itemcget(item, "fill"))

        # Save the new canvas as a PostScript file
        print_canvas.postscript(file=filename, colormode='color')

root = tk.Tk()
root.title("Scrollable Canvas")

canvas_frame = ScrollableCanvas(root)
canvas_frame.pack(fill=BOTH, expand=True)

# Draw a circle that extends outside the visible area
canvas_frame.draw_circle(500, 250, 100)

# Print the canvas content
canvas_frame.print_canvas("canvas_print.ps")

root.mainloop()

