import tkinter as tk
from tkinter import Canvas, Button
from myseqrecord import MySeqRecord
import gl


class EnhancedButton:
    def __init__(self, canvas, text, x, y, mySeqRecord:MySeqRecord, labelHeightPx):
        # Initialize button properties
        self.canvas = canvas
        self.text = text
        self.x = x
        self.y = y
        self.labelWidthPx = gl.canvasHorizontalMargin    
        self.mySeqRecord=mySeqRecord    
        # Create the label
        label = canvas.create_rectangle(x , y,   x + self.labelWidthPx , y + labelHeightPx,  fill="lightblue", outline="black", width=1)
        # Create text inside the rectangle to represent the label's text
        text=canvas.create_text(x+self.labelWidthPx//2, y+labelHeightPx//2, text=text, font=("Arial", 12), fill="black")
        from buttoncommands import clickOnSeqRecord
        canvas.tag_bind(label, "<Button-1>", lambda event: clickOnSeqRecord(event, canvas, self.mySeqRecord))
        canvas.tag_bind(text, "<Button-1>",  lambda event: clickOnSeqRecord(event, canvas, self.mySeqRecord))


def main():
    root = tk.Tk()
    root.title("Label as Button Inside Canvas")

    # Create the Canvas widget
    canvas = tk.Canvas(root, bg="lightgray", width=500, height=gl.canvasHeight)
    canvas.pack(padx=10, pady=10)
    # Sample drawing on the canvas
    canvas.create_rectangle(10, 50, 60, 150, fill="blue")
    # Create a red point at (20, 30)
    def action1(event: tk.Event) -> None:
        eb.handle_click(event)
    # Create labels using the createLabel method
    eb:EnhancedButton=EnhancedButton(canvas, "Action 1", 20, 50, None, action1, 12)
    canvas.create_line(0, 50, 200, 50, fill="red", width=1)  # Draw a 1-pixel point
    # Run the Tkinter event loop
    root.mainloop()


# if __name__ == "__main__":
#     main()