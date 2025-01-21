import tkinter as tk
from tkinter import Canvas, Button
import gl


class EnhancedButton:
    def __init__(self, canvas, text, x, y, action, labelHeightPx):
        # Initialize button properties
        self.canvas = canvas
        self.text = text
        self.x = x
        self.y = y
        self.action = action
        self.labelWidthPx = gl.canvasHorizontalMargin        
        # Create the label
        label = canvas.create_rectangle(x , y,   x + self.labelWidthPx , y + labelHeightPx,  fill="lightblue", outline="black", width=1)
        # Create text inside the rectangle to represent the label's text
        text=canvas.create_text(x+self.labelWidthPx//2, y+labelHeightPx//2, text=text, font=("Arial", 12), fill="black")
        canvas.tag_bind(label, "<Button-1>", self.handleClick)
        canvas.tag_bind(text, "<Button-1>", self.handleClick)

    def handleClick(self, event: tk.Event) -> None:
        # Get the coordinates of the click
        x, y = event.x, event.y
        # Find the bounding box of the text
        bbox = self.canvas.bbox("current")
        item = self.canvas.find_closest(x, y)  # Find the closest item to the mouse click
        if item and self.canvas.type(item)=="text" : # Get the type of the item    
            text = self.canvas.itemcget("current", "text")
            # If the click is within the bounding box, calculate the letter clicked
            if bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3]:
                # Find the character index in the text that was clicked
                char_index = int((x - bbox[0]) / (bbox[2] - bbox[0]) * len(text))
                clicked_char = text[char_index]
                print(f"Action 1 triggered: {text}, Clicked character: '{clicked_char}'")
        else:
                print(f"Action 1 was not over a character")


def main():
    root = tk.Tk()
    root.title("Label as Button Inside Canvas")

    # Create the Canvas widget
    canvas = tk.Canvas(root, bg="lightgray", width=500, height=400)
    canvas.pack(padx=10, pady=10)
    # Sample drawing on the canvas
    canvas.create_rectangle(10, 50, 60, 150, fill="blue")
    # Create a red point at (20, 30)
    def action1(event: tk.Event) -> None:
        eb.handle_click(event)
    # Create labels using the createLabel method
    eb:EnhancedButton=EnhancedButton(canvas, "Action 1", 20, 50, action1, 12)
    canvas.create_line(0, 50, 200, 50, fill="red", width=1)  # Draw a 1-pixel point
    # Run the Tkinter event loop
    root.mainloop()


# if __name__ == "__main__":
#     main()