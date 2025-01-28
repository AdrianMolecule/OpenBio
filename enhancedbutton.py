import tkinter as tk
from tkinter import Canvas, Button
from myseqrecord import MySeqRecord
import gl


class EnhancedButton:
    def __init__(self, canvas:Canvas, text, x, y, mySeqRecord:MySeqRecord, labelHeightPx):
        # Initialize button properties
        self.canvas = canvas
        self.text = text
        self.x = x
        self.y = y
        self.labelWidthPx = gl.canvasHorizontalMargin    
        self.mySeqRecord=mySeqRecord    
        from buttoncommands import clickOnSeqRecordToDelete
        labelDelete = canvas.create_rectangle(x , y+labelHeightPx-8,   x + self.labelWidthPx , y+labelHeightPx,  fill="red", outline="black", width=1)
        canvas.tag_bind(labelDelete, "<Button-1>", lambda event: clickOnSeqRecordToDelete(event, canvas, self.mySeqRecord))
        popupWindow = None
        def on_hover(event):
            global popupWindow
            textId = event.widget.find_withtag("current")[0]
            canvas.itemconfig(textDescription_id, fill="yellow")  # Change the text color to red
            current_text = event.widget.itemcget(textId, "text")
            popupWindow = tk.Toplevel()
            popupWindow.title("Info on selected sequence")            
            # Display the text in the popup window
            label = tk.Label(popupWindow, text=self.mySeqRecord.description, font=("Helvetica", 14))
            label.pack(padx=20, pady=20)   
            popupWindow.after(1000, popupWindow.destroy)  # Close the popup after 1 second
            return popupWindow  # Return the reference to the popup window                        

        # Function to reset the text color when mouse leaves
        def on_leave(event):
            text_id = event.widget.find_withtag("current")[0]
            canvas.itemconfig(text_id, fill="black")  # Change the text color to red   
            if popupWindow is not None and popupWindow.winfo_exists():
                print("destroy pop")
                popupWindow.destroy()  # Close the pop-up         
                            
        textDescription_id=canvas.create_text(x+self.labelWidthPx//2, y+labelHeightPx//2, text=text, font=("Arial", 8), angle=90, fill="black")
        canvas.tag_bind(textDescription_id, "<Enter>", lambda event, text_id=textDescription_id: on_hover(event))  # Mouse enters
        canvas.tag_bind(textDescription_id, "<Leave>", lambda event, text_id=textDescription_id: on_leave(event))  # Mouse leaves


def main():
    root = tk.Tk()
    root.title("Label as Button Inside Canvas")

    # Create the Canvas widget
    canvas = tk.Canvas(root, bg="lightgray", width=500, height=gl.canvasHeight)
    canvas.pack(padx=10, pady=10)
    # Sample drawing on the canvas
    canvas.create_rectangle(10, 50, 60, 150, fill="blue")
    canvas.create_rectangle(10, 50, 60, 150, fill="blue")
    # Create a red point at (20, 30)
    def action1(event: tk.Event) -> None:
        eb.handle_click(event)
    # Create labels using the createLabel method
    eb:EnhancedButton=EnhancedButton(canvas, "Action 1", 20, 50, None, action1, 10)
    canvas.create_line(0, 50, 200, 50, fill="red", width=1)  # Draw a 1-pixel point
    # Run the Tkinter event loop
    root.mainloop()


# if __name__ == "__main__":
#     main()