import tkinter as tk
from myseqrecord import MySeqRecord
import gl

class EnhancedButton:
    def __init__(self, text, x, y, mySeqRecord:MySeqRecord, labelHeightPx):
        # Initialize button properties
        y+=2 # to deal with some alignments between the 2 canvases
        # self.text = text 
        self.mySeqRecord=mySeqRecord    
        self.popupWindow = None
        # Create text inside the rectangle to represent the label's text
        labelDelete = gl.canvasLeft.create_rectangle(x , y+labelHeightPx-7,   x + gl.leftButtonsWidth , y+labelHeightPx,  fill="red", outline="black", width=1)
        from util import clickOnSeqRecordToDelete
        gl.canvasLeft.tag_bind(labelDelete, "<Button-1>", lambda event: clickOnSeqRecordToDelete(event,  self.mySeqRecord))

        gl.canvasLeft.create_rectangle(x , y,   x + gl.leftButtonsWidth , y + labelHeightPx-8,  fill="lightblue", outline="black", width=1)
        textDescription_id=gl.canvasLeft.create_text(x+gl.leftButtonsWidth//2, y+labelHeightPx-9, text=text, font=("Arial", 8), angle=90, fill="black", anchor="w")
        gl.canvasLeft.tag_bind(textDescription_id, "<Enter>", lambda event, text_id=textDescription_id: on_hover(event))  # Mouse enters
        gl.canvasLeft.tag_bind(textDescription_id, "<Leave>", lambda event, text_id=textDescription_id: on_leave(event))  # Mouse leaves
        gl.canvasLeft.tag_bind(textDescription_id, "<Button-1>", lambda event, text_id=textDescription_id: on_click(event))  # Mouse leaves
    
        def on_hover(event):
            showInfo(event)        
    
        def on_click(event):    
            showInfo(event)
        
        def showInfo(event):
            global popupWindow
            if self.popupWindow:
                return # don't open more popups
            textId = event.widget.find_withtag("current")[0]
            gl.canvasLeft.itemconfig(textDescription_id, fill="yellow")  # Change the text color to red
            current_text = event.widget.itemcget(textId, "text")            
            popupWindow = tk.Toplevel()
            popupWindow.geometry("1200x600")
            popupWindow.title(f"Info on {self.mySeqRecord.description} sequence")                      
            # Display the text in the popup window
            textBox = tk.Text(popupWindow, wrap=tk.NONE, width=40, height=10)            
            textBox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            scrollbarV = tk.Scrollbar(popupWindow, orient=tk.VERTICAL, command=textBox.yview)            
            scrollbarV.pack(side=tk.RIGHT, fill=tk.Y)
            scrollbarH = tk.Scrollbar(textBox, orient=tk.HORIZONTAL, command=textBox.xview)#for some reasons while the vertical can have master  is popupwindow this does not 
            scrollbarH.pack(side=tk.BOTTOM, fill=tk.X)            
            # textBox.config(yscrollcommand=scrollbarV.set, xscrollcommand=scrollbarH.set)            
            textBox.insert(tk.END, self.mySeqRecord.__str__())
            # label = tk.Label(popupWindow, text=self.mySeqRecord.__str__(), font=("Helvetica", 14),  anchor="nw", justify="left")
            # label.pack(fill=tk.BOTH, expand=True)
            #popupWindow.after(1000, popupWindow.destroy)  # Close the popup after 1 second
            return popupWindow  # Return the reference to the popup window                        
        
        # Function to reset the text color when mouse leaves
        def on_leave(event):
            text_id = event.widget.find_withtag("current")[0]
            gl.canvasLeft.itemconfig(text_id, fill="black")  # Change the text color to red   
            global popupWindow
            if popupWindow is not None and popupWindow.winfo_exists():
                popupWindow.destroy()  # Close the pop-up         
                            
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
    eb:EnhancedButton=EnhancedButton( "Action 1", 20, 50, None, action1, 10)
    # Run the Tkinter event loop
    root.mainloop()


# if __name__ == "__main__":
#     main()