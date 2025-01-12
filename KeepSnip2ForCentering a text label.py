import tkinter as tk
from tkinter import font as tkFont

def center_text(text_widget, canvas_width, text):
    # Get the current font of the Text widget
    current_font = tkFont.Font(font=text_widget.cget("font"))
    # Calculate the width of the text using the current font
    text_width = current_font.measure(text)
    print("canvas_width",canvas_width)
    print("text_width",text_width)
    # Calculate the number of spaces needed for centering
    num_spaces = (canvas_width - text_width) // (2 * current_font.measure(' '))
    spaces = ' ' * max(num_spaces, 0)  # Ensure non-negative space count
    print("spaces",max(num_spaces, 0))
    # Clear the current text and insert centered text
    text_widget.delete("1.0", tk.END)
    text_widget.insert("1.0", spaces + text)

# Create the main application window
root = tk.Tk()
root.title("Canvas with Centered Text Widget Example")

# Create a Canvas
windowWidth = 1600
canvas = tk.Canvas(root, width=windowWidth, height=200, bg='lightgrey')
canvas.pack()

# Create a Text widget
text_widget = tk.Text(canvas, width=20, height=5, wrap='word', font=("Arial", 12), bg='white', padx=5, pady=5)

# Insert the Text widget into the Canvas at position (0, 0)
canvas.create_window(0, 0,width=windowWidth, anchor='nw', window=text_widget)

# Insert centered text
center_text(text_widget,windowWidth, "Cent")

# Start the Tkinter main loop
root.mainloop()


#  current_font = tkFont.Font(font=text_widget.cget("font"))
    
#     # Calculate the width of the text using the current font
#     text_width = current_font.measure(text)
    
#     # Get the width of the Text widget
#     canvas_width = text_widget.winfo_width()