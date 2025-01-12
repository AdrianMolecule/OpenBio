import tkinter as tk

# Create the main application window
root = tk.Tk()
root.title("Scrolling Text Example with Border")

# Create a frame to hold the Text widget and the Scrollbar
frame = tk.Frame(root)
frame.pack(pady=20)

# Create a Text widget with a red border
text_widget = tk.Text(frame, width=10, height=5, wrap=tk.WORD, bd=1, relief="solid", highlightbackground="red")
text_widget.insert(tk.END, "HELLO WORLD AGAIN AND AGAIN")  # Insert the text
#text_widget.config(state=tk.DISABLED)  # Make it read-only

text_widget.pack()

# Start the Tkinter main loop
root.mainloop()
