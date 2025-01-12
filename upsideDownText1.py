import tkinter as tk
from tkinter import Canvas

class UpsideDownText(tk.Text):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Create a canvas to draw rotated text
        self.canvas = Canvas(self, bg=self["bg"], highlightthickness=0)
        self.canvas.pack(fill=tk.BOTH, expand=True)

    def _draw_rotated_text(self, text, x, y):
        # Rotate the text by 180 degrees and place it on the canvas
        self.canvas.delete("all")  # Clear previous drawings
        self.canvas.create_text(
            x, y, text=text, font=self.cget("font"), anchor=tk.CENTER, angle=180
        )

    def insert(self, index, text, *args):
        super().insert(index, text, *args)  # Insert the text into the Text widget
        # Draw the inserted text rotated
        self._draw_rotated_text(self.get("1.0", tk.END).strip(), 100, 50)

    def delete(self, index1, index2=None):
        super().delete(index1, index2)  # Delete in the Text widget
        # Redraw the text after deletion
        self._draw_rotated_text(self.get("1.0", tk.END).strip(), 100, 50)

# Create the main application window
root = tk.Tk()
root.title("Upside Down Text Example")

# Create an instance of the UpsideDownText class
# upside_down_text = UpsideDownText(root, wrap='word', width=400, height=200, font=("Arial", 14))

upside_down_text = UpsideDownText(root, width=10, height=5, wrap=tk.WORD, bd=1, relief="solid", highlightbackground="red")
upside_down_text.insert(tk.END, "HELLO WORLD AGAIN AND AGAIN")  # Insert the text

upside_down_text.pack(padx=10, pady=10)

# Insert some initial text to demonstrate functionality
upside_down_text.insert("1.0", "Hello, World!\nThis text is displayed upside down!")

# Start the Tkinter main loop
root.mainloop()
