import tkinter as tk
from tkinter import messagebox
import random
points=0
class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Math game")

        # Create StringVar for text boxes
        self.text1 = tk.StringVar()
        self.text2 = tk.StringVar()
        self.text3 = tk.StringVar()
        self.text4 = tk.StringVar()

        # Create text boxes
        self.createTextBoxes()
        # Create a Reset button
        self.createResetButton()
        self.createCheckButton()

    def createTextBoxes(self):
        # Create and pack 4 text boxes with labels
        tk.Label(self.root, text="First number").grid(row=0, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text1).grid(row=0, column=1, padx=5, pady=5)

        tk.Label(self.root, text="Second number").grid(row=1, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text2).grid(row=1, column=1, padx=5, pady=5)

        tk.Label(self.root, text="Guess").grid(row=2, column=0, padx=5, pady=5)
        self.text2_entry=tk.Entry(self.root, textvariable=self.text3)
        self.text2_entry.grid(row=2, column=1, padx=5, pady=5)
        self.text2_entry.bind("<Return>", self.check_on_enter)

        tk.Label(self.root, text="Points").grid(row=3, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text4).grid(row=3, column=1, padx=5, pady=5)
        self.text1.set(str((random.randint(2, 9))))
        self.text2.set(str((random.randint(2, 9))))
        self.text4.set(points)

    def copyToTextBox2(self, *args):
        # Copy the value from text1 to text2
        self.text2.set(self.text1.get())

    def copyToTextBox4(self, *args):
        # Copy the value from text3 to text4
        self.text4.set(self.text3.get())

    def createResetButton(self):
        # Create a reset button
        reset_button = tk.Button(self.root, text="Restart", command=self.resetAllTextBoxes)
        reset_button.grid(row=4, column=0, columnspan=2, pady=10)

    def createCheckButton(self):
        # Create a reset button
        reset_button = tk.Button(self.root, text="Check", command=self.check)
        reset_button.grid(row=4, column=1, columnspan=2, pady=10)

    def resetAllTextBoxes(self):
        # Reset all text boxes to 'xxx'
        self.text1.set(str((random.randint(2, 9))))
        self.text2.set(str((random.randint(2, 9))))
        self.text4.set(points)


    def check(self):
        correct=False
        if int(self.text3.get() )==int(self.text1.get())*int(self.text2.get()):
            print("Correct")
            correct=True
            self.showPopup("Correct")
            self.points = points + 1
        else:
            print("Incorrect")
            self.showPopup("Incorrect")
            self.points = points -1


    def showPopup(self, message):
        # Show a message box with a 1-second delay
        messagebox.showinfo("Info", message)


    def check_on_enter(self, event):
        self.check()

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
