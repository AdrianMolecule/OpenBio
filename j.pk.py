import tkinter as tk

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("Text Box Copy Example")

        # Create StringVar for text boxes
        self.text1 = tk.StringVar()
        self.text2 = tk.StringVar()
        self.text3 = tk.StringVar()
        self.text4 = tk.StringVar()

        # Create text boxes
        self.createTextBoxes()

        # Bind changes of text1 to copy to text2
        self.text1.trace("w", self.copyToTextBox2)

        # Bind changes of text3 to copy to text4
        self.text3.trace("w", self.copyToTextBox4)

        # Bind Enter key in text2 to call printTextBox1

        # Create Reset and Print buttons
        self.createResetButton()
        self.createPrintButton()

    def createTextBoxes(self):
        # Create and pack 4 text boxes with labels
        tk.Label(self.root, text="Text Box 1:").grid(row=0, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text1).grid(row=0, column=1, padx=5, pady=5)

        tk.Label(self.root, text="Text Box 2:").grid(row=1, column=0, padx=5, pady=5)
        self.text2_entry = tk.Entry(self.root, textvariable=self.text2)
        self.text2_entry.bind("<Return>", self.printTextBox1_on_enter)
        self.text2_entry.grid(row=1, column=1, padx=5, pady=5)

        tk.Label(self.root, text="Text Box 3:").grid(row=2, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text3).grid(row=2, column=1, padx=5, pady=5)

        tk.Label(self.root, text="Text Box 4:").grid(row=3, column=0, padx=5, pady=5)
        tk.Entry(self.root, textvariable=self.text4).grid(row=3, column=1, padx=5, pady=5)

    def copyToTextBox2(self, *args):
        # Copy the value from text1 to text2
        self.text2.set(self.text1.get())

    def copyToTextBox4(self, *args):
        # Copy the value from text3 to text4
        self.text4.set(self.text3.get())

    def createResetButton(self):
        # Create a reset button
        reset_button = tk.Button(self.root, text="Reset", command=self.resetAllTextBoxes)
        reset_button.grid(row=4, column=0, columnspan=2, pady=10)

    def createPrintButton(self):
        # Create a print button to print the value of text1
        print_button = tk.Button(self.root, text="Print Text Box 1", command=self.printTextBox1)
        print_button.grid(row=5, column=0, columnspan=2, pady=10)

    def resetAllTextBoxes(self):
        # Reset all text boxes to 'xxx'
        self.text1.set("xxx")
        self.text2.set("xxx")
        self.text3.set("xxx")
        self.text4.set("xxx")

    def printTextBox1(self):
        # Print the value of text1 to the console
        print("Value of Text Box 1:", self.text1.get())

    def printTextBox1_on_enter(self, event):
        # Call the printTextBox1 method when Enter is pressed in text2
        self.printTextBox1()

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
