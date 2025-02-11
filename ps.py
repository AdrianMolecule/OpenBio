import random
import string
import os


def generate_postscript(filename, text):
    """
    Generates a PostScript file containing the given text.

    Args:
        filename (str): The name of the output PostScript file.
        text (str): The text to be written to the file.
    """

    with open(filename, 'w') as file:
        file.write("%!PS-Adobe-3.0 EPSF-3.0\n")
        file.write("%%BoundingBox: 0 0 1000 100\n")  # Adjust bounding box as needed
        file.write("/Courier findfont 12 scalefont setfont\n")
        file.write("0 50 moveto\n")
        text="START",text,"END"
        file.write(f"({text}) show\n")
        # Draw a red circle
        file.write("500 50 50 0 360 arc\n")  # Circle centered at (500, 50) with radius 50
        file.write("0 setgray\n")  # Set fill color to black
        file.write("fill\n")

        file.write("showpage\n")

# Generate a random string of 1000 letters
random_string = "START aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa   END"

# Generate the PostScript file
#filePath=str(Path(__file__).resolve().parent)+"/"+
generate_postscript(os.getcwd()+"/"+"long_string.ps", random_string)
print(f"PostScript file 'long_string.ps' generated with a string of 1000 letters.")
