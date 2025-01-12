from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from pathlib import Path  

def write_string_to_pdf(pdf_path, text):
    """Writes a string to a PDF file.

    Args:
        pdf_path (str): The path to the PDF file to be created.
        text (str): The string to write to the PDF.
    """

    # Create a PDF canvas
    c = canvas.Canvas(pdf_path, pagesize=letter)

    # Set font and size
    c.setFont("Helvetica", 12)  # You can choose a different font if desired

    # Write the text to the PDF
    c.drawString(100, 700, text)  # Adjust position as needed

    # Save the PDF
    c.showPage()
    c.save()
    print("PDF created with string.")

# Specify the path to the PDF file
pdf_path = filePath=str(Path(__file__).resolve().parent)+"/output.pdf"

s="fffffffffffffffffffffffffffffffffffffAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA END "

# Write the string to the PDF
write_string_to_pdf(pdf_path, s)
