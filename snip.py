import random

# Create a list to store the multiplications
multiplications = []
multiplicationsWithAnswers = []

multPerPage=520
pages=1
timeInMinutes=round(multPerPage*pages*5/60,0)
print(f"This should take about {timeInMinutes}" + " minutes at 5 seconds per multiplication\n")
for _ in range(multPerPage*pages):
    num1 = random.randint(10, 99)  # Random two-digit number
    num2 = random.randint(10, 99)  # Random two-digit number
    result=(str(num1*num2).rjust(4))
    multiplication = f"{num1}x{num2}=        "  # Store only the multiplication with blank spaces after "="
    multiplications.append(multiplication)
    multiplicationWithAnswer = f"{num1}x{num2}={result}"  # Store only the multiplication with blank spaces after "="
    multiplicationsWithAnswers.append(multiplicationWithAnswer)
    cols=8
# Write the multiplications to a file called 'mul.txt'
mulString=(f"this should take about {timeInMinutes}" + " minutes\n")
for i in range(0, len(multiplications), cols):
    # Write 4 multiplications per line, separated by commas
    mulString+=(",\t".join(multiplications[i:i+cols]) + "\n")
mulStringWithAnswers=(f"this should take about {timeInMinutes}" + " minutes\n")
for i in range(0, len(multiplicationsWithAnswers), cols):
    # Write 4 multiplications per line, separated by commas
    mulStringWithAnswers+=(",\t".join(multiplicationsWithAnswers[i:i+cols]) + "\n")
print(mulString)
print(mulStringWithAnswers)
# Write the multiplications to a file called 'mul.txt'
with open("mul.txt", "w") as file:
    file.write(mulString)
with open("mulAnswers.txt", "w") as file:
    file.write(mulStringWithAnswers)
print("Multiplications have been written to 'mul.txt'.")    
