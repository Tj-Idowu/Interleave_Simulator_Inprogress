import random
import os

# Function to make, scramble, and output scrambled number list
def generate_and_scramble_sequence(start, end, output_filename="fastq_number_order.txt"):
    """
    Generates a list of numbers from 'start' to 'end' (inclusive), 
    scrambles the list, and writes the scrambled sequence to a file.
    """
    
    # 1. Generate the list of numbers (START_NUM to END_NUM)

    print(f"Generating list of numbers from {start} to {end}...")
    number_list = list(range(start, end + 1))
    
    # 2. Scramble the numbers
    print("Scrambling the order...")
    random.shuffle(number_list)
    
    # 3. Print the scrambled order to a file
    print(f"Writing scrambled order to '{output_filename}'...")
    
    # Open the file for writing
    try:
        with open(output_filename, 'w') as f:
            # Join the list of integers into a single string with newlines
            f.write('\n'.join(map(str, number_list)))
            
        print(f"Process complete! The file '{output_filename}' has been created.")
        print(f"File path: {os.path.abspath(output_filename)}")

    except IOError as e:
        print(f"Error writing to file: {e}")

# --- Execution ---
# Define the range (1 to 1,000,000)
START_NUM = 1
END_NUM = 40000000

# Run the function
generate_and_scramble_sequence(START_NUM, END_NUM)
