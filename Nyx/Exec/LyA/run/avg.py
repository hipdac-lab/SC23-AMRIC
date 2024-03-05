import numpy as np

def calculate_average_from_bin(file_path):
    # Load the binary file as a numpy array of type float64
    data = np.fromfile(file_path, dtype=np.float64)
    
    # Calculate the average
    average = np.mean(data)
    
    return average

# Example usage
file_path = 'gold-241.raw'
average = calculate_average_from_bin(file_path)
print(f"The average of the binary file is: {average}")