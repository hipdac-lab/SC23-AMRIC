import numpy as np

def read_and_add_bin_files(file1_path, file2_path, output_file_path):
    # Load the binary files as numpy arrays
    data1 = np.fromfile(file1_path, dtype=np.float64)
    data2 = np.fromfile(file2_path, dtype=np.float64)
    
    # Check if the files have the same size
    if data1.size != data2.size:
        raise ValueError("The binary files do not have the same size.")
    
    # Check for positions where both are zero or both are non-zero
    both_zero = np.logical_and(data1 == 0, data2 == 0)
    both_non_zero = np.logical_and(data1 != 0, data2 != 0)
    if np.any(both_zero) or np.any(both_non_zero):
        print("Found positions where both data points are zero or both are non-zero.")
    
    # Add the contents of the two files
    result = data1 + data2
    
    # Save the result to a new binary file
    result.tofile(output_file_path)

# Example usage
file1_path = 'test.raw'
file2_path = 'test-f.raw'
output_file_path = 'test-a.raw'
try:
    read_and_add_bin_files(file1_path, file2_path, output_file_path)
    print(f"Result saved to {output_file_path}")
except ValueError as e:
    print(f"Error: {e}")