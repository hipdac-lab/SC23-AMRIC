import numpy as np

def calculate_errors_and_values(file1_path, file2_path):
    # Load the binary files as numpy arrays
    data1 = np.fromfile(file1_path, dtype=np.float64)
    data2 = np.fromfile(file2_path, dtype=np.float64)
    
    # Ensure both files have the same number of elements
    if data1.size != data2.size:
        raise ValueError("Files have different sizes and cannot be compared.")
    
    # Skip positions where either data1 or data2 is zero for absolute error calculation
    valid_indices_abs = np.logical_and(data1 != 0, data2 != 0)
    abs_errors = np.abs(data1[valid_indices_abs] - data2[valid_indices_abs])
    if abs_errors.size == 0:
        max_abs_error = 0
        max_abs_error_values = (0, 0)
    else:
        max_abs_error_index = np.argmax(abs_errors)
        max_abs_error = abs_errors[max_abs_error_index]
        # Find the original index in the filtered array
        original_index_abs = np.nonzero(valid_indices_abs)[0][max_abs_error_index]
        max_abs_error_values = (data1[original_index_abs], data2[original_index_abs])
    
    # Calculate the point-wise absolute relative error, skipping positions where data2 is 0
    valid_indices_rel = data1 != 0  # Indices where data2 is not zero
    relative_errors = np.zeros(data1.shape)
    relative_errors[valid_indices_rel] = np.abs((data1 - data2) / data1)[valid_indices_rel]
    max_relative_error_index = np.argmax(relative_errors)
    max_relative_error = relative_errors[max_relative_error_index]
    max_relative_error_values = (data1[max_relative_error_index], data2[max_relative_error_index])
    
    return max_abs_error, max_abs_error_values, max_relative_error, max_relative_error_values

# Example usage
file1_path = 'test.raw'
file2_path = 'gold-241.raw'
try:
    max_abs_error, max_abs_error_values, max_relative_error, max_relative_error_values = calculate_errors_and_values(file1_path, file2_path)
    print(f"Maximum abs Error: {max_abs_error}, Values: {max_abs_error_values}")
    print(f"Maximum rel error: {max_relative_error}, Values: {max_relative_error_values}")
except ValueError as e:
    print(f"Error: {e}")
