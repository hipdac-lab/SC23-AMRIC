import numpy as np
import sys

if len(sys.argv) < 2:
    print("Usage: python script_name.py <file_name>")
    sys.exit(1)

file_name = sys.argv[1]

total_pre_time = 0.0
total_real_write_time = 0.0

iAvg = 0
outList = np.zeros(400)

with open(file_name, 'r') as file:
    data = file.readlines()

for line in data:
    if line.startswith("Write h5plotfile time ="):
        iAvg += 1
        number = float(line.split("Write h5plotfile time =")[1].strip().split(" ")[0])  # Extract the number and remove "seconds"
        print(f"AMRIC-SZ_Interp Total time = {number} seconds")

for line in data:
    if line.startswith("real write time ="):
        iAvg += 1
        number = float(line.split("real write time =")[1].strip().split(" ")[0])  # Extract the number and remove "seconds"
        total_real_write_time += number
for line in data:
    if line.startswith("pre time ="):
        iAvg += 1
        number = float(line.split("pre time =")[1].strip().split(" ")[0])  # Extract the number and remove "seconds"
        total_pre_time += number

print(f"AMRIC-SZ_Interp Preprocess time = {total_pre_time} seconds")
print(f"AMRIC-SZ_Interp Writing+Compression time = {total_real_write_time} seconds")
