import os
import sys

def main(file_a, file_b):
    if not os.path.isfile(file_a) or not os.path.isfile(file_b):
        print("Error: One or both files do not exist.")
        sys.exit(1)

    size_a = os.path.getsize(file_a)
    size_b = os.path.getsize(file_b)

    ratio = size_a / size_b
    print(f"CR is: {ratio:.2f}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python file_size_ratio.py <file_a> <file_b>")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

