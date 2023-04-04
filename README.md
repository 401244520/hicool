# hicool
HiCool
HiCool is a Python package that provides an object-oriented interface for working with HiCool and SCool data files. This package allows you to load, manipulate, and save HiCool and SCool files with ease.

Installation
pip install hicool
# Usage
from hicool import HiCool

# Load a HiCool or SCool file
hc = HiCool("path/to/file.hicool")

# Save the HiCool object as a HiCool file
hc.save_as("path/to/output.hicool")

# Get metadata information
hc.info()

# Convert the HiCool object to a Scanpy AnnData object
sce = hc.to_scanpy(embedding_name="my_embedding")

# Convert the HiCool object to a Higashi object
hig = hc.to_higashi()

# Convert the HiCool object to a fast Higashi object
fast_hig = hc.to_fast_higashi()

# compress_matrix
compress_matrix is a function that compresses a HiCool object into a 2D matrix. Each row in the matrix represents a cell, and each column represents a compressed feature.

hc = HiCool("path/to/file.hicool")
matrix, features = compress_matrix(hc, chrom="chr1", operation="bin_strength")
