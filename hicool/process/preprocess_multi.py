import os
from multiprocessing import Pool, cpu_count
import gzip
import re
import threading

def demultiplex_fastq(file1, file2, output_dir, barcodes):
    # Create output file objects for each cell barcode
    file_handles = {}
    for barcode in barcodes:
        file_name = os.path.join(output_dir, barcode + ".fastq.gz")
        file_handles[barcode] = gzip.open(file_name, "at")
    
    # Create a lock object for output file objects
    lock = threading.Lock()
    
    # Open input files
    with gzip.open(file1, "rt") as f1, gzip.open(file2, "rt") as f2:
        while True:
            # Read 4 lines at a time from input files
            try:
                seqid1, seq1, qualid1, qual1 = [f1.readline().strip() for _ in range(4)]
                seqid2, seq2, qualid2, qual2 = [f2.readline().strip() for _ in range(4)]
            except:
                break
            
            # Extract cell barcode from read 1
            barcode_match = re.search(r"(\+[ACGT]+)([ACGT]+)", qualid1)
            if not barcode_match:
                continue
            barcode = barcode_match.group(2)
            
            # Write reads to the corresponding output file object
            if barcode in barcodes:
                with lock:
                    file_handles[barcode].write(seqid1 + "\n" + seq1 + "\n" + qualid1 + "\n" + qual1 + "\n")
                    file_handles[barcode].write(seqid2 + "\n" + seq2 + "\n" + qualid2 + "\n" + qual2 + "\n")
            
   
