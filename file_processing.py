"""
Convert Fasta files to SAM
"""
import pyranges as pr
exons, cpg = pr.data.exons(), pr.data.cpg()

#select chromosome, start point and endpoint
exons["chrY", "-",  15591259:27197945]
