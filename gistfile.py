#!/usr/bin/env python

"""                                                                                 
%prog to extract sequences from a fasta file based on name or sequence identifier                                                   
"""
                                                                                
from Bio import SeqIO                                                               
import sys                                                                          
                                                                                    
wanted = [line.strip() for line in open(sys.argv[2])]                               
seqiter = SeqIO.parse(open(sys.argv[1]), 'fasta')                                    
SeqIO.write((seq for seq in seqiter if seq.id in wanted), sys.stdout, "fasta")
