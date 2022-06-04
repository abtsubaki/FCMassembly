"""                                                                                 
%prog to merge fasta files and only keep unique sequences                                                    
"""   


from Bio import SeqIO

with open('output.fa', 'a') as outFile:
    record_ids = list()
    for record in SeqIO.parse('input.fa', 'fasta'):
        if record.id not in record_ids:
            record_ids.append( record.id)
            SeqIO.write(record, outFile, 'fasta')

