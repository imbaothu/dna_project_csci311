# Implement edit distance algorithm

""" Edit distance - this is another interesting approach used in practice
(and also widely used in the field of computational linguistics). Given two sequences,
s and t, how many “edit operations” are required to transform one sequence into the 
other? The set of edit operations you would consider are insertion, deletion and substitution.
See https://en.wikipedia.org/wiki/Edit_distance for more information. """

import dna_project

def editDistance(queryfile, sequencesfile):
    # read the files and chop it up by the sequences
    query = dna_project.read_fasta(queryfile)
    sequences = dna_project.read_fasta(sequencesfile)

    # insert/delete/substitute for each sequence and count the operations performed
    print(sequences)
    # return sequence with the smallest number of operations

    return

editDistance("DNA_query.txt", "DNA_sequences.txt")