def dna_transcription(dna_sequences):
    """Convert a DNA sequence into an mRNA sequence by replacing Thymine (T) with Uracil (U)."""
    return dna_sequences.replace('T', 'U') 
print(dna_transcription("ATCGAAACTTTGCAAGGCTT"))

def translate_rna(rna_sequence):
    # Define the codon table with three-letter amino acid codes using dictionary
    Rna_codon = {
        'AUA':'Ile', 'AUC':'Ile', 'AUU':'Ile', 'AUG':'Met',
        'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACU':'Thr',
        'AAC':'Asn', 'AAU':'Asn', 'AAA':'Lys', 'AAG':'Lys',
        'AGC':'Ser', 'AGU':'Ser', 'AGA':'Arg', 'AGG':'Arg',
        'CUA':'Leu', 'CUC':'Leu', 'CUG':'Leu', 'CUU':'Leu',
        'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCU':'Pro',
        'CAC':'His', 'CAU':'His', 'CAA':'Gln', 'CAG':'Gln',
        'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGU':'Arg',
        'GUA':'Val', 'GUC':'Val', 'GUG':'Val', 'GUU':'Val',
        'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCU':'Ala',
        'GAC':'Asp', 'GAU':'Asp', 'GAA':'Glu', 'GAG':'Glu',
        'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGU':'Gly',
        'UCA':'Ser', 'UCC':'Ser', 'UCG':'Ser', 'UCU':'Ser',
        'UUC':'Phe', 'UUU':'Phe', 'UUA':'Leu', 'UUG':'Leu',
        'UAC':'Tyr', 'UAU':'Tyr', 'UAA':'Stop', 'UAG':'Stop',
        'UGC':'Cys', 'UGU':'Cys', 'UGA':'Stop', 'UGG':'Trp'
    }
  
   # Process the RNA sequence in groups of 3 (codons)
protein = ""   # this define the protein outside the for loop
rna_sequence = ["AUCGAAACUUUGCAAGGCUU"]
for i in range(0, len(rna_sequence), 3):
    codon = rna_sequence[i:i+3] # for the sliced sequence
    amino_acid = Rna_codon.get(codon, '?')  # Get the corresponding amino acid
    protein += RNA_codon[codon]
    if amino_acid == 'Stop':  # Stop translation if stop codon is found
            break
    protein += amino_acid  # Append amino acid to protein sequence
    
print(protein)