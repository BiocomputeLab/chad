#!/usr/bin/env python
# ========================================================
# Codon Harmonization And Diversification (CHAD) algorithm
# Authors: Shivang Joshi <shivang.joshi@bristol.ac.uk>,
#          Thomas E. Gorochowski <tom@chofski.co.uk>
# ========================================================


import random
import csv
import sys


# Load the codon data for the organism being used
def load_codon_data (filename):
    aa_codons = {}
    codon_freq = {}
    codon_rare = {}
    with open(filename, 'r') as file:
        csvreader = csv.reader(file)
        header = next(csvreader)
        for row in csvreader:
            if len(row) == 4:
                codon = row[0].strip()
                aa = row[1].strip()
                freq = float(row[2].strip())
                rare = int(row[3].strip())
                if aa not in aa_codons:
                    aa_codons[aa] = [codon]
                    codon_freq[aa] = [freq]
                else:
                    aa_codons[aa].append(codon)
                    codon_freq[aa].append(freq)
                if codon not in codon_freq:
                    codon_rare[codon] = [rare]
                else:
                    codon_rare[codon].append(rare)
            else:
                print("Malformed CSV file. Requires 4 columns per row.")
    return aa_codons, codon_freq, codon_rare


# Load the amino acid seqence list that needs requires diverse DNA sequences
def load_aa_seqs (filename):
    aa_seqs = {}
    cur_name = ""
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if len(line) > 0:
                if line[0] == ">":
                    cur_name = line[1:]
                else:
                    aa_seqs[cur_name] = line
    return aa_seqs


# Check is the current sequence is different enough from others generated
def check_seq_diverse (cur_seq, ex_seqs, max_homology):
    for idx in range(len(cur_seq)-(max_homology-1)):
        cur_seq_part = cur_seq[idx:idx+(max_homology-1)]
        for name, ex_seq in ex_seqs.items():
            if cur_seq_part in ex_seq:
                return False
    return True


# Check that current sequence is compatible with cloning enzymes
def check_cloning_compatible (cur_seq):
    enzymes = {"BsaI":["GGTCTC",  "GAGACC"],
               "BbsI":["GAAGAC",  "GTCTTC"],
               "SapI":["GCTCTTC", "GAAGAGC"]}
    for name, seq in enzymes.items():
        if seq[0] in cur_seq or seq[1] in cur_seq:
            return False
    return True


# Check that rare codons are not clustered or used.
def check_rare_codons_ok (cur_seq, codon_rare):
    # TODO: Rare codon checks not used at present
    return True


# Randomly sample a possible DNA sequence given an amino acid sequence and codon frequencies
def sample_dna_seq (aa_seq, aa_codons, codon_freq):
    dna_seq = ""
    for aa_idx in range(len(aa_seq)):
        aa = str(aa_seq[aa_idx])
        # `choices' returns a list so need to take first element
        codon = random.choices(aa_codons[aa], weights=codon_freq[aa], k=1)[0]
        dna_seq += str(codon)
    return dna_seq


# Cycle through amino acid sequences and generate diverse set of corrisponding DNA sequences
def generate_diverse_seqs (aa_seqs, aa_codons, codon_freq, codon_rare, max_homology):
    print("Codon Harmonization And Diversification (CHAD) algorithm")
    dna_seqs = {}
    dna_seq = ""
    for aa_name, aa_seq in aa_seqs.items():
        MAX_TRIES = 500000
        try_num = 1
        while try_num < MAX_TRIES:
            dna_seq = sample_dna_seq(aa_seq, aa_codons, codon_freq)
            if (check_seq_diverse(dna_seq, dna_seqs, max_homology) and
                check_cloning_compatible(dna_seq) and
                check_rare_codons_ok(dna_seq, codon_rare)):
                break
            try_num += 1
        if try_num >= MAX_TRIES:
            print("Could not generate compatible sample for " + str(aa_name))
        else:
            dna_seqs[aa_name] = dna_seq
            print("Generated sequence for: " + aa_name)
    return dna_seqs


# Write dict of name -> seq pairs to FASTA format file
def write_to_fasta (seqs, filename):
    with open(filename, 'w') as file:
        for name, seq in seqs.items():
            file.write(">" + str(name) + "\n")
            file.write(str(seq) + "\n\n")


# Print the usage of the command
def print_usage ():
    print("Codon Harmonization And Diversification (CHAD) algorithm")
    print("Usage: chad MAX_HOMOLOGY_BP CODON_DATA_CSV INPUT_AA_FASTA OUTPUT_DNA_FASTA")


# Handle command line args and generate diversified sequences
def main():
    # Ensure results are reproducible so use a fixed seed
    random.seed(123)
    # Handle command line arguments
    if len(sys.argv) == 2 and (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print_usage()
        sys.exit()
    elif len(sys.argv) != 5:
        print("CHAD requires 4 arguments. Type `chad -h` for help")
        sys.exit(1)
    else:
        max_homology = int(sys.argv[1])
        codon_data_filename = sys.argv[2]
        aa_fasta_filename = sys.argv[3]
        output_dna_fasta_filename = sys.argv[4]
        aa_codons, codon_freq, codon_rare = load_codon_data(codon_data_filename)
        aa_seqs = load_aa_seqs(aa_fasta_filename)
        dna_seqs = generate_diverse_seqs(aa_seqs, aa_codons, codon_freq, codon_rare, max_homology)
        write_to_fasta(dna_seqs, output_dna_fasta_filename)


if __name__ == "__main__":
    main()
