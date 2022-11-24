import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import Align
from IPython.display import display

DNA_bases = ['A', 'C', 'G', 'T']


def get_triplicate_length(given_length):
    neareast_divisible_by_3 = int(given_length/3)
    nearest_3_multiple = neareast_divisible_by_3 * 3 
    return nearest_3_multiple





def randomly_generate_DNA_sequences(DNA_min_length, DNA_seq_max_length):

    #randomly generate a number between 1 and 200 
    DNA_length = np.random.randint(DNA_min_length, DNA_seq_max_length)
    DNA_length = get_triplicate_length(DNA_length)
    DNA_sequence = []
    for base_pair in range(DNA_length):
        base_pair = np.random.choice(DNA_bases)
        #print(base_pair)
        DNA_sequence.append(base_pair)
    return DNA_sequence


def create_mulitple_DNA_sequences(number_DNA_sequences, min_DNA_length, max_DNA_length):
    list_of_DNA_sequences= []
    for sequence in range(number_DNA_sequences):
        DNA_sequence = randomly_generate_DNA_sequences(min_DNA_length, max_DNA_length)
        list_of_DNA_sequences.append(DNA_sequence)

    return(list_of_DNA_sequences)


def replace_thymidine_with_uracil(list_of_DNA_sequences):
    for DNA_sequence in list_of_DNA_sequences:
        for base in range(len(DNA_sequence)):
            if DNA_sequence[base] ==  'T':
                DNA_sequence[base] =  'U'
    return list_of_DNA_sequences

def transcribe_with_DNA_bases(list_of_DNA_sequences):
    for DNA_sequence in list_of_DNA_sequences:
        for base in range(len(DNA_sequence)):
            if DNA_sequence[base] ==  'T':
                DNA_sequence[base] =  'A'
            elif DNA_sequence[base] ==  'C':
                DNA_sequence[base] =  'G'
            elif DNA_sequence[base] ==  'G': 
                DNA_sequence[base] =  'C'
            elif DNA_sequence[base] ==  'A':
                DNA_sequence[base] =  'T'
    return list_of_DNA_sequences

def transcribe_DNA(list_of_DNA_sequences):
    list_transcribed_with_DNA_bases = transcribe_with_DNA_bases(list_of_DNA_sequences)
    converted_to_uracil = replace_thymidine_with_uracil(list_transcribed_with_DNA_bases)
    return(converted_to_uracil)

codon_seqs = {      ("A", "A", "A"): "K",   ("A","A","C"):"N",    ("A","A","G"):"K", ("A","A","U"):"N", 
                ("A","C","A"):"T", ("A","C","C"):"T",  ("A","C","G"):"T",   ("A","C","U"):"T", 
                ("A","G","A"):"R", ("A","G","C"):"S",  ("A","G","G"):"R",   ("A","G","U"):"S", 
                ("A","U","A"):"I", ("A","U","C"):"I",  ("A","U","G"):"M",   ("A","U","U"):"I", 

                ("C","A","A"):"Q", ("C","A","C"):"H",  ("C","A","G"):"Q",   ("C","A","U"):"H", 
                ("C","C","A"):"P", ("C","C","C"):"P",  ("C","C","G"):"P",   ("C","C","U"):"P", 
                ("C","G","A"):"R", ("C","G","C"):"R",  ("C","G","G"):"R",   ("C","G","U"):"R", 
                ("C","U","A"):"L", ("C","U","C"):"L",  ("C","U","G"):"L",   ("C","U","U"):"L", 

                ("G","A","A"):"E", ("G","A","C"):"D",  ("G","A","G"):"E",   ("G","A","U"):"D", 
                ("G","C","A"):"A", ("G","C","C"):"A",  ("G","C","G"):"A",   ("G","C","U"):"A", 
                ("G","G","A"):"G", ("G","G","C"):"G",  ("G","G","G"):"G",   ("G","G","U"):"G", 
                ("G","U","A"):"V", ("G","U","C"):"V",  ("G","U","G"):"V",   ("G","U","U"):"V", 

                ("U","A","A"):"_", ("U","A","C"):"Y",  ("U","A","G"):"_",   ("U","A","U"):"T", 
                ("U","C","A"):"S", ("U","C","C"):"S",  ("U","C","G"):"S",   ("U","C","U"):"S", 
                ("U","G","A"):"_", ("U","G","C"):"C",  ("U","G","G"):"W",   ("U","G","U"):"C", 
                ("U","U","A"):"L", ("U","U","C"):"F",  ("U","U","G"):"L",   ("U","U","U"):"F"}



def encode_mRNA(mRNA_seqs):
    list_of_codons= []
    list_of_codon_seqs = []
    for seq in mRNA_seqs:
        seq_len = len(seq)
        num_bp = seq_len / 3
        #print(f'num_bp = {num_bp}')
        seq_array = np.array(seq)
        list_of_arrays = np.split(seq_array, num_bp)
        for codon in list_of_arrays:
            
            list_of_codons = list(map(tuple, list_of_arrays))
            #print(codon)
            
            #print(f'list of codons ={list_of_codons}')

        # #print(f'list_of_codons:{list_of_codons}')
        list_of_codon_seqs.append(list_of_codons)
        #print(list_of_codon_seqs)
    return(list_of_codon_seqs)
#
def turn_codons_into_peptide_sequence(list_of_lists_of_codons):
    #amino_acid_seq = []
    list_peptide_seqs = []
    
    for sequence in list_of_lists_of_codons:
        amino_acid_seq = []
        #print(f'sequence = {sequence}')
        for codon in sequence:
            
            #print(codon)
            codon = codon_seqs[codon]
            #print(codon)
            amino_acid_seq.append(codon)
            #print(amino_acid_seq)
        list_peptide_seqs.append(amino_acid_seq)
        #print(list_peptide_seqs)
        #print(f'list_of_list_of_codons={list_of_lists_of_codons}')

    return list_peptide_seqs


def will_the_peptide_translate(list_of_peptides):
    list_of_translated_peptides = []

    
    for peptide in list_of_peptides:
        #print(peptide)
        
        if '_' in peptide and peptide[0] == 'M':
            stop_loc = peptide.index('_')
            translated_peptide = peptide[0:stop_loc]
            if len(translated_peptide) > 6:
                list_of_translated_peptides.append(translated_peptide)
            #print((f'translated_peptide ={translated_peptide}'))
            
    return list_of_translated_peptides


def list_to_string(list_of_lists_for_stringing_up):
    list_of_seq_strings = []
    for sequence in list_of_lists_for_stringing_up:
        sequence = ''.join(sequence)
        #print(sequence)
        list_of_seq_strings.append(sequence)

    return list_of_seq_strings


def search_protein(list_of_generated_sequences, path_to_database):
    length_seq = []
    sequences = pd.read_excel(path_to_database)
    #display(sequences)  
    aligner = Align.PairwiseAligner(match_score = 1.0) 
    list_of_hits = []
    counter = 0
    for seq1 in list_of_generated_sequences:
        for seq in sequences["Sequence"]:
            #ensures a complete match of the whole sequence
            if len(seq) == len(seq1):
                align_score = aligner.score(seq1, seq)
                if align_score != 0:
                    score = len(seq) / align_score
                if score == 1:
                    print('score!')
                    protein = sequences.iloc[counter, 1]
                    index_protein_hitnumber =(counter, protein, len(list_of_hits))
                    list_of_hits.append(index_protein_hitnumber)
                    
            counter = counter + 1
    print(list_of_hits)
    return list_of_hits
# def BLAST():
    
#     inst = NCBIWWW.qblast('blastp', 
#     'nr', 
#     'QIKDLLVSSSTDLDTTLVLVNAIYFKGMWKTAFNAEDTREMPFHVTKQESKPVQMMCMNNSFNVATLPAEKMKILELPFASGDLSMLVLLPDEVSDLERIEKTINFEKLTEWTNPNTMEKRRVKVYLPQMKIEEKYNLTSVLMALGMTDLFIPSANLTGISSAESLKISQAVHGAFMELSEDGIEMAGSTGVIEDIKHSPESEQFRADHPFLFLIKHNPTNTIVYFGRYWSP', megablast=False,
#     expect=1000,
#     word_size=7,
#     nucl_reward=1,
#     nucl_penalty=-3,
#     )
#     blast_records = NCBIXML.parse(inst)
#     return blast_records
def protein_length_stuff():
    list1 = []
    sequences = pd.read_excel(r"C:\Users\marko\Downloads\uniprot-download_true_fields_accession_2Cprotein_name_2Csequence_for-2022.11.18-14.09.30.99.xlsx")
    for seq in sequences["Sequence"]:
        len_seq = len(seq)
        #print(len_seq)
        list1.append(len_seq)
    return list1
list1 = protein_length_stuff()
from scipy.stats import shapiro 
#SD = np.std(list1)
mean = np.mean(list1)
minus_3sd, plus_3sd = [np.mean(list1) - 3 * np.std(list1), np.mean(list1) + 3 * np.std(list1)]
range = max(list1) - min(list1)
print(minus_3sd, plus_3sd, range, mean)
#it's not normally distributed, it looks to roughly follow a gamma distribution
is_it_normally_distr = shapiro(list1)
print(is_it_normally_distr)
#ecluding the higher outliers (above 3 sigma) which seem to cause issues with graphing
final_list = [x for x in list1 if (x < mean + 1497.1556427819269)]

#create gamma distribution using rsolution between points (1), mean + 3 sigma, number of examples to create 
def create_gamma_dist(resolution, range, num_samples):
    gamma_dist = np.random.gamma(resolution, range, num_samples)
    return gamma_dist

def save_protein_match_score(file_name, data):
    with open(f'{file_name}.json', 'w') as outfile:
    outfile.write(json_string)

gamma_dist = create_gamma_dist(3, mean + 1497.1556427819269, 205788)


y = gamma_dist
plt.hist(gamma_dist, bins =100)
plt.xlabel('length')
plt.ylabel('Frequency')

plt.show() 



# DNA_sequences = create_mulitple_DNA_sequences(10, 2, 32767)
# sequences_transcribed = transcribe_DNA(DNA_sequences)
# into_codons = encode_mRNA(sequences_transcribed)
# #print(translated_into_codons)
# peptide_seqs = turn_codons_into_peptide_sequence(into_codons)
# #print(peptide_seqs)
# translated_peptides = will_the_peptide_translate(peptide_seqs)
# strings = list_to_string(translated_peptides)
# search = search_protein(strings, (r"C:\Users\marko\Downloads\uniprot-download_true_fields_accession_2Cprotein_name_2Csequence_for-2022.11.18-14.09.30.99.xlsx"))
# print(search)






