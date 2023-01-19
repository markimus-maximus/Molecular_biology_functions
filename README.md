# Molecular_biology_functions

This project aims to simulate biological molecular biology processes using Python code. An overview of the project workflow is in the below image:
![image](https://user-images.githubusercontent.com/107410852/211831557-46748cc1-a487-456d-8693-f54b107d2adf.png)

## DNA sequence generation

The first task was to generate a function which produces DNA sequences. The `def randomly_generate_DNA_sequence` function takes the maximum DNA sequence length as a positional argument. There is also an argument to declare whether to include methionine DNA as the first 3 bases, and a further argument for stop codon as the last 3 bases. This function can finally take into consideration the natural biases of bases. 

~~~
def randomly_generate_DNA_sequence(DNA_seq_max_length, met_first=False, stop_last = False, bp_weighting=False):

    #randomly generate a number between range
    
        
    DNA_length = np.random.randint(153, high=DNA_seq_max_length)
    #print(f'DNA length: {DNA_length}')
    
    DNA_length_3 = get_triplicate_length(DNA_length)
    #print(f'new DNA length: {DNA_length_3}')
    DNA_sequence = []
    if bp_weighting==False: 
        for base_pair in range(1, DNA_length_3):
            base_pair = np.random.choice(DNA_bases)
            #print(base_pair)
            DNA_sequence.append(base_pair)
    else:
        for base_pair in range(1, DNA_length_3):
            base_pair = np.random.choice(DNA_bases, p=(0.293, 0.207, 0.20, 0.30))
            #print(base_pair)
            DNA_sequence.append(base_pair)
    if met_first == True:
        DNA_sequence[0:2] = ['T', 'A', 'C']
        #print(DNA_sequence)
    if stop_last == True:
        DNA_sequence[-3:] = ['A', 'C', 'T']
    #print(DNA_sequence)
    #print(f'len DNA sequ ={len(DNA_sequence)}')
    return DNA_sequence
~~~   
An accessory function was written to convert a randomly generated DNA length to a triplicate length, in order to reflect the correct lengths of DNA.

~~~
def get_triplicate_length(given_length):
    nearest_divisible_by_3 = round(given_length/3)
    #print(f'nearest_div = {nearest_divisible_by_3}')
    nearest_3_multiple = nearest_divisible_by_3 * 3 
    #print(f'nearest_multiple = {nearest_3_multiple}')
    return nearest_3_multiple
~~~

To get an idea of the length distribution of protein sequences in the human peptide database, `get_protien_length_list` function was written which returns a list of protein lengths. Inspection of the lengths of peptide sequences revealed that there was a distribution resembling a gamma distribution (see figure below).
![image](https://user-images.githubusercontent.com/107410852/211844300-b85bf911-1cc5-4d5d-856d-65754692c092.png)

Accordingly, a function `create_DNA_sequence_with_gamma_lengths` was written to generate a list of sequences with lengths representative of the natural lengths of peptides. This function uses much of the same code as the `randomly_generate_DNA_sequence` function, wth the addition of taking arguments for defining the shape and scale and using these as arguments for `np.random.gamma(shape,scale)`.

The `create_multiple_DNA_sequences` function combines the above functions to create a number of DNA sequences according to the `number_DNA_sequences` argument.

The below figure with gamma distribution looks comparable in distribution to the peptide lengths above.

![image](https://user-images.githubusercontent.com/107410852/211847130-0b7bf42c-8d90-4fc6-9e3b-9829a29ce72e.png)

## Transcribing DNA sequences

Two functions `replace_thymidine_with_uracil(list_of_DNA_sequences)` and `transcribe_DNA_bases(list_of_DNA_sequences)` allow for transcription of the DNA into RNA. The returned in a list of RNA sequences which are the complimentary strand of the DNA sequences.

## Translating RNA into amino acid sequence

The function `encode_mRNA(mRNA_seqs)` takes an RNA sequence and splits it into chunks of 3 with the following code, where `num_amino_acids` is the number of splits:
~~~
for seq in mRNA_seqs:
        seq_len = len(seq)
        num_amino_acids = seq_len / 3
        #print(f'num_amino_acids = {num_amino_acids}')
        seq_array = np.array(seq)
        list_of_arrays = np.split(seq_array, num_amino_acids)
~~~
Next, a list of codons is returned with:
~~~
list_of_codons = list(map(tuple, list_of_arrays))
list_of_codon_seqs.append(list_of_codons)
~~~


Code was next required to convert lists of codons into peptide sequences, using the codon sequence as a dictionary key to get the associated amino acid (in single letter form annotation). This single letter is then appended to a list of peptide sequences. To achieve this the `turn_codons_into_peptide_sequence` function was created which takes a list of codon lists as an argument and iteratively works through each codon list to convert the codon into an amino acid, according to the below dictionary:

~~~
                ("A","A","A"):"K", ("A","A","C"):"N",  ("A","A","G"):"K",   ("A","A","U"):"N", 
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
~~~
In addition to using DNA sequences to create peptide sequences, a further approach was to use codons as the starting point for p for peptide sequence generation. The advantage to this approach is that codons can be chosen according to their natural biases. In approaching this, `make_seq_from_codons` function was written. Similarly to creating DNA sequences above, a positional argument is the maximum DNA length. In addition, there are arguments available to define the first codon as a methionine reside (met_first argument), and an option to remove stop codons (to prevent premature truncation of sequences, no_stops argument). Finally, as with DNA sequence generation (above) there is an option to use gamma distribution for sequence lengths (gamma_dist argument). Interestingly, having stop codons available to choose from when constructing the peptide sequences generated a gamma distribution without the need for applying gamma distribution. It is therefore tempting to say that the random insertion of stop codons which occurs naturally results in the gamma distribution observed of peptide sequences.

## Processing peptide sequences and searching for matches
The function `will_the_peptide_translate` was written to determine if the protein would translate by assessing if there was a methionine in residue position 1 and stop codon in the last residue position. Room for improvement here could be to search for an in frame methionine if it is not the first residue position, and generate the protein sequence from here. 

To enable search functions with protein sequences the generated peptide sequences with `list` datatype were converted to `str` sequence with the `list_to_string` function. After creating a list of DNA sequence string, the search_protein function was written. This function takes the list of protein sequences and a path to a protein database to return a list of matches. To achieve the sequence alignments, the `Align` class from the `Bio` library was utilised. Later, the needleman-wunsch algorithm was written which can create aligignments between proteins.



    
 
