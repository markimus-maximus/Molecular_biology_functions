# Biology_fun

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

Inspection of the lengths of peptide sequences in the human peptide database revealed that there was a distribution resembling a gamma distribution (see figure below).
![image](https://user-images.githubusercontent.com/107410852/211844300-b85bf911-1cc5-4d5d-856d-65754692c092.png)

Accordingly, a function `create_DNA_sequence_with_gamma_lengths` was written to generate a list of sequences with lengths representative of the natural lengths of peptides. This function uses much of the same code as the `randomly_generate_DNA_sequence` function, wth the addition of taking arguments for defining the shape and scale and using these as arguments for `np.random.gamma(shape,scale)`.

The `create_multiple_DNA_sequences` function combines the above functions to create a number of DNA sequences according to the `number_DNA_sequences` argument.

The below figure with gamma distribution looks comparable in distribution to the peptide lengths above.

![image](https://user-images.githubusercontent.com/107410852/211847130-0b7bf42c-8d90-4fc6-9e3b-9829a29ce72e.png)

## Transcribing DNA sequences

Two functions `replace_thymidine_with_uracil(list_of_DNA_sequences)` and `transcribe_DNA_bases(list_of_DNA_sequences)` allow for transcription of the DNA into RNA.



    
    
 
