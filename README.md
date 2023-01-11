# Biology_fun

This project aims to simulate biological molecular biology processes using Python code. An overview of the project workflow is in the below image:
![image](https://user-images.githubusercontent.com/107410852/211831557-46748cc1-a487-456d-8693-f54b107d2adf.png)

## DNA sequence generation

The first task was to generate a function which produces DNA sequences. The `def randomly_generate_DNA_sequence` function takes the maximum DNA sequence length as a positional argument. There is also an argument to declare whether to include methionine DNA as the first 3 bases, and a further argument for stop codon as the last 3 bases. This function can finally take into consideration the natural biases of bases. 

~~~def randomly_generate_DNA_sequence(DNA_seq_max_length, met_first=False, stop_last = False, bp_weighting=False):

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
    return DNA_sequence~~~
    
 
