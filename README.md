# pyboneucleic-motifs

## Requirements

This project requires Python 3; required packages will be covered in 'Getting Started'

## Getting started

To get started, you can simply clone this git repository. To install the required dependencies, run:

```
pip3 install -r requirements.txt
```

This step is not necessary if you don't wish to run the algorithms on data in FASTA format

To run the included tests, run any of the test files with Python:

```
python3 tests.py
python3 fastaTest.py
```

```tests.py``` generates random sequences and inserts a motif randomly into each of them. Each of the implemented algorithms is then run on the generated sequences, and the resulting PWM is printed to the console.

```fastaTest.py``` reads a set of sequences in FASTA format and builds a data array out of them. Each of the implemented algorithms is then run on the input sequences, and the resulting PWM is printed to the console, for varying motif sizes (from 4 to 13). Sequences of 1000bp in promoter regions from the human genome is included (obtained from <http://genome.ucsc.edu/cgi-bin/hgTables> and masked using [RepeatMasker](https://repeatmasker.org))

## Main Library Functions

There are two main library functions, and they have similar signatures:

```
motifEMOOPS(sequences: List[bytearray], k: int, bgFreqs: Dict[int, float]) -> List[Dict[int, float]]

motifGibbsOOPS(sequences: List[bytearray], k: int, bgFreqs: Dict[int, float]) -> List[Dict[int, float]]
```

```sequences``` is a list of sequences, which are each a bytearray. Strings can be converted using the ```bytearray(string)``` function. 

```k``` is the size of motif to search for

```bgFreqs``` is a dictionary that contains the background frequencies for each character in the input sequences. Since the input sequences are byte arrays, the input dictionary maps from ascii values for each character (integers) to their probability.

Both return a list of length ```k```, where each element of the list is a dictionary that represents the probability of each character to appear at that position. The format is similar to the format of ```bgFreqs```

## Printing Generated Motifs

A helper function is provided to print a PWM representation as returned from either of the motif searching functions:

```
printMotif(pwm, alpha, k)
```

Where ```pwm``` is the result of one of the above functions, alpha is either a dictionary in the same format as ```bgFreqs``` or alternately an iterable list of bytes or integers that represent the keys for the dictionary at each position of ```pwm```, and ```k``` is the length of the motif.