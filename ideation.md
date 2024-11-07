## Codon optimization
The problem: trying to minimize a bunch of different variables, to give the global min 


Using the sliding windows option

Windows are created from previous 

Within each window check:
1. Forbidden sequences
    - Should accept input from user about what restriction enzymes they are using
    - RNAse cleavage sites, nucleases etc
    - internal start, termination sequences
    - internal ribosome binding sites

transcription termination
global GC content
length of longest rare-codon stretch (needs: what codons are rare for this organism)
- Not using rare codons (for the organism being expressed in) as much
- Having a CAI beforehand??
homopolymeric runs
20+ bp regions of high/low GC content
Secondary structure: repetitive sequences that can recombine

Doing this with dependency injection: trying to start with general methods, and plugging in information about the organism we are working in.

Starting with a 



protein letter to id
id to protein letter

Doing these operations in matrices?


## implementation
selecting window
monte carlo search over the codons in the window, minimizing score
selecting codons
moving window


How many before the codon should I consider?
How many after the codon should I consider?


selecting codon at random from a probability distribution (CAI)
-> generalize later into CC score??

Need to start small -> get larger
Setting up dependency injection now, modularity


RBSChooser: finds best RBS for the sequence by checking first 6 AA of the CDS
MCTS: Creates optimal CDS for a given peptide sequence

Returns: Transcript object that has an RBSOption, the peptide sequence (as a string), and a list of the codons


already have checkers, just need to use them. 


Checking for forbidden sequences between the chosen RBS and the designed CDS

everything gets run on proteome_benchmarker.




# Actual
Don't have to do all the bs stuff above, just move the window, randomly generate the codons from the list, and run the checkers on the sequence. 

For my custom checker:
- doing global/local gc content checker?
- asdf

have to copy/paste my RBSChooser into this one.


guided random codon choice algorithm
- input: taking a sequence of amino acids, and a sequence of dna
    - we get 9 bp previous, 9 bp that the algorithm chooses, and the next 6 amino acids from the window
    - but i need to make this more general, so no assumptions about length
- output: 9bp for the desired codons

given this sequence of aas, find a sequence of bps that work

loop in here: 
- random generation
- checking for absolute nos

master checker function

sliding window loop
- dealing with first codons (theres no 9bp previous)

## After solution
I can go to real monte carlo tree search algorithms
I wonder how complex it would be for an actual program, how long it would take to run. Could I make optimizations?


## MODIFYING CODON PERCENTAGES SO THEY ADD TO 1
I got lazy... I can go back and rewrite this function if needed to make a proper distribution. Add the last column to the lists like now, then find sum of the list and divide each entry by that sum. Then i'm sure I can find some rounding function that makes them all 1, if the float values weren't working out. `

Increased
AGG 0.01% from 0.03 to 0.04
ATC 0.01% from 0.40 to 0.41

Decreased:
AGC 0.01% from 0.27 to 0.26


## 20241030
How should I refactor everything?

- Code is independent of organism
    - Don't need heuristics from E. coli to run, using just the supplied files to make conclusions about everything.
- Figure out sampling
    - How do I increase diversity above 0.5 and keep rare codons below 3, while also maximizing CAI
    - "These are the bounds, find the highest CAI"
    






Found values using the top 187 most expressing genes:
avg diversity: 0.24393129387900409
avg CAI: 0.36561386565057086
avg rare codons: 0.49732620320855614
Num True (good sequences): 5
Num False (bad sequences): 182


Using new checker (len(codon_counts) / 62 )
For the most expressing genes:
avg diversity: 0.7157150250129384
avg CAI: 0.36561386565057086
avg rare codons: 0.49732620320855614
Num True: 167, Num False: 20


