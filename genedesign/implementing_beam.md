# BeamSearch
need to make all the checkers work on a single pass through the string




## Functions
### initiate

### Run method
Inputs:
peptide (amino acid string), ignores (set of RBSOptions to ignore)

returns:
tuple with the best RBSOption and a list of codons

### score

### violates_constraints
Checks for forbidden sequences, internal promoters


## Facts
### code stuff
- RBSOptions are generated using the first 6 aas of the generated sequence
- RBSChooser already takes into account secondary structure in the UTR when calculating strength

### Expression level
- The folding energy of the entire mRNA was not significantly correlated with fluorescence
    - Folding energy of the first third of the mRNA was strongly correlated
    - Region from -4 to +37 relative to start predicted folding energy explained 44% of the variation based on T7 promoter set
- forbidden sequences anywhere are an instant no-go
- RNAse E cleavage sites predicts 4.7% of variation in expression
    - means I can include the num of predicted sites in the scoring function
- Secondary structure in the rbs (dealt with by RBSChooser secondary structure scoring)
- CAI, codon usage, rare codon counts in the CDS
    - done by CodonChecker
