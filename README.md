# LEAF - (phytopLasmas Effector prediction via rAndom Forest/tobe changed)
--description--

## Usage
LEAF can be used as a stand-alone script or with a singularity3.7 container (recommended)
### LEAF stand-alone 
The required python3.8.10 libraries are:
- biopython
- pandas
- joblib

To properly use LEAF you can download the directory and run the LEAF1.0.sh file
```
LEAF1.0 help

.LEAF.sh.swp
-dtf   do the feature table (if you do not have yet a well formatted feature table for the sequences of interest, use this command to make LEAF do it for you)
sequence.fasta
```


## Feature Considered
Protein features considered for classification are:
- Protein length
- Signal Peptide (SignalP 4.1)
- Transmembrane Region (TMHMM)
 - AA in transmembrane domain
 - first 60 AA
 - Probability of N-in
 - Warning signal sequence 
-  Intrinsically disordered regions (MobiDB-lite)
- Motifs/Profiles in AA sequence (Prosite)
- CLUMPs (MOnSTER)
- bin of sequence and presence of CLUMPs

## Classification
- Methods: **Random Forest**, **XGBoost**, **Gaussian Naive Bayes**, **Multinomial Naive Bayes**
- Overall data-set partition and validation: **5-fold cross validation** 
- Goodness of model:
  - Accuracy
  - F-measure
  - Precision-Recall
  - Feature importance (SHAP)
- number of variables (features): 30




