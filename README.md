# LEAF - (phytopLasmas Effector prediction via rAndom Forest/tobe changed)
--description--




## Usage
LEAF can be used as a stand-alone script or with the available singularity3.7 container (recommended)

### LEAF from container
To properly use LEAF you can download the directory and execute the LEAF1.0.sh file in the container

**---LEAF.sh help---**
```
-first argument: -dft/pre_computed_feature_table.tsv	use "-dft" (do feature table) if no pre-computed feature table is available (LEAF will start the prediction of the features and feature table build-up)
							use /path/to/pre_computed_feature_table.tsv otherwise
							(the column in feature_table must be in the same order of those in training_feature_tables.tsv) 
if -dft:
	-second argument: path/to/protein_sequence.fasta	the input file in FASTA format containing AA sequences (can be a selection of proteins or an entire proteome)
	-third argument: path/to/output_directory			ouput directory in which to save both feature predictions/feature table/LEAF putative effector prediction 
	-fourth argument: suffix/prefix						to distinguish the current run of LEAF (e.g. strain name, CaPmali_AT)

otherwise:
	-second argument: path/to/output_directory			ouput directory in which to save LEAF putative effector prediction 
	-third argument: prefix								to distinguish the current run of LEAF (e.g. strain name, CaPmali_AT)
```
**Container Usage**
```
singularity exec -B binding/dirs LEAF1.0.simg /opt/LEAF.sh -dft /path/to/aa_sequences.fasta /path/to/output_dir suffix/prefix
```
**or**
```
singularity exec -B /binding/dirs LEAF1.0.simg /opt/LEAF.sh /path/to/feature_table.tsv /path/to/output_dir prefix
```
### LEAF stand-alone 
The required python3.8.10 libraries are:
- biopython
- pandas
- joblib

The required software to be installed separately or used by other containers (e.g. singularity), are:
- SignalP - v4.1
- TMHMM - v2.0
- MobiDB-lite

To properly use LEAF you can download the directory and execute the following steps:

**Predict features** 
```
signalp4.1 -f long -s notm -t gram+ -T /path/to/tmpdir /path/to/aa_sequences.fasta > signalp_out.txt
```
```
tmhmm /path/to/aa_sequences.fasta > tmhmm_out.txt
```
```
python mobidb-lite.py -bin mobidb-lite/binx -l /path/to/aa_sequences.fasta -o mobidb_out.txt
```
**or**
```
singularity exec -B /binding/dirs /path/to/signalp4.1.simg signalp4.1 -f long -s notm -t gram+ -T /path/to/tmpdir /path/to/aa_sequences.fasta > signalp_out.txt
```
```
singularity exec -B /binding/dirs /path/to/tmhmm2.0.simg tmhmm /path/to/aa_sequences.fasta > tmhmm_out.txt
```
```
singularity exec -B /binding/dirs /path/to/mobidb-lite.simg python /opt/mobidb-lite.py -bin /opt/mobidb-lite/binx -l /path/to/aa_sequences.fasta -o mobidb_out.txt
```

**Build-up the feature table**
```
python3.8.10 ./build_feature_table.py -i /path/to/aa_sequences.fasta\
					-o /path/to/output_dir/feature_table_name\
					-sp /path/to/feature_prediction_dir/signalp_out.txt\
					-tm /path/to/feature_prediction_dir/tmhmm_out.txt\
					-mb /path/to/feature_prediction_dir/mobidb_out.txt\
					-pr ./pre_feature_prediction/prosite_eff.fasta\
					-fte ./training_feature_tables/feature_table_eff_std.tsv\
					-prm ./pre_feature_prediction/prosite_motifs_profiles_eff.txt\
					-ms ./pre_feature_prediction/monster_score_eff.tsv\
					-mmc ./pre_feature_prediction/df_motif_CLUMPs_eff.tsv
```
**Predict putative effector proteins with LEAF**
```
python3.8.10 ./LEAF1.0.py -ft /path/to/feature_table_name.tsv -o /path/to/output_dir -px distinguishable_name 
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




