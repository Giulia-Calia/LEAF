#!/bin/bash

# 1 --> command line parameter to ask for the build up of the feature_table.tsv or not '-dft'/'path/to/pre_computed_feature_table.tsv', meaning 'asking yes'/'asking not'

if [ $1 == "-dft" ]
then

# 2 --> 'path/to/protein_sequence.fasta' | the input file in FASTA format containing AA sequences (can be a selection of proteins or an entire proteome)
# 3 --> 'path/to/output_directory' | ouput directory in which to save both output and metadata from LEAF application 
# 4 --> 'suffix'/'prefix' to distinguish the current run of LEAF (e.g. strain name, CaPmali_AT) 

	mkdir $3/tmp
	mkdir $3/feature_prediction
	signalp4.1 -f long -s notm -t gram+ -T $3/tmp $2 > $3/feature_prediction/signalp_$4.txt
	tmhmm $2 > $3/feature_prediction/tmhmm_$4.txt
	python /opt/mobidb-lite.py -bin /opt/mobidb-lite/binx -l $2 -o $3/feature_prediction/mobidb_$4.txt
	python3.8 /opt/build_feature_table.py -i $2 -o $3/feature_table_$4 -sp $3/feature_prediction/signalp_$4.txt -tm $3/feature_prediction/tmhmm_$4.txt -mb $3/feature_prediction/mobidb_$4.txt -pr /opt/pre_feature_prediction/prosite_eff.fasta -fte /opt/pre_feature_prediction/feature_table_eff.tsv -prm /opt/pre_feature_prediction/prosite_motifs.txt -ms /opt/pre_feature_prediction/monster_score.tsv -mmc /opt/pre_feature_prediction/df_motif_CLUMPs.tsv
	python3.8 /opt/LEAF1.0.py -ft $3/feature_table_$4.tsv -o $3 -px $4

else

# 2 --> output directory for prediction.tsv file
# 3 --> prefix used to distinguish the current run of LEAF
	python3.8 /opt/LEAF1.0.py -ft $1 -o $2 -px $3
	
fi
