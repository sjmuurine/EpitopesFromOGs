# README

This workflow finds epitopes from conserved orthogroups

Scripts for jobs to be run are in a folder “999scripts”. 
Below are explanations and how to run them 

Before that, create a working environment for you. 
A folder for scripts and another for data. 
Upload the file ConservedInGramposGramnegBoth2.fa to the data folder 

# The workflow:

First, fasta sequences are arranged into files according to the orthogroups:

(cd directory/where/fastas/are/ )
python3 /scratch/directory/directory/999scripts/split_orthogroups_v2.py

Then sequences within the orthogroups are aligned with MAFFT and parameter --auto
This allows automatic algorithm selection based on dataset size
in our case it was / will be L-INS-i (Probably most accurate, very slow) 
the script records log files from the alignmet to "mafft_errors.log"
and aligned fastas / orthogroup in a folder "alignments" 

run the alignment like this:
(cd directory/where/arranged/fastas/are/ )
sbatch /scratch/directory/directory/999scripts/batch_mafft_orthogroups2.sh

Now find conserved regions from the alignments.
Folder 999scripts has a script called batch_conservation_v2.py.
The line below runs another bash script, Conserved60minL8merg10.sh, 
that uses batch_conservation_v2.py with parameters min conservation% 60,
min peptide length 8 and merging regions together that are separated 
less than 10 aa. 
You can change the parameters by editing Conserved60minL8merg10.sh  

Run it like this
bash /scratch/directory/directory/999scripts/Conserved60minL8merg10.sh

Now the results are in a folder "conservation_results_60pct_8aa_merge10"
Each orthogroup has its own file and there is also a summary file 
Summary file can be used for identifying orthogroups that can be 
analyzed visually (the aligned fastas, that is)

