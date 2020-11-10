#This script requires two files: 1) Mutation file 2) Protein database 
#Usage:python3 MutantDatabaseCreation.py <Fasta file> <mutation file>
#Note that this script will not work on single letter amino acids codes
#The input file  Mutation file should be in following format
#Mutation file
gene1	Cys33Arg
gene2	Gly26Val
gene3	Trp343Arg
###########
#Protein database
#This script will only work on uniprot database. The database file should be in following format
>header1	seq1
>header2	seq2
>header3	seq3