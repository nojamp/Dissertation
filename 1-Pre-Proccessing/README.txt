1-Pre-Proccessing.R

Input:
- phenotype data
- data file linking genetic IDs and phenotype IDs
- swapped.sample file (genetic data IDs)

Output:
- proccessed phenotype data
- genetic IDs file of eligible individuals
- summary of data

Function:
- unfies units for age
- creates season variable
- codes NAs for negative or zero metabolite values (and all other vars)
- Only keeps individuals with clinical, metabolite, and genetic data from at least one clinical
- summarises variables