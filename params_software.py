#Working directory of search files
p="D:/Dropbox/PTMExchange/Documents/Benchmarking paper/GitHub/Example_Search_Files/PXD008355"
#Working directory of each search to be compared
TPP_wd = p + "/TPP" #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
Peaks_wd = p + "/PEAKS"
MQ_wd = p + "/MaxQuant"
PD_wd = p + "/PD"
#Location of .mgf search files
search_files_location = "Z:/hlkramsb/Athaliana_data/PXD008355"
#location of search database
database = "Z:/hlkramsb/Search_databases/Araport11_crap_targetdecoy.fasta"
#PXD ID
PXD = "PXD008355"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories
w = [TPP_wd, Peaks_wd, MQ_wd, PD_wd]
#list all software names
s = ["TPP", "PEAKS", "MaxQuant", "PD"]