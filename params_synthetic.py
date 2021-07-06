#Working directory of search files
p="D:/Dropbox/PTMExchange/Documents/Benchmarking paper/GitHub/Example_Search_Files/PXD007058"
#Working directory of each search to be compared
TPP_wd = p + "/TPP" #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
Peaks_wd = p + "/PEAKS"
MQ_wd = p + "/MaxQuant"
PD_wd = p + "/PD"
#Location of .mgf search files, if searching PEAKS. Can be left blank. 
search_files_location = "Z:\hlkramsb\PXD007058\Peaks\Input files"
#location of search database
database = "Z:/hlkramsb/Search_databases/Araport11_crap_targetdecoy.fasta"
#Synthetic peptide fasta with pools
s_pep_pools = "D:/Dropbox/PTMExchange/Documents/Benchmarking paper/GitHub/Benchmarking_paper_only/Synthetics_comparison/pool_ALL_name.fasta"
#Synthetic peptide fasta without pools
s_pep = "D:/Dropbox/PTMExchange/Documents/Benchmarking paper/GitHub/Benchmarking_paper_only/Synthetics_comparison/pool_ALL_name_no_pool.fasta"
#PXD ID
PXD = "PXD007058"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories
w = [TPP_wd, Peaks_wd, MQ_wd, PD_wd]
#list all software names
s = ["TPP", "PEAKS", "MaxQuant", "PD"]