#Working directory
p="D:/Dropbox/PTMExchange/Documents/Benchmarking paper/GitHub/Example_Search_Files/PXD008355"
#Working directory of each search to be compared
TPP_wd = p + "\TPP" #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
wd = p + "\TPP_pGly"
wd2 = p + "\TPP_pLeu"
wd3 = p + "\TPP_pAsp"
wd4 = p + "\TPP_pGlu"
wd5 = p + "\TPP_pPro"
#location of search database
database = "Z:\hlkramsb\Search_databases\Araport11_crap_targetdecoy.fasta"
#PXD ID
PXD = "PXD000853"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories
w = [TPP_wd, wd, wd2, wd3, wd4, wd5]
#list all search names
s = ["pA","pG","pL","pD","pE","pP"]