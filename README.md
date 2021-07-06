# PhosphoFLR
Scripts for estimating false localisation rates (FLR) using two methods, model FLR and Decoy amino acid FLR. A further FLR estimation can be calculated, Answer Key FLR, while searching the snthetic dataset.

# Usage

Module files contained in *Benchmarking_paper_only* are used to compare the software invesitgated in the manuscript(TPP, Peaks, MaxQuant and Mascot). For the analysis of the Synthetic Peptide Dataset(PXD007058), module files are contained within subfolder *Synthetics_comparison*. 

*TPP_reusable* folder contains supported module code for the analysis of Datasets using TPP, not including the Synthetic Peptide dataset. 

For the analysis of datasets searching with different software (PEAKS, MaxQuant and/or Mascot/ptmRS):
	When searching PEAKS results, before running the FLR comparison script, PEAKS A-Scores muust first be mapped to probabilities using *Benchmarking_paper_only/map_scores_to_prob.R* R script. 

	    $py Software_comparison params_software
    
    The params file (params_software) must be given and must containt the directories of each of the search file locations, including the mgf search file location when searching PEAKS. 
    Pyteomics must be installed - https://pypi.org/project/pyteomics/
	

For the analysis of the Synthetic Peptide Dataset (PXD007058):

      $py Synthetic_comparison params_synthetic
    
    The params file (params_synthetic) must be given and must contain the directories of each of the search file locations, including location of synthetic pepetide fasta files (pool_ALL_name.fasta and pool_ALL_name_no_pool.fasta). 
    If searching with PEAKS, again, the scores must first be mapped to probabilities using the Benchmarking_paper_only/map_scores_to_prob.R R script and the mgf location must be specified in the params file.
    Pyteomics must be installed - https://pypi.org/project/pyteomics/
	
  
For searching all other datasets, with TPP only:

      $py TPP_comparison parmas_TPP	
    The params file (params_TPP) must be given and must contain the directories of each of the search file locations. 


*Benchmarking_paper_only* folder also conains some aditional analysis scripts used in the manuscript:

minus_pos_score.py - gives plots of average score for each amino acid in the +1/-1 positions for each search file

	  $py minus_pos_score.py output_minus.jpg output_positive.jpg comparison_file_STY.csv comparison_file_pAla.csv comparison_file_pLeu.csv comparison_file.pGly.csv
	
PSM_mod_counts.py - gives boxplots of PSM coutns of targets and decoys at p>0.95 and p<0.95 and boxplots of mod counts between targets and decoys at p>0.95 and p<0.95.

	  $py PSM_mod_counts.py comaprison_file.csv

pSTY_X_directional.py - Comparison of STY and decoy amino acid distribution. STY search and output location specified in files_direction.

	  $py pSTY_X_directional.py input_comparison_file.csv output.jpg
    
# Example files

Example files are contained within */Example_Search_files* and are called within the params files for each of the example searches
https://www.dropbox.com/sh/iyvra02fphecx1j/AABA56MNAfSA2kQYkdNHfOOta?dl=0
