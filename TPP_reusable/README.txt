For searching all datasets, excluding PXD007058 Synthetic Dataset, with TPP only:
	$py TPP_comparison parmas_TPP

	The params file (params_TPP) must be given and must contain the directories of each of the search file locations.


TPP_comparison.py calls the following modules:

TPP_reusable/convert_xml_sax
    Convert ".ptm.pep.xml" TPP PTMprophet results output to .csv file for downstream analysis
    output="interact.ptm.pep.csv"
	
TPP_reusable/FDR
    Calculate false detection rate statistics and plots
    -Calls Extract_DF
        Loads CSV file to dataframe for FDR calculations
    Writes folder for downstream analysis eg. "FDR_0.01_PTM_score_0" FDR filter 0.01, no PTM score filter
    output="FDR_output.csv","FDR.jpg","FDR_filter.jpg","FDR_score.jpg"
	
TPP_reusable/Confident_PTMs
    Filters for FDR cutoff specified in "params_TPP"
    -Calls Extract_peptides_DB
        Creates a dictionary of peptides within the search database (specified in "parmas_TPP")
    Gives PTM positions on proteins for PSMs passing FDR threshold
    output="PTM_confident.csv"
	
TPP_reusable/Unique_PTMs
    Gives collapsed outputs of unique hits
        a) collapse by mass shift - collapsing PSMs by peptide+mass shift, keeping the best scoring: "Peptide_mass_confident_PTM_unique.csv"
        b) collapse by peptide  - collapsing PSMs by peptide, keeping best scoring: "Peptide_confident_PTM_unique.csv"
        c) collapse by protein position - collapsing modification by position on protein, keeping best scoring: "Site_confident_PTM_unique.csv"
    Also gives none-collapsed in same format: "All_confident_no_collapse.csv"
	
TPP_reusable/Comparison
    Gives a spectrum by spectrum comparison of all search files specified in "params_TPP"
    Creates a subfolder within "Spectrum_Comparisons" with the search names (specified within "params_TPP") separated by "_" (eg. pA_pG_pL_pD_pE_pP)
    output="FDR_[cutoff]_PTM_score_0_spectrum_comparison_[searchnames].csv" (where each row represents each spectrum and gives the peptide, protein, score, mods et. found in each search for that spectrum),"All_confident_PTM_no_collapse_Site-based_spectrum_match.csv" (binary spectrum match column, do search comparisons provide same spectra results)
		
TPP_reusable/Post_analysis
    Expands PSMs to site-based format, where one row is one site on a peptide. Calculates and generates figures for FLR calculations.
	Output="All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR.csv", can also be used for filtering for no-choice peptides ("All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_filtered.csv") or sorting by PTM score rather than combined probability ("All_confident_PTM_no_collapse_Site-based_spectrum_match_FLR_PTM_sort.csv") or collapsing by peptide+mod ("All_confident_PTM_no_collapse_Site-based_spectrum_match_new_FLR_collapse.csv") or collapsed by protein position ("Site_confident_PTM_unique_FLR.csv") 
	These output files are then processed by "pAla" (see below) to provide the decoy amino acid FLR calculations
	Final FLR calculations are then passed through either "plot_FLR_comparisons" or "plot_FLR_comparisons_PTM_prob" functions 
	Output folder is generated "Comparisons/[software names, seperated by "_"]/" and multiple comparison plots, from each of the ordering/filtering options seen above, are generated here
	
TPP_reusable/pAla
    Calculates decoy amino acid FLR 
	Takes input from "FLR files" generated from Post_analysis and calculates Decoy amino acid FLR. 
	output=[input]_p[Decoy amino acid].csv
	
	