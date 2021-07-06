import pandas as pd
import os
import time
import glob

import Benchmarking_paper_only.FDR as FDR
import TPP_reusable.Confident_PTMs as Confident_PTMs
import Benchmarking_paper_only.Synthetics_comparison.Unique_PTMs_SP as Unique_PTMs_SP
import Benchmarking_paper_only.extract_scan as extract_scan
import TPP_reusable.convert_xml_sax as convert_xml_sax
import TPP_reusable.Comparison as Comparison
import Benchmarking_paper_only.Synthetics_comparison.Post_analysis_SP as Post_analysis_SP
import Benchmarking_paper_only.Synthetics_comparison.SP_comparison as SP_comparison
import TPP_reusable.pAla as pAla
import sys
from params_synthetic import *

start_time = time.time()
Ascore_cutoff = 0
sub = "FDR_" + str(FDR_cutoff) + "_PTM_score_" + str(Ascore_cutoff)
decoy_list = []
for wd, software in zip(w, s):
    os.chdir(wd)
    if "PEAKS" in software:
        input = "DB search psm.csv"
        results_file = "DB search psm_edit.csv"
        # extract scan from peaks output
        if os.path.isfile(results_file):
            print("Scan already found")
        else:
            extract_scan.scan(search_files_location, input, results_file, wd)
            print("Scan found")
        search_software = "PEAKS"
    elif "TPP" in software:
        # Run pepXML converter - if csv not in dir
        results_file = "interact.ptm.pep.csv"
        if os.path.isfile(results_file):
            print("XML converted file exists")
        else:
            xml = glob.glob("*.ptm.pep.xml")
            try:
                print("Converting XML file")
                convert_xml_sax.convert(xml[0])
                print("XML converted")
            except ValueError:
                raise ValueError("Please provide TPP '.ptm.pep.xml' results file")
        search_software = "TPP"
    elif "MaxQuant" in software:
        results_file = "evidence.txt"
        if not os.path.isfile(results_file):
            print("Please provide location for 'evidence.txt' MaxQuant results file")
        search_software = "MaxQuant"
    elif "PD" in software:
        text_file = glob.glob("*PSMs.txt")
        if len(text_file) != 1:
            raise ValueError('should be only one PSM txt file in the PD directory')
        results_file = text_file[0]
        search_software = "PD"
    # calculate FDR
    if not os.path.exists(sub):
        os.mkdir(sub)
    FDR_output = sub + "/FDR_output.csv"
    if "pL" in software or "pL" in wd or "pLeu" in wd:
        search = "STYL"
        decoy = "pLeu"
    elif "pG" in software or "pG" in wd or "pGly" in wd:
        search = "STYG"
        decoy = "pGly"
    elif "pD" in software or "pD" in wd or "pAsp" in wd:
        search = "STYD"
        decoy = "pAsp"
    elif "pE" in software or "pE" in wd or "pGlu" in wd:
        search = "STYE"
        decoy = "pGlu"
    elif "pP" in software or "pP" in wd or "pPro" in wd:
        search = "STYP"
        decoy = "pPro"
    else:
        search = "STYA"
        decoy = "pAla"
    decoy_list.append(decoy)
    FDR.calculateFDR(results_file, search_software, FDR_output, PXD, wd, search)
    print("FDR done")
    print("--- %s seconds ---" % (time.time() - start_time))
    # filter for FDR (and A-score) cutoff
    confident_PTM_output = sub + "/PTM_confident.csv"
    Confident_PTMs.confident(FDR_output, database, FDR_cutoff, Ascore_cutoff, confident_PTM_output)
    print("Confident PTMs done")
    print("--- %s seconds ---" % (time.time() - start_time))
    # Collapse for best scoring for each peptide, protein site, mass shift on protein or no collapse
    unique_peptide = sub + "/Peptide_confident_PTM_unique.csv"
    unique_site = sub + "/Site_confident_PTM_unique.csv"
    unique_mass = sub + "/Peptide_mass_confident_PTM_unique.csv"
    non_collapse = sub + "/All_confident_PTM_no_collapse.csv"
    # If synthetic dataset, use alternative analysis (include pool etc.)
    if PXD == "PXD007058":
        Unique_PTMs_SP.unique(confident_PTM_output, unique_peptide, unique_site, unique_mass, non_collapse)
    else:
        sys.exit("Use Non Synthetics Analysis Pipeline, either Software comparison or TPP comparison")
    print(FDR_cutoff, software, "--- %s seconds ---" % (time.time() - start_time))

# Compare spectrum between different searches
files = []
for i in w:
    files.append(i + "/" + sub + "/PTM_confident.csv")
Comparison.comparison(p, sub, files, s)
print("Spectrum Comparison Done --- %s seconds ---" % (time.time() - start_time))

if PXD == "PXD007058":
    # Compare synthetics search to known answer key
    SP_comparison.answer_key_comparison(s_pep_pools, s_pep, w, sub)
    file_input_list = ["All_confident_PTM_no_collapse_Synthetic_Peptides_w-pool_Site-based.csv"]
else:
    sys.exit("Use Non Synthetics Analysis Pipeline, either Software comparison or TPP comparison")
files = []
for i in w:
    files.append(i + "/" + sub + "/")
software_list_edit = ""
for i in s:
    software_list_edit += i + "_"
os.chdir(p)
if not os.path.exists("Comparisons"):
    os.mkdir("Comparisons")
# load in list of hits where spectrum comparisons gave match between all searches
spectrum_match_list = []
spectrum = p + "\Spectrum_Comparisons/" + software_list_edit[
                                          :-1] + "/" + sub + "_spectrum_comparison_" + software_list_edit[:-1] + ".csv"
df = pd.read_csv(spectrum, dtype=str)
for i in range(len(df)):
    if df.loc[i, 'Match'] == "TRUE":
        spectrum_match_list.append(df.loc[i, 'Spectrum'])

print("Starting FLR calculations --- %s seconds ---" % (time.time() - start_time))
# Post analysis - FLR calulations and plots
files = []
for i in w:
    for f in file_input_list:
        files.append(i + "/" + sub + "/" + f)

        # print(i + "/" + sub)
for file in files:
    if PXD == "PXD007058":
        Post_analysis_SP.full_comparison(file)
        output = file.replace(".csv", "_final.csv")
        Post_analysis_SP.spectrum_comparison(output, spectrum_match_list)
        plot_input = output.replace(".csv", "_spectrum_match.csv")
        Post_analysis_SP.model_FLR(plot_input)
        Post_analysis_SP.PTM_sort_FLR(plot_input)
    else:
        sys.exit("Not Synthetic Analysis. Please use either Software comparison or TPP comparison")
if PXD == "PXD007058":
    final_file_list = []
    final_PTM_file_list = []
    for i, decoy in zip(w, decoy_list):
        pAla.calulate_decoy_FLR(
            i + "/" + sub + "/All_confident_PTM_no_collapse_Synthetic_Peptides_w-pool_Site-based_final_spectrum_match_FLR.csv",
            decoy)
        pAla.calulate_decoy_FLR(
            i + "/" + sub + "/All_confident_PTM_no_collapse_Synthetic_Peptides_w-pool_Site-based_final_spectrum_match_FLR_PTM_sort.csv",
            decoy)
        final_file_list.append(
            i + "/" + sub + "/All_confident_PTM_no_collapse_Synthetic_Peptides_w-pool_Site-based_final_spectrum_match_FLR_" + decoy + ".csv")
        final_PTM_file_list.append(
            i + "/" + sub + "/All_confident_PTM_no_collapse_Synthetic_Peptides_w-pool_Site-based_final_spectrum_match_FLR_PTM_sort_" + decoy + ".csv")
    Post_analysis_SP.plot_FLR_comparisons(sub, final_file_list, s, p)
    Post_analysis_SP.plot_FLR_comparisons_PTM_prob(sub, final_PTM_file_list, s, p)
else:
    sys.exit("Not Synthetic Analysis. Please use either Software comparison or TPP comparison")
print("FLR calculations done --- %s seconds ---" % (time.time() - start_time))
