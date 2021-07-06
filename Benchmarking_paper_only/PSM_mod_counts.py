import pandas as pd
import matplotlib.pyplot as plt
import sys

#file="C:/Users/krams/Dropbox/PTMExchange/PXD000612/pA/FDR_0.01_PTM_score_0/All_confident_PTM_no_collapse_Site-based_spectrum_match_new_FLR_pAla.csv"
file=list(sys.argv)[1]

def spectrum_counts(input):
    df = pd.read_csv(input)
    df = df.loc[df['PTM_final_prob'] >= 0.75]
    for i in range(len(df)):
        peptide = df.loc[i, 'Peptide_mod']
        df.loc[i, 'Alt_mod_count'] = peptide.count('[') - peptide.count('[Phosphorylation]')
        df.loc[i, 'Deam_count'] = peptide.count('[Deamidation]')
        df.loc[i, 'ox_count'] = peptide.count('[Oxidation]')
        df.loc[i, 'pyro_count'] = peptide.count('[Pyro')
        df.loc[i, 'acetyl_count'] = peptide.count('[Acetyl]-')

    df1 = df.loc[df['PTM_final_prob'] >= 0.95]
    df2 = df.loc[df['PTM_final_prob'] < 0.95]

    ax = df1.groupby("DecoyP").boxplot(column=['PSM_count'], showfliers=True)
    ax[0].set_xticklabels(['PSM count'])
    ax[1].set_xticklabels(['PSM count'])
    ax[0].set_title('Target p>=0.95')
    ax[1].set_title('Decoy p>=0.95')
    output = input.replace(".csv", "_p0.95_PSM_outliers.jpg")
    plt.tight_layout()
    plt.savefig(output, bbox_inches = "tight", dpi=300)

    ax2 = df2.groupby("DecoyP").boxplot(column=['PSM_count'], showfliers=True)
    ax2[0].set_xticklabels(['PSM count'])
    ax2[1].set_xticklabels(['PSM count'])
    ax2[0].set_title('Target p<0.95')
    ax2[1].set_title('Decoy p<0.95')
    output = input.replace(".csv", "_p0.75_PSM_outliers.jpg")
    plt.tight_layout()
    plt.savefig(output, bbox_inches = "tight", dpi=300)

    ax2 = df2.groupby("DecoyP").boxplot(column=['Alt_mod_count', 'Deam_count', 'ox_count', 'pyro_count', 'acetyl_count'], showfliers=True)
    ax2[0].set_xticklabels(['Other mod count', 'Deamidation count', 'Oxidation count', 'Pyro- count', 'Acetyl count'],rotation=90)
    ax2[1].set_xticklabels(['Other mod count', 'Deamidation count', 'Oxidation count', 'Pyro- count', 'Acetyl count'],rotation=90)
    ax2[0].set_title('Target p<0.95')
    ax2[1].set_title('Decoy p<0.95')
    output = input.replace(".csv", "_p0.75_mods.jpg")
    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight", dpi=300)

    ax = df1.groupby("DecoyP").boxplot(column=['Alt_mod_count', 'Deam_count', 'ox_count', 'pyro_count', 'acetyl_count'], showfliers=True)
    ax[0].set_xticklabels(['Other mod count', 'Deamidation count', 'Oxidation count', 'Pyro- count', 'Acetyl count'],rotation=90)
    ax[1].set_xticklabels(['Other mod count', 'Deamidation count', 'Oxidation count', 'Pyro- count', 'Acetyl count'],rotation=90)
    ax[0].set_title('Target p>=0.95')
    ax[1].set_title('Decoy p>=0.95')
    output = input.replace(".csv", "_p0.95_mods.jpg")
    plt.tight_layout()
    plt.savefig(output, bbox_inches="tight", dpi=300)
    
spectrum_counts(file)