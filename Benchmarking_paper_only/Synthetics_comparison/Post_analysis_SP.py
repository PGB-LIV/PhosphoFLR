import pandas as pd
import matplotlib.pyplot as plt
import os

#Replace columns with Full match and Partial match (True/False)
def full_comparison(file):
    df = pd.read_csv(file, dtype={'Synthetic_peptide_match': str,'Synthetic_peptide_site_match':str,'Synthetic_peptide_pool-site_match':str})
    all_dict=[]
    partials=[]
    for i in range(len(df)):
        if df.loc[i,'Synthetic_peptide_match']=="TRUE" and df.loc[i,'Synthetic_peptide_pool-site_match']=="TRUE":
            all_dict.append("TRUE")
        else:
            all_dict.append("FALSE")
        if df.loc[i,'Synthetic_peptide_match']=="FALSE" and df.loc[i,'Synthetic_peptide_pool-site_match']=="TRUE":
            partials.append("TRUE")
        elif df.loc[i,'Synthetic_peptide_match']=="FALSE" and df.loc[i,'Synthetic_peptide_pool-site_match']=="FALSE":
            partials.append("FALSE")
        else:
            partials.append("N/A")

    df['Synthetic_full_match']=all_dict
    df['Synthetic_partial_match']=partials
    df2=df.drop(['Synthetic_peptide_match', 'Synthetic_peptide_site_match','Synthetic_peptide_pool-site_match'],axis=1)
    output=file.replace(".csv","_final.csv")
    df2.to_csv(output,index=False)
    return df2

#Add column for if spectrum matches between comparisons
def spectrum_comparison(input,spectrum_match_list):
    df=pd.read_csv(input).fillna('N/A')
    match_counts=[]
    for i in range(len(df)):
        match_count=0
        USI_list=df.loc[i,'All_USI'].split(";")
        for s in USI_list:
            source=s.split(":")[2].strip(" ")
            scan=s.split(":")[4]
            charge=s.split("/")[1].strip(" ")
            spectrum=source+"."+scan+"."+scan+"."+charge
            spectrum=spectrum.replace(".0",".")
            if spectrum in spectrum_match_list:
                match_count+=1
        match_counts.append(match_count)
    df['Spectrum_matches']=match_counts
    output=input.replace(".csv","_spectrum_match.csv")
    df.to_csv(output,index=False)

#Syntehtic FLR: false localisation count/Count
#model FLR from combined probability: 1-(PSM prob*PTM prob)/Count
def model_FLR(input):
    df=pd.read_csv(input,dtype={'Synthetic_full_match': str,'PTM_best_score':float}).fillna('N/A')
    df = df[df.Use_for_SP_analysis == "Y"]
    df = df.reset_index(drop=True)
    for i in range(len(df)):
        if df.loc[i,'Synthetic_full_match'] == "False":
            df.loc[i,'Synthetic_full_match_False']=1
        else:
            df.loc[i, 'Synthetic_full_match_False']=0
    df['PTM_final_prob']=df['Score']*df['PTM_best_score']
    df=df.sort_values(by=(['PTM_final_prob']),ascending=[False])
    df = df.reset_index(drop=True)
    df['Count']=(df.index) + 1
    df['False_count']=df['Synthetic_full_match_False'].cumsum()
    df['True_count']=df['Count']-df['False_count']
    df['FLR']=df['False_count']/df['Count']
    df['q_value'] = df['FLR']
    df['q_value'] = df.iloc[::-1]['FLR'].cummin()
    df['final_temp']=1-df['PTM_final_prob']
    df['PTM_final_prob_FLR']=df['final_temp'].cumsum()/df['Count']
    df=df.drop(['final_temp'],axis=1)
    df['PTM_final_prob_q_value'] = df['PTM_final_prob_FLR']
    df['PTM_final_prob_q_value'] = df.iloc[::-1]['PTM_final_prob_FLR'].cummin()
    df = df.reset_index(drop=True)
    df.to_csv(input.replace(".csv", "_FLR.csv"))

#Sort by PTM score only - Syntehtic FLR: false localisation count/Count
def PTM_sort_FLR(input):
    df=pd.read_csv(input,dtype={'Synthetic_full_match': str,'PTM_best_score':float}).fillna('N/A')
    df = df[df.Use_for_SP_analysis == "Y"]
    df = df.reset_index(drop=True)
    for i in range(len(df)):
        if df.loc[i,'Synthetic_full_match'] == "False":
            df.loc[i,'Synthetic_full_match_False']=1
        else:
            df.loc[i, 'Synthetic_full_match_False']=0
    df=df.sort_values(by=(['PTM_best_score']),ascending=[False])
    df = df.reset_index(drop=True)
    df['Count']=(df.index) + 1
    df['False_count']=df['Synthetic_full_match_False'].cumsum()
    df['True_count']=df['Count']-df['False_count']
    #df['FLR']= df['False_count']/df['True_count']
    df['FLR']=df['False_count']/df['Count']
    lim=df['FLR'].min()-0.01
    df['q_value'] = df['FLR']
    df['q_value'] = df.iloc[::-1]['FLR'].cummin()
    df.to_csv(input.replace(".csv","_FLR_PTM_sort.csv"))

#Plot FLR comparisons - Syntehtic FLR and pX_FLR for each software	
def plot_FLR_comparisons(filter,file_list,software_list,output_folder):
    import seaborn as sns;sns.set_palette("colorblind")
    folder=""
    legend=[]
    for i,j in zip(software_list,file_list):
        folder+="_"+i
        legend.append(i + " Answer Key FLR")
        legend.append(i + " Model FLR")
        if "pL" in j or "pLeu" in j:
            legend.append(i + " pLeu Decoy FLR")
        elif "pGlu" not in j and ("pG" in j or "pGly" in j):
            legend.append(i + " pGly Decoy FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(i + " pAsp Decoy FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(i + " pGlu Decoy FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(i + " pPro Decoy FLR")
        else:
            legend.append(i + " pAla Decoy FLR")
    output=output_folder+"\Comparisons/"+folder[1:]
    if not os.path.exists(output):
        os.mkdir(output)

    x_list=['PTM_final_prob','Count']
    for x in x_list:
        lines = ["blue", "green", "red", "cyan", "magenta", "orange"]
        counter = 0
        max = 0
        for a, b in zip(file_list,lines):
            if "pGly" in a:
                p="pGly_q_value"
            elif "pLeu" in a:
                p="pLeu_q_value"
            elif "pAsp" in a:
                p="pAsp_q_value"
            elif "pGlu" in a:
                p="pGlu_q_value"
            elif "pPro" in a:
                p="pPro_q_value"
            else:
                p="pAla_q_value"
            df=pd.read_csv(a)
            if df[x].max()>max:
                max=df[x].max()
            if x == "PTM_final_prob":
                lim = (1, 0)
            else:
                lim=(0,max)
            c = "C" + str(counter + 1)
            if counter==0:
                ax=df.plot.line(x=x, y="q_value", linestyle="-.",xlim=lim, color=c, figsize=(14,7))
				#add in model FLR if needed
                df.plot.line(x=x, y="PTM_final_prob_q_value", linestyle=":", color=c, ax=ax, xlim=lim)
                df.plot.line(x=x, y=p, linestyle="-", color=c, xlim=lim, ax=ax)
                #ax=df.plot.line(x=x, y=p, linestyle="-", linewidth=3, color=c, xlim=lim, figsize=(14,7))
            else:
                df.plot.line(x=x, y="q_value", ax=ax, linestyle="-.", color=c, xlim=lim)
                df.plot.line(x=x, y="PTM_final_prob_q_value", linestyle=":", color=c, ax=ax, xlim=lim)
                df.plot.line(x=x, y=p, linestyle="-", color=c, ax=ax, xlim=lim)
                #linewidth=3
            counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        ax.set_xlabel("Count of Sites")
        if x=="PTM_final_prob":
            ax.set_xlabel("Combined probability")
        if "Site_confident_PTM_unique" in file_list[0]:
            plt.savefig(output+"/"+filter+"_"+x+"_All_software_Site_Confident_"+p+".jpg",dpi=300)
        else:
            plt.savefig(output + "/" + filter + "_" + x + "_All_software_" + p + ".jpg", dpi=300)


#Plot FLR comparisons for PTM score ordering - Syntehtic FLR and pX_FLR for each software	
def plot_FLR_comparisons_PTM_prob(filter,file_list,software_list,output_folder):
    import seaborn as sns; sns.set_palette("colorblind")
    folder=""
    legend=[]
    for i,j in zip(software_list,file_list):
        folder+="_"+i
        legend.append(i + " Answer Key FLR")
        if "pL" in j or "pLeu" in j:
            legend.append(i + " pLeu Decoy FLR")
        elif "pG" in j or "pGly" in j:
            legend.append(i + " pGly Decoy FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(i + " pAsp Decoy FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(i + " pGlu Decoy FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(i + " pPro Decoy FLR")
        else:
            legend.append(i + " pAla Decoy FLR")
    output=output_folder+"/Comparisons/"+folder[1:]
    if not os.path.exists(output):
        os.mkdir(output)
    x_list = ['PTM_best_score','Count']
    for x in x_list:
        #lines = ["-", ":", '-.', '--']
        lines = ["blue", "green", "red", "cyan", "magenta", "orange"]
        counter = 0
        max = 0
        for a, b in zip(file_list,lines):
            if "pGly" in a:
                p="pGly_q_value"
            elif "pLeu" in a:
                p="pLeu_q_value"
            elif "pAsp" in a:
                p="pAsp_q_value"
            elif "pGlu" in a:
                p="pGlu_q_value"
            elif "pPro" in a:
                p="pPro_q_value"
            else:
                p="pAla_q_value"
            df=pd.read_csv(a)
            if df[x].max()>max:
                max=df[x].max()
            if x == "PTM_best_score":
                lim = (1, 0)
            else:
                lim=(0,max)
            c = "C" + str(counter + 1)
            if counter==0:
                ax=df.plot.line(x=x, y="q_value", linestyle="-.", color=c, xlim=lim, figsize=(14,7))
                df.plot.line(x=x, y=p, linestyle="-", color=c, ax=ax, xlim=lim)
            else:
                df.plot.line(x=x, y="q_value", ax=ax, linestyle="-.", color=c, xlim=lim)
                df.plot.line(x=x, y=p, linestyle="-", color=c, ax=ax, xlim=lim)
            counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        ax.set_xlabel("Count of Sites")
        if x=="PTM_best_score":
            ax.set_xlabel("PTM probability")
        plt.savefig(output+"/"+filter+"_"+x+"_All_software_"+p+"_PTM_Q.jpg",dpi=300)
