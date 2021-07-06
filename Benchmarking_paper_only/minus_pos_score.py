import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

#usage: python minus_pos_score.py minus_output_name positive_output_name arg1 arg2... where arg are comparison files
#eg :python minus_pos_score.py RAPA_minus.jpg RAPA_pos.jpg C:/Users/krams/Dropbox/PTMExchange/CSV_files/RAPA_STY.csv C:/Users/krams/Dropbox/PTMExchange/CSV_files/RAPA_pAla.csv C:/Users/krams/Dropbox/PTMExchange/CSV_files/RAPA_pLeu.csv C:/Users/krams/Dropbox/PTMExchange/CSV_files/RAPA_pGly.csv
file_list=list(sys.argv)[3:]
output_files=list(sys.argv)[1:3]

print(file_list)
print(output_files)


positions=['Minus_pos','Positive_pos']
for p, output in zip(positions,output_files):
    fig, axs=plt.subplots(nrows=2,ncols=2,sharey=False,sharex=False,figsize=(28,14),gridspec_kw={'hspace': 0.25})
    counter=0
    if "Minus" in p:
        label="-1 position"
    else:
        label="+1 position"
    for file in file_list:
        df=pd.read_csv(file) 
        #df=df.head(20)
        if "STY" in file:
            df['PTM_final_prob']=df['PTM_score']*df['Score']
            df=df.loc[df['PTM_final_prob']>=0.75]
            df = df.reset_index(drop=True)
        else:
            # df=df.loc[df['PTM_final_prob']>=0.75]
            if "pAla" in file:
                df=df.loc[df['pAlaq_value']<=0.05]
            if "pLeu" in file:
                df=df.loc[df['pLeu_q_value']<=0.05]
            if "pGly" in file:
                df=df.loc[df['pGlyq_value']<=0.05]
            #print(df['PTM_final_prob'].min())
            df = df.reset_index(drop=True)
             #print(df.columns)
        for i in range(len(df)):
            if df.loc[i,'PTM_positions']==1:
                df.loc[i,'Minus_pos']=" N-term"
            elif df.loc[i,'PTM_positions']>=len(df.loc[i,'Peptide']):
                df.loc[i,'Positive_pos']=" C-term"
            else:
                df.loc[i,'Minus_pos']=df.loc[i,'Peptide'][df.loc[i,'PTM_positions']-2]
                df.loc[i,'Positive_pos']=df.loc[i,'Peptide'][df.loc[i,'PTM_positions']]
            #print(df.loc[i,'Minus_pos'],df.loc[i,'Peptide'],df.loc[i,'PTM_positions'])
            #print(str(i)+"/"+str(len(df)))
        #print(df)    
        
        df2= df.groupby([p], as_index=False)[['PTM_final_prob']].mean()
        if counter==0:
            colors = tuple(np.where(df2[p]=="A", "tab:blue", "tab:blue"))
            df2.plot.bar(x=p,y="PTM_final_prob", legend=False,ax=axs[0,0], color=colors)
            axs[0,0].set_title("a)",loc="left")
            axs[0,0].set_title("STY",loc="center")
            axs[0,0].set_ylim(0.8,1.0)
            axs[0,0].set_ylabel("Average combied probability")
            axs[0,0].set_xlabel(label)
        elif counter==1:
            colors = tuple(np.where(df2[p]=="A", "tab:red", "tab:blue"))
            df2.plot.bar(x=p,y="PTM_final_prob", legend=False,ax=axs[0,1], color=colors)
            axs[0,1].set_title("b)",loc="left")
            axs[0,1].set_title("STYA",loc="center")
            axs[0,1].set_ylim(0.8,1.0)
            axs[0,1].set_ylabel("Average combied probability")
            axs[0,1].set_xlabel(label)
        elif counter==2:
            colors = tuple(np.where(df2[p]=="L", "tab:red", "tab:blue"))
            df2.plot.bar(x=p,y="PTM_final_prob", legend=False,ax=axs[1,0], color=colors)
            axs[1,0].set_title("c)",loc="left")
            axs[1,0].set_title("STYL",loc="center")
            axs[1,0].set_ylim(0.8,1.0)
            axs[1,0].set_ylabel("Average combied probability")
            axs[1,0].set_xlabel(label)
        elif counter==3:
            colors = tuple(np.where(df2[p]=="G", "tab:red", "tab:blue"))
            df2.plot.bar(x=p,y="PTM_final_prob", legend=False,ax=axs[1,1], color=colors)
            axs[1,1].set_title("d)",loc="left")
            axs[1,1].set_title("STYG",loc="center")
            axs[1,1].set_ylim(0.8,1.0)
            axs[1,1].set_ylabel("Average combied probability")
            axs[1,1].set_xlabel(label)        
        counter+=1
    
    plt.savefig(output, dpi=300)
    

    