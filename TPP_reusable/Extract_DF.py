import pandas as pd
import re
import numpy as np

#Extract PTM prophet output to df
def extract_PTMprophet_df(input,PXD):
	file = open(input, 'r')

	peptides = []
	proteins = []
	scores = []
	PTMs = []
	PTM_scores=[]
	all_positions=[]
	PTM_info = []
	spectrum = []
	usi_all=[]
	sources = []
	all_peptide_mods=[]
	counter = 1
	for line in file:
		if counter == 1:
			counter += 1
		else:
			PTM_score_temp = ""
			data_row = line.strip().split(',')
			peptide = (data_row[1]).replace('"', '')
			peptide=peptide.replace(',',':')
			if peptide !="-":
				peptides.append(peptide)
			if peptide == "-":
				peptides.append(data_row[5])
			scores.append(data_row[0])
			spectrum.append(data_row[2].replace(".mzML",""))
			protein = (data_row[6]).replace('"', '')
			protein = protein.replace(",", ":")

			proteins.append(protein)
			mod_peptide=data_row[11]
			PTM_temp=data_row[10]

			spectrum_temp = (data_row[2]).split('.', 1)
			sources.append(spectrum_temp[0])

			all_scores=[]
			all_scores_temp=[]
			PTM_pos_list=[]
			PTM_score_list_temp=[]
			positions_temp=0
			if peptide != "-":
				pep_temp=peptide.split(":")
				mods_temp=PTM_temp.split(";")
				for p,m in zip(pep_temp,mods_temp):
					pos_temp = 0
					temp = re.split('[()]', p)
					while len(temp) > 1:
						PTM_score = (temp[1])
						positions_temp+=float(PTM_score)
						all_scores_temp.append(PTM_score)
						temp2 = temp[0]
						pos_temp = pos_temp + (len(temp2))
						pos = "x" + str(pos_temp)
						PTM=m
						PTM_score_list_temp.append(PTM_score)
						PTM_pos_list.append(pos_temp)
						PTM_score_temp_2 = (";" + pos + ":" + PTM + ":" + PTM_score[:-1])
						PTM_score_temp = PTM_score_temp + PTM_score_temp_2
						del temp[0:2]
				for i in sorted(all_scores_temp, reverse=True)[:round(positions_temp)]:
					all_scores.append(i)
			else:
				mod_peptide="-"
				PTM_temp="-"

			if mod_peptide[0]=="n":
				mod_peptide=mod_peptide[1:]

			PTM_score_temp = str(PTM_score_temp)
			PTM_score_temp = PTM_score_temp.replace(",'", "")
			PTM_info.append(PTM_score_temp[1:])

			temp = data_row[2].split('.')
			source_temp = temp[0]
			scan_temp = temp[2]
			z = data_row[12]

			position_list=""
			PTM_list=""
			score_list=""
			PTM_position_temp=0
			counter_pos=1
			if "[" not in mod_peptide:
				PTM_list+="-"
			if '[Acetyl]-' in mod_peptide:
				PTM_list+=";Acetyl"
				score_list += ";0"
				position_list+=";0"
			for i in mod_peptide.replace('[Acetyl]-','').split("["):
				if counter_pos==1:
					PTM_position_temp+=len(i)
					counter_pos+=1
				else:
					PTM_list+=";"+(i.split("]")[0])
					position_list+=";"+str(PTM_position_temp)

					if i.split("]")[0]=="Phosphorylation" or i.split("]")[0]=="Oxidation":
						for j,k in zip(PTM_pos_list,PTM_score_list_temp):
							if PTM_position_temp == j:
								score_list+=";"+k
					else:
						score_list+=";0"

					PTM_position_temp+=len(i.split("]")[1])

			USI = "mzspec:" + PXD + ":" + source_temp + ":scan:" + scan_temp + ":" + mod_peptide + "/" + z

			usi_all.append(USI)
			all_peptide_mods.append(mod_peptide)
			PTMs.append(PTM_list[1:])
			all_positions.append(position_list[1:])
			PTM_scores.append(score_list[1:])

	df = pd.DataFrame(
		{"Peptide_mod": all_peptide_mods, "Peptide" : peptides, "Protein": proteins, "Source": sources, "Score": scores, "PTM": PTMs, "PTM_positions":all_positions, "PTM Score": PTM_scores, "PTM_info":PTM_info,"Spectrum": spectrum, "USI": usi_all})
	return df

#Write PTM_info in format position:mod:score
def extract_PTM_info(score_temp,mod,PTM_info):
	score_list=[]
	scores_top=[]
	mod_count=0
	for i in score_temp.split("("):
		for b in i.split(")"):
			if "." in b or "1" in b:
				score_list.append(b)
				mod_count+=float(b)
	scores=""
	for i in sorted(score_list,reverse=True)[:round(mod_count)]:
		scores+=i
		scores_top.append(i)
	for i in scores_top:
		pos=0
		for a in score_temp.split(")"):
			temp=a.split("(")
			pos+=len(temp[0])
			if len(temp)==2:
				if temp[1]==i:
					if round(mod_count)>0:
						PTM_info+=";"+temp[0][-1:]+str(pos)+":"+mod+":"+str(i)
						mod_count+=-1.0
	return(PTM_info)
