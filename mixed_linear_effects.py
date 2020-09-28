### Katharina Sielemann ###
### kfrey@cebitec.uni-bielefeld.de ###
### v1 ###

#all list.txt files contain (I) the dispensability score in the first column, separated by a comma from (II) the second column containing either gene length or exon number or distance to closest TE gene

#imports
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

#------------ FAMILIES
ap2 = pd.read_csv("ath_ap2_scores_exonnumber_list.txt", header = None)
ap2.rename(columns = {0 : "ds", 1 : "exon_number"}, inplace=True)
ap2["family"] = "ap2"

WRKY = pd.read_csv("ath_wrky_scores_exonnumber_list.txt", header = None)
WRKY.rename(columns = {0 : "ds", 1 : "exon_number"}, inplace=True)
WRKY["family"] = "WRKY"

MYB = pd.read_csv("ath_myb_scores_exonnumber_list.txt", header = None)
MYB.rename(columns = {0 : "ds", 1 : "exon_number"}, inplace=True)
MYB["family"] = "MYB"

all_exonnumbers = ap2.append(WRKY)
all_exonnumbers = all_exonnumbers.append(MYB)


ap2 = pd.read_csv("ath_ap2_scores_length_list.txt", header = None)
ap2.rename(columns = {0 : "ds", 1 : "length"}, inplace=True)
ap2["family"] = "ap2"

WRKY = pd.read_csv("ath_wrky_scores_length_list.txt", header = None)
WRKY.rename(columns = {0 : "ds", 1 : "length"}, inplace=True)
WRKY["family"] = "WRKY"

MYB = pd.read_csv("ath_myb_scores_length_list.txt", header = None)
MYB.rename(columns = {0 : "ds", 1 : "length"}, inplace=True)
MYB["family"] = "MYB"

all_lengths = ap2.append(WRKY)
all_lengths = all_lengths.append(MYB)

all_family_df = all_exonnumbers.merge(all_lengths, on = ["ds", "family"])
all_family_df.to_csv("families_data.csv")

md = smf.mixedlm("ds ~ exon_number+length", all_family_df, groups=all_family_df["family"])

mdf = md.fit()
print(mdf.summary())

#---- ap2
ap2_df = all_family_df.query("family == 'ap2'").copy()

md = smf.mixedlm("ds ~ exon_number+length", ap2_df, groups=ap2_df["family"])

mdf = md.fit()
print(mdf.summary())

#---- wrky
wrky_df = all_family_df.query("family == 'WRKY'").copy()

md = smf.mixedlm("ds ~ exon_number+length", wrky_df, groups=wrky_df["family"])

mdf = md.fit()
print(mdf.summary())

#---- myb
myb_df = all_family_df.query("family == 'MYB'").copy()

md = smf.mixedlm("ds ~ exon_number+length", myb_df, groups=myb_df["family"])

mdf = md.fit()
print(mdf.summary())


#---------- ALL GENES
scores_distanceTE = pd.read_csv("scores_distanceTE_list.txt", header=None)
scores_distanceTE.rename(columns = {0:"score", 1:"distance_TE"}, inplace=True)

scores_length = pd.read_csv("scores_genelength_list.txt", header=None)
scores_length.rename(columns = {0:"score", 1:"length"}, inplace=True)

scores_exonnumber = pd.read_csv("scores_exonnumber_list.txt", header=None)
scores_exonnumber.rename(columns = {0:"score", 1:"exon_number"}, inplace=True)


#--------- distance to TE (27247 genes)
scores_distanceTE["group"] = 1

md = smf.mixedlm("score ~ distance_TE", scores_distanceTE, groups=scores_distanceTE["group"])
mdf = md.fit()
print(mdf.summary())


#-------- gene length and exon number (30125 genes)
data_df = scores_length.merge(scores_exonnumber, on = ["score"])
data_df["group"] = 1

md = smf.mixedlm("score ~ length+exon_number", data_df, groups=data_df["group"])
mdf = md.fit()
print(mdf.summary())
