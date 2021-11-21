import sys
import pandas as pd

def is_unlocalized(x):
    if str(x)=="nan":
        return False
    for v in x.split("|"):
        if v[0]=="x":
            return True
    return False

def get_mass(x):
    if x["is_unlocalized"] == False:
        return x["peptide_mass"]
    return x["precursor_mass"]

def remove_extension(x):
    return ".".join(x.split(".")[:-1])

d = []
for fn in sys.argv[1:]:
    data = pd.read_csv(fn)
    data = data[(data["q-value"] <= 0.01) & (data["database"] == "T")]
    data["is_unlocalized"] = data["modifications"].apply(is_unlocalized)
    data["Full Sequence"] = data["matched_peptide"] + data["modifications"].astype(str)

    proteins = pd.read_csv(fn.replace(".csv",".proteins.csv"))
    proteins = proteins[proteins["protein_group_q-value"]<=0.01]
    proteins = proteins[proteins["is_shared_peptide"]==False]
    proteins = proteins[["ionbot_match_id","protein_group"]].drop_duplicates()
    data = data.merge(proteins,on="ionbot_match_id",how="left")

    tmp = pd.DataFrame()
    tmp["Scan Retention Time"] = data["observed_retention_time"] / 60
    tmp["Precursor Charge"] = data["charge"]
    tmp["Base Sequence"] = data["matched_peptide"]
    tmp["Full Sequence"] = data["Full Sequence"]
    tmp["Peptide Monoisotopic Mass"] = data.apply(get_mass,axis=1)
    tmp["Protein Accession"] = data["protein_group"]
    tmp["File Name"] = data["spectrum_file"].apply(remove_extension)

    d.append(tmp)

# For PSMs with an unlocalized modification ('x') I use the median
# of all precursor masses for each such peptidoform.
# FlashLFQ does not accept peptidoforms with different masses.
data = pd.concat(d)
data["Peptide Monoisotopic Mass"] = data.groupby("Full Sequence")["Peptide Monoisotopic Mass"].transform('median')

data.to_csv("flashlfq.tsv",sep="\t",index=False)
                