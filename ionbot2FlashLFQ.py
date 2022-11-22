from dataclasses import replace
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
    return x["precursor_mass"]

def get_proteins(x):
    tmp = x.split('||')
    buf = []
    for p in tmp:
        pp = p.split('((')[0]
        if pp in selected_proteins:
            buf.append(pp)
    if len(buf) == 0:
        return "empty"
    return ";".join(buf)       

def is_not_shared(x):
    return x in not_shared_psms

def remove_extension(x):
    return ".".join(x.split(".")[:-1])

def add_delta_mass(x):
    if "unknown modification" in x["modifications"]:
        delta_mass = x["peptide_mass"] - x["precursor_mass"]
        return x["modifications"].replace("unknown modification","unknown_modification_%s"%str(delta_mass))
    else:
        return x["modifications"]

d = []
for fn in sys.argv[1:]:
    data = pd.read_csv(fn)
    data = data[(data["q-value"] <= 0.01) & (data["database"] == "T")]
    data["modifications"].fillna("",inplace=True)
    data["modifications"] = data.apply(add_delta_mass,axis=1)
    data["Full Sequence"] = data["matched_peptide"] + data["modifications"].astype(str)
    print(len(data))


    # the proteins inference result file is used to:
    # - remove peptides shared between protein groups
    # - remove proteins not identified at 1% FDR
    proteins = pd.read_csv(fn.replace(".csv",".proteins.csv"))
    proteins = proteins[proteins["protein_group_q-value"]<=0.01]
    proteins = proteins[proteins["is_shared_peptide"]==False]
    not_shared_psms = list(proteins["ionbot_match_id"].unique())
    not_shared_psms = dict.fromkeys(not_shared_psms,True)
    selected_proteins = list(proteins["protein"].unique())
    selected_proteins = dict.fromkeys(selected_proteins,True)

    data["not_shared"] = data["ionbot_match_id"].apply(is_not_shared)
    data = data[data["not_shared"]==True]
    print(len(data))

    tmp = pd.DataFrame()
    tmp["Scan Retention Time"] = data["observed_retention_time"] / 60
    tmp["Precursor Charge"] = data["charge"]
    tmp["Base Sequence"] = data["matched_peptide"]
    tmp["Full Sequence"] = data["Full Sequence"]
    tmp["Peptide Monoisotopic Mass"] = data.apply(get_mass,axis=1)
    tmp["Protein Accession"] = data["proteins"].apply(get_proteins)
    tmp["File Name"] = data["spectrum_file"].apply(remove_extension)

    d.append(tmp)
    break

# For PSMs with an unlocalized modification ('x') I use the median
# of all precursor masses for each such peptidoform.
# FlashLFQ does not accept peptidoforms with different masses.
data = pd.concat(d)
data["Peptide Monoisotopic Mass"] = data.groupby("Full Sequence")["Peptide Monoisotopic Mass"].transform('median')

data.to_csv("flashlfq.tsv",sep="\t",index=False)
                
