import sys
import os
import argparse
import pandas as pd

def is_unlocalized(x):
    if str(x)=="nan":
        return False
    for v in x.split("|"):
        if v[0]=="x":
            return True
    return False

def get_mass(x):
    return x["peptide_mass"]

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

def contains_unknown_modification(x):
    if str(x) == "nan":
        return False
    return "unknown modification" in x

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser(description='ionbot.quant: conver ionbot result files to FlashLFQ input file')
parser.add_argument('-x', action="store_true", dest="do_not_filter",
                    help="don't filter peptide matches with unexpected modification")
parser.add_argument('folder', type=dir_path, nargs='+',help='folder with ionbot results')
args = parser.parse_args()

d = []
for fn in args.folder:
    #ionbot results are filtered for identified matches (q-value<=0.01) that don't contain an unknown modification
    data = pd.read_csv(fn+"/ionbot.first.csv")
    data = data[(data["q-value"] <= 0.01) & (data["database"] == "T")]
    data["contains_unknown_modification"] = data["modifications"].apply(contains_unknown_modification)
    data = data[data["contains_unknown_modification"]==False]
    print(len(data))
    if args.do_not_filter == False:
        data = data[data["unexpected_modification"].isnull()]
    print(len(data))

    data["modifications"].fillna("",inplace=True)
    data["Full Sequence"] = data["matched_peptide"] + data["modifications"].astype(str)

    # the proteins inference result file is used to:
    # - remove peptides shared between protein groups
    # - remove proteins not identified at 1% FDR
    proteins = pd.read_csv(fn+"/ionbot.first.proteins.csv")
    proteins = proteins[proteins["protein_group_q-value"]<=0.01]
    proteins = proteins[proteins["is_shared_peptide"]==False]
    not_shared_psms = list(proteins["ionbot_match_id"].unique())
    not_shared_psms = dict.fromkeys(not_shared_psms,True)
    selected_proteins = list(proteins["protein"].unique())
    selected_proteins = dict.fromkeys(selected_proteins,True)

    data["not_shared"] = data["ionbot_match_id"].apply(is_not_shared)
    data = data[data["not_shared"]==True]

    tmp = pd.DataFrame()
    tmp["Scan Retention Time"] = data["observed_retention_time"] / 60
    tmp["Precursor Charge"] = data["charge"]
    tmp["Base Sequence"] = data["matched_peptide"]
    tmp["Full Sequence"] = data["Full Sequence"]
    tmp["Peptide Monoisotopic Mass"] = data.apply(get_mass,axis=1)
    tmp["Protein Accession"] = data["proteins"].apply(get_proteins)
    tmp["File Name"] = data["spectrum_file"].apply(remove_extension)

    d.append(tmp)

data = pd.concat(d)

data.to_csv("flashlfq.tsv",sep="\t",index=False)
                
