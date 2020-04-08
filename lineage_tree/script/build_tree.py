import pandas as pd
data = pd.read_csv("../data/clusterDD.csv")
data.columns = ["cdr3", "v", "count", "isotype", "nt", "sample", "cluster", "file"]




#%%
# data = data.dropna()
data = data.drop_duplicates()
data["unique_cluster_id"] = data["sample"].str.replace("_sc","") +  data["cluster"].astype("int").astype("str").str.zfill(5)
clusters = data["unique_cluster_id"].unique()  # restore CDR3, V, COUNT, ISOTYPE, NT
# clusters =

"""
cols in data:
CDR3 V counts isotype nt  sample  cluster  file  sample_cluster
0    1    2      3     4     5       6      7          8
"""
#%%
for cluster in clusters:

    #merge the same seqs from scXCR and XCR-Seq together respectively
    seqs = data[data.iloc[:, 8] == cluster][["cdr3", "v", "isotype", "nt","count"]]
    seqs = seqs.groupby(["cdr3", "v", "isotype", "nt"]).sum().reset_index() #单细胞和混合的合一块
    fasta = open("../data/fasta/cluster" + cluster + ".fasta", "w")

    # lens = seqs.iloc[:, 4].map(lambda x: len(x))
    #
    # # check the number of seqs in each cluster
    # if len(lens.unique()) > 1 :
    #     print("【warning】 cluster {}".format(cluster))
    # else :
    #     print("The lens of elements in cluster {} is {}".format(cluster, lens.to_list()))

    for index, row in seqs.iterrows():

        name_col_selected = ["cdr3","v", "count", "isotype"]
        seq_name = list(row[name_col_selected].apply(lambda x: str(x)))
        seq_name.append("Seq_" + str(index))
        name = ">" + "_".join(seq_name)

        fasta.writelines(name + "\n")
        fasta.writelines(row["nt"] + "\n")
    fasta.close()

print("Done")


#%%
from pathlib import Path
dir = Path("../data/fasta")

import time, subprocess
fastas = list(dir.rglob("cluster*.fasta"))

errLog = open("../data/clustalw/log", "w")

for fasta in fastas:

    cluster = fasta.name.lstrip(fasta.suffix).zfill(5)

    # print("d【{}】{} start...".format(time.ctime(), cluster))

    out_file = "../data/clustalw/" + cluster + ".msa.nt.fa"
    cmd = ["../source/clustal-omega",
           "--force",
           "-i", fasta,
           "-o", out_file]
    # print("【{time}】{cmd}".format(time=time.ctime, cmd=cmd))

    subprocess.run(cmd, stderr=errLog)
    # with open("../data/clustalw/{}.msa.stdout".format(cluster), "w") as stdout:
    #     subprocess.run(cmd, stderr = errLog, stdout = stdout)

errLog.close()


#%%
from pathlib import Path
import os
import time, subprocess


try:
    os.makedirs("./data/igphyml")
except :
    pass

os.chdir("./data/igphyml")

dir = Path("../clustalw")
fastas = list(dir.rglob("cluster*.fa"))

errLog = open("log", "w")

count = 0

total_count = len(fasta)

for fasta in fastas:


    cmd = ["../source/igphyml/src/igphyml",
           "-i", fasta,
           "-m", "GY",
           "--threads","4"]
    # print("【{time}】{cmd}".format(time=time.ctime, cmd=cmd))

    with open("{}.igphyml.stdout".format(fasta), "w")as stdout: #这个文件会在这里的
        subprocess.run(cmd, stderr=errLog, stdout=stdout)

    print("【{}】{}/{} cluster {} is done".format(time.ctime(), count, total_count, fasta))

errLog.close()

