import pandas as pd
from _functools import reduce
from os.path import exists
data = pd.read_csv("../data/clusterDD.csv")
data.columns = ["cdr3", "v_gene", "count", "isotype", "cdr3_nt", "sample", "cluster", "file"]

data = data.drop_duplicates()
data["unique_cluster_id"] = data["sample"].str.replace("_sc","") +  data["cluster"].astype("int").astype("str").str.zfill(5)
clusters = data["unique_cluster_id"].unique()  # restore CDR3, V, COUNT, ISOTYPE, NT
# clusters =

"""
cols in data:
CDR3 V counts isotype nt  sample  cluster  file sample_cluster
0    1    2      3     4     5       6      7       8
"""
scSeqs_from_GM_list = []
with open("../data/Tree_with_GM.txt","w") as f:
    for cluster in clusters:
        #merge the same seqs from scXCR and XCR-Seq together respectively

        seqs = data[data.iloc[:, 8] == cluster] #[["v_gene", "isotype", "cdr3_nt", "sample", "cluster"]]
        isotypes = " ".join(seqs["isotype"].unique())

        if "IGHG" in isotypes and "IGHM" in isotypes: #
            # isotype 里面只要有G和M的，把这个cluster中的单细胞数据留下
            scSeqs_from_GM_list.append(seqs[seqs["sample"].str.endswith("sc")])

            print("cluster {} has {}".format(cluster,isotypes))
            f.write("cluster {} has {}\n".format(cluster,isotypes))
#%%
#同时带有IGHM和IGHG的cluster数据
scSeqs_from_GM = reduce(lambda x, y :x.append(y), scSeqs_from_GM_list)
scSeqs_from_GM["sample"] = scSeqs_from_GM["sample"].str.replace("__sc","") #只找单细胞的因为只有单细胞的才有注释，同时规范样本名


# sc_table_dir = Path("./data/scBCR")
# sc_table_dir.rglob("./*BCR_filtered_contig_annotations.csv")

clusters = scSeqs_from_GM["unique_cluster_id"].unique()

cluster_anno = pd.read_table("../data/scBCR/cluster_sample_rename.txt",sep="\t")

mark_file_with_plasma = open("../data/Tree_with_GM_plasma.txt","w")
for cluster in clusters:

    samples = scSeqs_from_GM[scSeqs_from_GM["unique_cluster_id"]==cluster]["sample"].unique()

    for sample in samples:

        sub_table = scSeqs_from_GM[(scSeqs_from_GM["unique_cluster_id"] == cluster)& (scSeqs_from_GM["sample"] == sample)]

        if sample == "S4":
            sub_table["sample"] = "R4"
            sample = "R4"

        sub_table = sub_table.drop_duplicates()
        sub_table["Seq_ID"] = range(len(sub_table))

        sc_table_file = "../data/scBCR/{}/{}_BCR_filtered_contig_annotations.csv".format(sample,sample)

        if exists(sc_table_file):
            sc_table = pd.read_csv(sc_table_file)
            sc_table = sc_table[(sc_table["high_confidence"]==True) & (sc_table["productive"]==True) & (sc_table["is_cell"]==True)]
            bc_sc_table = pd.merge(sub_table, sc_table, left_on=["v_gene", "cdr3_nt"], right_on=["v_gene", "cdr3_nt"], how="left")


            #预处理
            bc_sc_table["barcode"] = bc_sc_table["sample"] +"-" + bc_sc_table["barcode"].str.split("-",expand=True)[0]

            annnoted_sc = pd.merge(bc_sc_table, cluster_anno, left_on="barcode", right_on="cell", how="left")

            if str(27) in annnoted_sc["new_clusterID"].unique().astype("str"):
                mark_file_with_plasma.write(str(cluster)+"\n")
            annnoted_sc.to_csv("{}_{}_annotated".format(cluster, sample))

        else:
            print(RuntimeWarning("{} not found".format(sc_table_file)))

    print("{} Done".format(cluster))

mark_file_with_plasma.close()


#