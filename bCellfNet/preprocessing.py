import pandas as pd
from pyvis.network import Network
from copy import deepcopy
from math import *
from scipy.stats import fisher_exact, ttest_ind
from functools import reduce
from multipy.fdr import qvalue

def draw_net(package):
    sample, data = package
    net = Network()
    nodes = targeted_isotpye
    colors = list(map(lambda x: isotype_color[x], nodes))
    # nodes =
    # sizes = data.groupby("isotype")["cluster"].c  ount()[targeted_isotpye][nodes].to_list()
    # net.add_nodes(nodes, color = colors, size = sizes)
    net.add_nodes(nodes, color=colors)

    tmp = deepcopy(nodes)

    edges = []

    sample_num = len(data)
    max = 0
    while len(tmp) > 0:
        element_from = tmp.pop(0)
        for element_to in tmp:
            if set([element_from, element_to]) == set(["IGHA1", "IGHA2"]):
                continue
            set_x = set(data[data["isotype"] == element_from]["cluster"])
            set_y = set(data[data["isotype"] == element_to]["cluster"])

            x_and_y =  len(set_x.intersection(set_y))
            x_and_Ny = len(set_x - set_y)
            Nx_and_y = len(set_y - set_x)

            Ny_and_Nx = sample_num - len(set_x.union(set_y))
            odd, p = fisher_exact([[x_and_y,Nx_and_y],[x_and_Ny,Ny_and_Nx]])

            # p = -log10(p)
            # p = 2 if p < -log10(0.01) else p  #卡阈值为0.01
            width = x_and_y
            # width = 1 if width > 0 and width < 1 else width

            edges.append([element_from, len(set_x), element_to, len(set_y), width, p])

    edges = pd.DataFrame(edges,columns = ["from", "from_size","to", "to_size","weight", "p"])
    edges["qValue"] = qvalue(edges["p"])[1]
    edges["Sig"] = ["*" if y <0.05 else "" for y in edges["qValue"]]
    edges["sample"] = sample
    edges = edges.set_index("sample", append=True)
    # edges.to_csv("./data/{}.csv".format(sample))
    # edges = edges.set_index(["index",sample])
    return edges
    # edges = list(map(lambda x: x[:-1]+[1+10*x[-1]/max], edges))
    # net.add_edges(edges)
    # net.show("test.html")

allCluster = pd.read_csv("../data/allClusters.csv")

allCluster = allCluster.drop(0).drop_duplicates()

# allCluster["V7"] = pd.Categorical(allCluster["V7"], order = allCluster.sort_values()..astype("int").unique())

allCluster_simplified = allCluster[["V4","V6","V7"]]

allCluster_simplified.columns = ["isotype", "sample", "cluster"]

allCluster_simplified.loc[:,"sample"] = allCluster_simplified["sample"].str.split("_",expand=True)[0]

allCluster_simplified = allCluster_simplified[allCluster_simplified["isotype"].astype("str") != "nan"]

allCluster_simplified["sample"].unique()

grouped = allCluster_simplified.groupby("sample")

group_test = grouped.get_group("S1")

isotype_color = dict(zip(
    ['IGHG1', 'IGHA1', 'IGHGP', 'IGHM', 'IGHG3', 'IGHG2','IGHA2', 'IGHD', 'IGHG4', 'None'],
    ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"]))

targeted_isotpye = [ "IGHM", "IGHG3", "IGHG1", "IGHG2", "IGHA1", "IGHA2"] # 按发生的先后顺序来

summary = list(map(draw_net, grouped))

summary = reduce(lambda x,y: x.append(y), summary)


summary.to_csv("./outs/cs_network.csv")

isotypes_tmp = deepcopy(targeted_isotpye)

without_bacterial = ["S4", "S7", "S12", "S14", "S16"]

sig_result = []

class_switch_order = [ "IGHM", "IGHG3", "IGHG1", "IGHG2", "IGHA1", "IGHA2"]

while len(isotypes_tmp) > 0:
    element_from = isotypes_tmp.pop(0)
    for element_to in isotypes_tmp:
        if set([element_from,element_to]) == set(["IGHA1", "IGHA2"]):
            continue
        # if element_from > element_to:
        #     tmp = element_from
        #     element_from = element_to
        #     element_to = tmp
        values = summary[(summary["from"]==element_from) &(summary["to"]==element_to)].reset_index()
        values_with_bac = values[~values["sample"].isin(without_bacterial)]["p"]
        values_without_bac = values[values["sample"].isin(without_bacterial)]["p"]
        var, p = ttest_ind(values_with_bac,values_without_bac)
        sig_result.append([element_from,  element_to, p])

sig_result = pd.DataFrame(sig_result,columns=["from","to","p"])
sig_result["qValue"] = qvalue(sig_result["p"])[1]
sig_result = sig_result.sort_values(["qValue","p"])
sig_result.to_csv("./outs/classSwitch_sig_bacterial.csv",index=False)