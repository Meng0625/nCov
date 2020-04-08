#%%
from pathlib import Path
from ete3 import *
import math
dir = Path("../data/clustalw")
forest = list(dir.rglob("*_igphyml_tree.txt"))

isotype_color = dict(zip(
    ['IGHG3', 'IGHA1',  'IGHM','IGHG1', 'IGHG2','IGHA2', 'IGHD', 'IGHG4', 'None'],
    ["#E64B35","#4DBBD5", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#4B0082", "#7E6148", "#303030"]))

# isotype_color = {
#     'IGHM':"#696969",
#     'IGHG3':"#b22222",
#     'IGHG1':"#C71585",
#     'IGHG2':"#FF4500",
#     'IGHG4':"#FFA500",
#     'IGHA1':"#006400",
#     'IGHA2':"#32CD32",
#     'IGHD':"#7B68EE",
#     'IGHGP':"#800080",
#     'None':"#FFF0F5"
# }


def layout(node:TreeNode):

    node.img_style["size"] = 0

    if node.is_leaf():
        node.dist = 0.1
        node.img_style["shape"] = "circle"

        try:
            node_color = isotype_color[node.isotype]
        except:
            print("wanring {} occurs".format(node))
            node_color = "#FFFFFF"

        node.img_style["fgcolor"] = node_color
        node.img_style["shape"] = "circle"
        if hasattr(node, "hz_width"):
            node.img_style["size"] = node.hz_width
        node.img_style["vt_line_width"] = 2
        # node.img_style["hz_line_width"] = node.hz_width
        node.img_style["hz_line_width"] = 2

        if node.aa == annotation_node_showed_as:
            node.aa = node.aa.rjust(max_aa_len)
            node.img_style["hz_line_color"] = "#FFFFFF"
            node.img_style["vt_line_color"] = "#FFFFFF"
            # anno_sequence = "".join(list(map(lambda x: str(math.ceil(x/3)) if (x%3==1 and math.floor(x/3)%2 ==0 and math.ceil(x/3)<10) else "…" if x%6==1 else " " if math.floor(x/3)%2==1  else "\t", range(max_nt_len))))
            anno_sequence = "".join(list(map(lambda x: " " if math.floor(x/3)%2==1  else " ", range(max_nt_len))))
            node.add_feature("sequence",anno_sequence)

        if hasattr(node, "aa"):
            faces.add_face_to_node( faces.TextFace("{} ".format(node.aa),fsize=15 , fgcolor = node_color,bold=True), node, 0,position="aligned" )

        if hasattr(node,"sequence"):
            SequenceFace =  faces.SequenceFace(node.sequence,"nt",18)
            faces.add_face_to_node(SequenceFace, node, 1, aligned=True)

    else:
        node.dist = 0.005 if node == t or node.is_root() else 0.05
        if node.is_root() :
            node.img_style["hz_line_color"] = "#FFFFFF"
            node.img_style["vt_line_color"] = "#FFFFFF"

    if hasattr(node, "in_GMA_path") and node.in_GMA_path:
        node.img_style["hz_line_color"] = "#FF0000"
        node.img_style["vt_line_color"] = "#FF0000"



def find_closest_ans(node:Tree):
    if node.is_root():
        return None
    parent = node.up
    aunts = parent.get_sisters()
    if parent.name != "" and hasattr(parent,"sequence"):
        return parent.sequence
    elif aunts != []:
        for aunt in aunts:
            if aunt.name != "" and hasattr(aunt, "sequence"):
                return aunt.sequence
    return find_closest_ans(parent)


mutation_rank = {}

for tree_file in forest:

    fa_file = Path(tree_file.as_posix().rstrip("_igphyml_tree.txt"))
    cluster = tree_file.name.split(".")[0].lstrip("cluster")

    print("{} :working wisht {}".format(time.ctime(), cluster))

    t = Tree(tree_file.as_posix())

    aligned_seqs = SeqGroup(fa_file.as_posix())
    PhyloTree.link_to_alignment(alignment = aligned_seqs, self = t)

    #
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    ts.show_scale = False
    # ts.title.add_face(TextFace("The color denotes the isotpye.\nThe size denotes frequency.".format(cluster),bold=True), column=1)
    ts.branch_vertical_margin = 0
    ts.legend_position = 3


    T = Tree()
    T.add_child(t)
    T.add_child(name="annotation")
    annotation_node_showed_as = "Translation Frame"
    anno_node =  T&"annotation"
    anno_node.add_feature("freq", int(0))
    anno_node.add_feature("aa", annotation_node_showed_as)
    anno_node.add_feature("isotype", "None")
    anno_node.add_feature("hz_width", 0)
    max_aa_len = 0
    max_nt_len = 0

    #calculating mutations...
    nor_mutations_sum = 0
    for node in t.traverse():
        if node.is_leaf():
            closest_ans_seq = find_closest_ans(node)
            if closest_ans_seq != None and hasattr(node,"sequence"):
                mis_match = 0
                for i in range(len(node.sequence)):
                    if node.sequence[i] != closest_ans_seq[i]:
                        mis_match = mis_match + 1
                nor_mutations_sum  =  nor_mutations_sum + mis_match/len(node.sequence)
            else:
                if not node.up.is_root:
                    print("ERROR: node {} has no parents".format(node.name))
                #     print("{}'s parent is root".format(node.name))
                # else:

            (aa, gene_v, freq, isotype, *seqNum) = node.name.split("_")
            node.add_feature("freq",int(freq))
            node.add_feature("aa", aa)
            node.add_feature("isotype", isotype)
            node.add_feature("hz_width", math.log2(int(freq)))
            node.add_feature("in_GMA_path",False)


            if node.isotype == "IGHGP":
                print("tree {} has IGHGP".format(tree_file))
                node.delete()

            if len(node.aa) > max_aa_len:
                max_aa_len = len(node.aa)

            if len(node.sequence) > max_nt_len:
                max_nt_len = len(node.sequence)

        node.img_style["hz_line_color"] = "#DDDDDD"
        node.img_style["hz_line_width"] = 2
        node.img_style["vt_line_width"] = 2
        node.img_style["vt_line_color"] = "#DDDDDD"



    legend_size_count = 0
    for size in [10000, 1000, 100, 10]:
        site_num = legend_size_count * 2
        ts.legend.add_face(CircleFace(math.log2(size), color = "#DCDCDC"), column = site_num)
        ts.legend.add_face(TextFace(str(size).ljust(8), bold=True), column = (site_num+1))
        legend_size_count = legend_size_count + 1

    legend_isotpye_count = 0
    for isotype in ['IGHM', 'IGHD','IGHG3','IGHG1', 'IGHG2', 'IGHG4','IGHA1','IGHA2', 'None']:
        site_num = legend_isotpye_count * 2
        if legend_isotpye_count >=legend_size_count:
            ts.legend.add_face(TextFace(" ",  bold=True),  column=(site_num))
            ts.legend.add_face(TextFace(8 * " ", bold=True), column=(site_num + 1))
        ts.legend.add_face(CircleFace(2, color = isotype_color[isotype]), column = site_num)
        printed_isotype = "Others" if isotype == "None" else isotype
        ts.legend.add_face(TextFace(printed_isotype.ljust(8), fgcolor = isotype_color[isotype], bold=True), column = (site_num+1))
        legend_isotpye_count = legend_isotpye_count + 1
    mutation_rank[cluster] =  nor_mutations_sum
    if nor_mutations_sum ==0:
        print("take care of cluster {}".format(cluster))

    T.render("../data/tree/cluster{}.png".format(str(cluster).zfill(5)), tree_style=ts, dpi=300)




# Figure 2
    draw_flag = False
    for node in t.traverse():
        if node.is_leaf():
            if node.isotype == "IGHM":
                childs = node.up.get_leaves()
                childs_with_GA = []
                for child in childs:
                    if child.isotype.startswith("IGHG") or child.isotype.startswith("IGHA"):
                        draw_flag = True
                        node.up.img_style["vt_line_color"] = "#FF0000"
                        #
                        node.img_style["hz_line_color"] = "#FF0000"
                        node_tmp = child
                        while node_tmp != node.up:
                            #
                            node_tmp.img_style["hz_line_color"] = "#FF0000"
                            node_tmp.in_GMA_path = True
                            node_tmp = node_tmp.up
        else:
            freq = 0
            hz_width = 0
            leaves = node.get_leaves()
            for leaf in leaves:
                if hasattr(leaf,"freq"):
                    freq = freq + leaf.freq
                # hz_width = hz_width + leaf.hz_width
            node.add_feature("freq",freq)
            node.add_feature("hz_width", math.log2(freq))

    if draw_flag:
        T.render("../data/tree/Tree_with_M-G_A_order/cluster{}.png".format(str(cluster).zfill(5)), tree_style=ts, dpi = 300)


sorted_rank = sorted(mutation_rank.items(), key=lambda item: item[1])

with open("../data/tree/mutationRank.txt","w") as f:
    for element in sorted_rank:
        f.write("{} : {}\n".format(*element))