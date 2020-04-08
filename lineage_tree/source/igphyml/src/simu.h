/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PHYML :  a program that  computes maximum likelihood  phylogenies from
DNA or AA homologous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#ifndef CURR_H
#define CURR_H

#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "models.h"
#include "free.h"
#include "spr.h"
#include "pars.h"
#include "alrt.h"

void Simu_Loop(option* io);
int Simu(t_tree *tree,int n_step_max);
void Select_Edges_To_Swap(t_tree *tree,t_edge **sorted_b,int *n_neg);
void Update_Bl(t_tree *tree,phydbl fact);
void Make_N_Swap(t_tree *tree,t_edge **b,int beg,int end);
int Make_Best_Swap(t_tree *tree);
int Mov_Backward_Topo_Bl(t_tree *tree,phydbl lk_old,t_edge **tested_b,int n_tested);
void Unswap_N_Branch(t_tree *tree,t_edge **b,int beg,int end);
void Swap_N_Branch(t_tree *tree,t_edge **b,int beg,int end);
void Check_NNI_Scores_Around(t_node *a, t_node *d, t_edge *b, phydbl *best_score);
int Mov_Backward_Topo_Pars(t_tree *tree, int pars_old, t_edge **tested_b, int n_tested);

#endif
