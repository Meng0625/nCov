/*
IgPhyML: a program that computes maximum likelihood phylogenies under
non-reversible codon models designed for antibody lineages.

Copyright (C) Kenneth B Hoehn. Sept 2016 onward.

built upon

codonPHYML: a program that  computes maximum likelihood phylogenies from
CODON homologous sequences.

Copyright (C) Marcelo Serrano Zanetti. Oct 2010 onward.

built upon

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

** spr.h: Header file for the SPR routines.
**
** Wim Hordijk   Last modified: 28 August 2006
*/

#ifndef _SPR_H_
#define _SPR_H_

#include "utilities.h"
#include "models.h"
#include "lk.h"
#include "free.h"
#include "optimiz.h"
#include "alrt.h"
#include "pars.h"
#include "simu.h"

#define ALL   1
#define BEST  2
#define ONE   3

/*
** _move_: Structure for holding the relevant information for candidate SPR moves.
*/

typedef struct
{
  t_node   *v_prune, *u_prune, *v_n, *v_nx1, *u_n, **path;
  t_edge   *e_prune, *e_regraft;
  phydbl  l_connect, l_est[3], delta_lk, d_L, d_up_v, d_un_v;
  int     dist, rgrft_rank, optim_rank, globl_rank;
} _move_;



void Init_SPR          (t_tree *tree);
void Clean_SPR         (t_tree *tree);
void Optim_SPR         (t_tree *tree, int max_size, int method);
int  Perform_SPR_Moves (t_tree *tree, int max_size);
int  Perform_Best_SPR  (t_tree *tree, int max_size);
int  Perform_One_SPR   (t_tree *tree, int max_size);

void Calc_Tree_Length (t_edge *e_prune, t_node *v_prune, t_tree *tree);
void Tree_Length      (t_node *v_prune, t_node *u_prune, t_node *v_n, t_node *v_n_1,
		       t_node *v_nx1, t_node *v_0, t_node *u_n, phydbl d_up_v_1,
		       phydbl d_uu, phydbl d_L_1, int n, t_tree *tree);
int  Est_Lk_Change    (t_edge *e_prune, t_node *v_prune, t_tree *tree);
int  Best_Lk_Change   (t_edge *e_prune, t_node *v_prune, t_tree *tree);
void Make_Move        (_move_ *move, int type, t_tree *tree);
int  Find_Optim_Local (t_tree *tree);
int  Find_Optim_Globl (t_tree *tree);
void Prune            (t_edge *e, t_node *v, t_edge **e_connect, t_edge **e_avail,
		       t_tree *tree);
void Regraft          (t_edge *e, t_node *v, t_edge *avail, t_tree *tree);
void PostOrder_v      (t_tree *tree, t_node *v, t_edge *e);
void PostOrder_w      (t_tree *tree, t_node *v, t_edge *v_e, t_node *w, t_edge *e);





void Speed_Spr(t_tree *tree, int max_cycles);
//void Speed_Spr_Loop(t_tree* tree);
void Speed_Spr_Loop(option* io);
void Make_Spr_List(t_tree *tree);
void Init_One_Spr(spr *a_spr);
spr *Make_One_Spr(t_tree *tree);
int Spr(phydbl init_lnL, t_tree *tree);
int Spr_Recur(t_node *a, t_node *d, t_tree *tree);
int Test_All_Spr_Targets(t_edge *pulled, t_node *link, t_tree *tree);
void Randomize_Spr_List(t_tree *tree);
void Test_One_Spr_Target_Recur(t_node *a, t_node *d, t_edge *pulled, t_node *link, t_edge *residual, int *best_found, t_tree *tree);
phydbl Test_One_Spr_Target(t_edge *target, t_edge *arrow, t_node *link, t_edge *residual, t_tree *tree);
int Try_One_Spr_Move_Triple(spr *move, t_tree *tree);
int Try_One_Spr_Move_Full(spr *move, t_tree *tree);
void Make_Best_Spr(t_tree *tree);
void Include_One_Spr_To_List_Of_Spr(spr *move, t_tree *tree);
void Reset_Spr_List(t_tree *tree);
int Evaluate_List_Of_Regraft_Pos_Triple(spr **spr_list, int list_size, t_tree *tree);
void Best_Spr(t_tree *tree);
int Check_Spr_Move_Validity(spr *this_spr_move, t_tree *tree);
void Spr_Subtree(t_edge *b, t_node *link, t_tree *tree);
void Spr_Pars(t_tree *tree);
void Print_Trace(t_tree *tree); //Added by Ken 9/2/2017


#endif  /* _SPR_H_ */


/*
** EOF: spr.h
*/
