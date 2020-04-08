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

#include "simu.h"
#include "io.h"

/*********************************************************/

void Simu_Loop(option* io){
  phydbl lk_old;
  int i;
  if(io->mod->opt_heuristic_manuel==YES)Lazy_Exit("Heuristic optimisation for topology search",__FILE__,__LINE__);

  io->mod->update_eigen=YES;
  io->both_sides = 0;
  if(io->mod->optDebug)printf("here1\n");
  Lk_rep(io);
  if(io->mod->optDebug)printf("here\n");

  if((io->mod->s_opt->print) && (!io->mod->quiet)) PhyML_Printf("\n. Maximizing likelihood (using NNI moves)...\n");
  
  Print_Lk_rep(io,"[Initial trees and parameters]");
  
  //added by Ken 8/3/2017
  /* Optimise parameters and branch lengths */
  Optimiz_All_Free_Param(io,(io->mod->quiet)?(0):(io->mod->s_opt->print),0);
  io->both_sides=1;
  	  //optimze branch lengths initially
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree)
#endif
  For(i,io->ntrees){
	  t_tree* tree=io->tree_s[i];
  	  int startnode = 0;
	  if(tree->mod->whichrealmodel<=HLP17){
		  startnode=tree->mod->startnode;
		  tree->both_sides = 1;
		  Lk(tree);
		  Update_Ancestors_Edge(tree->noeud[tree->mod->startnode],tree->noeud[tree->mod->startnode]->v[0],tree->noeud[tree->mod->startnode]->b[0],tree); //added by Ken 7/11
		  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
	  }
	  Optimize_Br_Len_Serie(tree->noeud[startnode],
			tree->noeud[startnode]->v[0],
			tree->noeud[startnode]->b[0],
			tree,
			tree->data);
	  tree->both_sides = 1;
	  Lk(tree);
	  if(tree->mod->whichrealmodel<=HLP17){
		  Get_UPP(tree->noeud[startnode], tree->noeud[startnode]->v[0], tree);
	  }
	  if(tree->mod->print_trace){
		 // Print_Trace(tree);
	  }
  }

  //Start NNI loop while optimizing substitution parameters
  int j;
  int first=0;
  do{
	  lk_old = io->replnL;
	  if(first != 0)Optimiz_All_Free_Param(io,(io->mod->quiet)?(0):(io->mod->s_opt->print),0);
	  first++;
	  io->threads=0;
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree) private(j)
#endif
	  For(i,io->ntrees){//do this in parallel
		  t_tree* tree;
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
		  {
    		  tree=io->tree_s[io->threads++];
    		 // printf("\nSimu Now on %d %d",io->threads,tree->mod->num);
    	 }
	     if(!Simu(tree,10)){
	    	  Check_NNI_Five_Branches(tree);
	     }
	  }
  	  io->replnL=0.0;
  	  For(i,io->ntrees)io->replnL+=io->tree_s[i]->c_lnL; //CHANGE BACK - JUST USED FOR DEBUGGING
  }while(io->replnL > lk_old + io->mod->s_opt->min_diff_lk_global);

  //Start NNI loop without optimizing substitution parameters
  io->threads=0;
#if defined OMP || defined BLAS_OMP
#pragma omp parallel for if(io->splitByTree) private(j)
#endif
  For(i,io->ntrees){
  	  t_tree* tree;// = io->tree_s[i];
#if defined OMP || defined BLAS_OMP
#pragma omp critical
#endif
    	 {
    		  tree=io->tree_s[io->threads++];
    		  //printf("\nNNI Now on %d %d",io->threads,tree->mod->num);
    	 }
    	 do{
    		 if(!Check_NNI_Five_Branches(tree)) break;
    	 }while(1);
    }
}

/*********************************************************/

int Simu(t_tree *tree, int n_step_max)
{
  phydbl old_loglk,n_iter,lambda;
  int i,n_neg,n_tested,n_without_swap,n_tot_swap,step,it_lim_without_swap;
  t_edge **sorted_b,**tested_b;
  int opt_free_param;
  int recurr;

  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  tested_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  
  old_loglk           = UNLIKELY;
  tree->c_lnL         = UNLIKELY;
  n_iter              = 1.0;
  it_lim_without_swap = (tree->mod->invar)?(8):(5);
  n_tested            = 0;
  n_without_swap      = 0;
  step                = 0;
  lambda              = .75;
  n_tot_swap          = 0;
  opt_free_param      = 0;
  recurr              = 0;
  Update_Dirs(tree);
      
  if(tree->lock_topo){
      PhyML_Printf("\n. The tree topology is locked.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
  }
  do{
      ++step;
      old_loglk = tree->c_lnL;	    
      tree->both_sides = 1;
      Lk(tree);
      if(tree->c_lnL < old_loglk){
    	  if((tree->mod->s_opt->print) && (!tree->mod->quiet)) PhyML_Printf("\n\n. Moving backward\n");
    	  if(!Mov_Backward_Topo_Bl(tree,old_loglk,tested_b,n_tested))
    		  Exit("\n. Err: mov_back failed\n");
    	  if(!tree->n_swap) n_neg = 0;
	  
    	  For(i,2*tree->n_otu-3) tree->t_edges[i]->l_old = tree->t_edges[i]->l;
    	  tree->both_sides = 1;
    	  Lk(tree);
      }

      if(step > n_step_max) break;

      /*if(tree->mod->print_trace){
    	  PhyML_Fprintf(tree->mod->fp_out_trace,"[%f]%s\n",tree->c_lnL,Write_Tree(tree)); fflush(tree->mod->fp_out_trace);
    	  if(tree->mod->print_site_lnl) Print_Site_Lk(tree,tree->mod->fp_out_lk); fflush(tree->mod->fp_out_lk);
    	  Print_Trace(tree);
      }*/

      if((tree->mod->s_opt->print) && (!tree->mod->quiet)) Print_Lk(tree,"[Topology           ]");

      /*********************************************************/
      if((FABS(old_loglk-tree->c_lnL) < tree->mod->s_opt->min_diff_lk_global) ||
    		  (n_without_swap > it_lim_without_swap)) break;

      if(tree->mod->whichrealmodel <= HLP17)
    	  Get_UPP(tree->noeud[tree->mod->startnode], tree->noeud[tree->mod->startnode]->v[0], tree);

      Fill_Dir_Table(tree);
      Fix_All(tree);
      n_neg = 0;
      For(i,2*tree->n_otu-3){
      	  if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax)){
      		 NNI(tree,tree->t_edges[i],0);
      	  }
      }
      
      Select_Edges_To_Swap(tree,sorted_b,&n_neg);	    	  
      Sort_Edges_NNI_Score(tree,sorted_b,n_neg);	    
      Optimiz_Ext_Br(tree);	  	    
      Update_Bl(tree,lambda);
	  	        
      n_tested = 0;
      For(i,(int)CEIL((phydbl)n_neg*(lambda)))
      	  tested_b[n_tested++] = sorted_b[i];

      Make_N_Swap(tree,tested_b,0,n_tested);
      n_tot_swap += n_tested;
      
      if(n_tested > 0) n_without_swap = 0;
      else             n_without_swap++;
      
      n_iter+=1.0;
    }
  while(1);
    
/*   Round_Optimize(tree,tree->data); */

  free(sorted_b);
  free(tested_b);

  return n_tested;
}

/*********************************************************/

void Select_Edges_To_Swap(t_tree *tree, t_edge **sorted_b, int *n_neg)
{
  int i;
  t_edge *b;
  int min;
  phydbl best_score;

  *n_neg = 0;
  min = 0;
  
  For(i,2*tree->n_otu-3){
      b = tree->t_edges[i];
      best_score = b->nni->score;

      if((!b->left->tax) && (!b->rght->tax) && (b->nni->score < -tree->mod->s_opt->min_diff_lk_move)){
    	  Check_NNI_Scores_Around(b->left,b->rght,b,&best_score);
    	  Check_NNI_Scores_Around(b->rght,b->left,b,&best_score);
    	  if(best_score < b->nni->score) continue;
    	  sorted_b[*n_neg] = b;
    	  (*n_neg)++;
      }
  }
}

/*********************************************************/

void Update_Bl(t_tree *tree, phydbl fact)
{
  int i;
  t_edge *b;

  For(i,2*tree->n_otu-3){
      b = tree->t_edges[i];
      b->l = b->l_old + (b->nni->l0 - b->l_old)*fact;
    }
}

/*********************************************************/

void Make_N_Swap(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  int dim;

  dim = 2*tree->n_otu-2;

  tree->n_swap = 0;
  for(i=beg;i<end;i++){
      /* we use t_dir here to take into account previous modifications of the topology */
      Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
	   b[i]->nni->swap_node_v2,
	   b[i]->nni->swap_node_v3,
	   b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
	   tree);
      b[i]->l = b[i]->nni->best_l;
      tree->n_swap++;
    }
}

/*********************************************************/

int Make_Best_Swap(t_tree *tree)
{
  int i,j,return_value;
  t_edge *b,**sorted_b;
  int dim;

  dim = 2*tree->n_otu-2;

  sorted_b = (t_edge **)mCalloc(tree->n_otu-3,sizeof(t_edge *));
  
  j=0;
  For(i,2*tree->n_otu-3) if((!tree->t_edges[i]->left->tax) &&
			    (!tree->t_edges[i]->rght->tax))
                             sorted_b[j++] = tree->t_edges[i];

  Sort_Edges_NNI_Score(tree,sorted_b,tree->n_otu-3);

  if(sorted_b[0]->nni->score < -0.0){
      b = sorted_b[0];
      return_value = 1;
      Swap(b->nni->swap_node_v2->v[tree->t_dir[b->nni->swap_node_v2->num*dim+b->nni->swap_node_v1->num]],
	   b->nni->swap_node_v2,
	   b->nni->swap_node_v3,
	   b->nni->swap_node_v3->v[tree->t_dir[b->nni->swap_node_v3->num*dim+b->nni->swap_node_v4->num]],
	   tree);
      b->l = b->nni->best_l;
  }
  else return_value = 0;

  free(sorted_b);

  return return_value;
}

/*********************************************************/

int Mov_Backward_Topo_Bl(t_tree *tree, phydbl lk_old, t_edge **tested_b, int n_tested)
{
  phydbl *l_init;
  int i,step,beg,end;
  t_edge *b;


  l_init = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));

  For(i,2*tree->n_otu-3) l_init[i] = tree->t_edges[i]->l;
  
  step = 2;
  tree->both_sides = 0;
  do{
      For(i,2*tree->n_otu-3){
    	  b = tree->t_edges[i];
	  	  b->l = b->l_old + (1./step) * (l_init[i] - b->l_old);
      }

      beg = (int)FLOOR((phydbl)n_tested/(step-1));
      end = 0;
      Unswap_N_Branch(tree,tested_b,beg,end);
      beg = 0;
      end = (int)FLOOR((phydbl)n_tested/step);
      Swap_N_Branch(tree,tested_b,beg,end);
      
      if(!end) tree->n_swap = 0;
      
      tree->both_sides = 0;
      Lk(tree);
      
      step++;

  }while((tree->c_lnL < lk_old) && (step < 1000));

  
  if(step == 1000){
      if(tree->n_swap)  Exit("\n. Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");
      For(i,2*tree->n_otu-3){
    	  b = tree->t_edges[i];
    	  b->l = b->l_old;
      }
      tree->both_sides = 0;
      Lk(tree);
    }

  free(l_init);

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3){
      if(tree->t_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->t_edges[i]->nni->score = +1.0;
  }
  
  if(tree->c_lnL > lk_old)                                return  1;
  else if((tree->c_lnL > lk_old-tree->mod->s_opt->min_diff_lk_local) && 
	  (tree->c_lnL < lk_old+tree->mod->s_opt->min_diff_lk_local)) return -1;
  else                                                    return  0;
}

/*********************************************************/

int Mov_Backward_Topo_Pars(t_tree *tree, int pars_old,
		t_edge **tested_b, int n_tested){
  int i,step,beg,end;
  
  step = 2;
  tree->both_sides = 0;
  do{
    beg = (int)FLOOR((phydbl)n_tested/(step-1));
    end = 0;
    Unswap_N_Branch(tree,tested_b,beg,end);
    beg = 0;
    end = (int)FLOOR((phydbl)n_tested/step);
    Swap_N_Branch(tree,tested_b,beg,end);
      
    if(!end) tree->n_swap = 0;
      
    tree->both_sides         = 0;
    Pars(tree);
      
    step++;

  }while((tree->c_pars > pars_old) && (step < 1000));

  
  if(step == 1000){
      if(tree->n_swap)Exit("\n. Err. in Mov_Backward_Topo_Bl (n_swap > 0)\n");

      tree->both_sides = 0;
      Pars(tree);
  }

  tree->n_swap = 0;
  For(i,2*tree->n_otu-3){
      if(tree->t_edges[i]->nni->score < 0.0) tree->n_swap++;
      tree->t_edges[i]->nni->score = +1.0;
  }

  
  if(tree->c_pars < pars_old)       return  1;
  else if(tree->c_pars == pars_old) return -1;
  else                              return  0;
}

/*********************************************************/

void Unswap_N_Branch(t_tree *tree, t_edge **b, int beg, int end){
  int i;
  int dim;

  dim = 2*tree->n_otu-2;

  if(end>beg){
      for(i=beg;i<end;i++){

    	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
	       tree);

    	  b[i]->l = b[i]->l_old;
      }
  }else{
      for(i=beg-1;i>=end;i--){
    	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
	       tree);

	  b[i]->l = b[i]->l_old;
	}
  }
}

/*********************************************************/

void Swap_N_Branch(t_tree *tree,t_edge **b, int beg, int end)
{
  int i;
  int dim;
  dim = 2*tree->n_otu-2;
  if(end>beg){
      for(i=beg;i<end;i++){
    	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
	       tree);
    	  b[i]->l = b[i]->nni->best_l;
      }
  }else{
      for(i=beg-1;i>=end;i--){
    	  Swap(b[i]->nni->swap_node_v2->v[tree->t_dir[b[i]->nni->swap_node_v2->num*dim+b[i]->nni->swap_node_v1->num]],
	       b[i]->nni->swap_node_v2,
	       b[i]->nni->swap_node_v3,
	       b[i]->nni->swap_node_v3->v[tree->t_dir[b[i]->nni->swap_node_v3->num*dim+b[i]->nni->swap_node_v4->num]],
	       tree);
    	  b[i]->l = b[i]->nni->best_l;
      }
  }
}

/*********************************************************/

void Check_NNI_Scores_Around(t_node *a, t_node *d, t_edge *b, phydbl *best_score){
  int i;
  For(i,3){
      if((d->v[i] != a) && (!d->v[i]->tax)){
    	 if((d->b[i]->nni->score > *best_score-1.E-10) &&
    		(d->b[i]->nni->score < *best_score+1.E-10)){
	           d->b[i]->nni->score = *best_score+1.;
    	 }
    	 if(d->b[i]->nni->score < *best_score){
	      *best_score = d->b[i]->nni->score;
	    }
      }
  }
}

/*********************************************************/
