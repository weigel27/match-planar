/************************************************************************/
/*                                                                      */
/*             ROUTINE FOR WEIGHTED PERFECT MATCHING PROBLEMS           */
/*                                                                      */
/*  Written by:  A. Rohe                                                */
/*  Date:  November 7, 1995  (rohe)                                     */
/*         November 12, 1995 (rohe)                                     */
/*         February 7, 1996  (rohe - short look at program)             */
/*         July, 1996        (rohe - change from tree to forest)        */
/*         October, 1996     (bico - change to concorde format)         */
/*                                                                      */
/*                                                                      */
/* EXPORTED FUNCTION:                                                   */
/*   int perfect_match (int ncount, CCdatagroup *dat, int ecount,       */
/*          int **elist, int **elen, char *blo_filename,                */
/*          char *mat_filename, int just_fractional, int no_fractional, */
/*          int use_all_trees, int use_tree, int partialprice,          */
/*          double *totalzeit)                                          */
/*     COMPUTES a minumum weight perfect matching.                      */
/*         -ncount is the number of nodes in the graph                  */
/*         -dat gives the data to generate edge lengths (if the         */
/*            matching is over the complete graph), it can be NULL      */
/*         -ecount is the number of edges in the initial edge set       */
/*         -elist is a pointer to an array of edges in end end format,  */
/*            this array will be freed by perfect_match to save on      */
/*            memory                                                    */
/*         -elen is a pointer to an array of edge lengths for the edges */
/*            in elist, it will be freed by the perfect_match           */
/*         -blo_filename is the name of a file (it can be NULL), if it  */
/*            is not NULL then the blossoms and weights will be written */
/*            to it                                                     */
/*         -mat_filename is the name of a file (it can be NULL), if it  */
/*            is not null, then the matching will be writtien to it     */
/*         -just_fractional should be 0 to compute a matching, and 1    */
/*            to only compute a fractional matching                     */
/*         -if no_fractional is nonzero, then a fractional jumpstart    */
/*            will not be used                                          */
/*         -if use_all_trees is 1, then a common delta is               */
/*            computed for all trees in the forest; if it is 2, then    */
/*            a single tree is grown                                    */
/*         -if partialprice is > 0, then the pricing routine will first */
/*            price the partialprice nearest neighbors before moving to */
/*            the complete edge set                                     */
/*         -totalzeit will return the time used by the matching routine */
/*            (excluding the time needed for a final check that the     */
/*            matching is indeed optimal, if this checking is turned    */
/*            on)                                                       */
/*                                                                      */
/*  NOTES:                                                              */
/*                                                                      */
/*    Returns 0 if it worked and 1 otherwise (for example, when one     */
/*    of the mallocs failed). The nodes in the graph should be named    */
/*    0 through #nodes - 1.                                             */
/*                                                                      */
/************************************************************************/

#define NDEBUG   /* remove this to turn assertions on */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "machine_specific.h"

#define CARELESS_REPAIRS 
/* #define BALL_DERIGS_REPAIRS */
#define GREEDY_DUAL_CHANGE 
#define SWITCH_LEVEL 32
#define MAX_BAD_PORTION 0.25

#define PRINT_LEVEL 0
#define SHRINKS_SIZE 1000
/* infinity */
#define INTEGER_MIN -999999999
#define INTEGER_MAX  999999999
/* for labels */
#define NIX 0
#define PLUS 1
#define MINUS 2
/* for status */
#define UNMATCHED 0
#define HALVES 1
#define MATCHED 2
/* for x */
#define HALF 2

#define CC_SAFE_MALLOC(nnum,type)                                          \
    (type *) malloc (((unsigned int) (nnum)) * sizeof (type))

#define CC_FREE(object,type) {                                             \
    free ((void *) (object));                                              \
    object = (type *) NULL;                                                \
}

#define CC_IFFREE(object,type) {                                           \
    if ((object)) CC_FREE ((object),type);                                 \
}

/*
  #define OTHEREND(e,n,list)  ((list) + (e)->nod1 + (e)->nod2 - ((n) - (list)))
*/
#define OTHEREND(e,n,list) (((e)->nod1 + (list)) == (n) ?                \
                             (e)->nod2 + (list) : (e)->nod1 + (list))
#define OTHEREND_INT(e,n) ((e)->nod1 + (e)->nod2 - (n))

#define FIND_SURFACE(n) {                                                \
            while((n)->blossom_parent != -1) {                           \
                (n)=G->nodelist+(n)->blossom_parent;                     \
            }                                                            \
        }

#define SHRINK_SORT_CONST 3
#define BADCOUNT_CONST 10000
#define PNODE(n) (G->nodelist + (n))
#define PEDGE(e) (G->edgelist + (e))

#define INODE(n) ((n) - G->nodelist)
#define IEDGE(e) ((e) - G->edgelist)

typedef struct edge {
  int   slack;
  char  mark;  
  char  x; 
  int   ptrs[2];
  int   nod1;
  int   nod2;
  int   orig_nod1;
  int   orig_nod2;
} edge;

typedef struct node {
  int   edg_list;
  int   matched_edge;

  /* These are for the augmenting path tree */
  int   child;
  int   sibling;
  int   parent;
  int   parent_edge;

  /* These are for the blossom_tree */
  int   blossom_next;
  int   blossom_parent;
  int   blossom_root;
  int   blossom_next_edge;
  int   penultimate;

  int   pi;
  int   mark;      /* unused pseudo nodes are linked by mark */
  int   tree_root;

  char  status;
  char  label;
  char  hit;

  int   dummy;    /* just to pad to multiple of 8 bytes */
} node;

typedef struct nodeptr {
  int next;
  int surf;
  int delta;
  int dad;
  int sum;
  int status;
#ifdef GREEDY_DUAL_CHANGE
  char tree_processed;
  int gnext;      /* used to pad to multiple of 8 bytes and in reordering */
#endif
} nodeptr;

typedef struct shrink_ary_T {
  int node_i;
  int node_j;
  int node_k;
  int edge;
  int size;
  int next;
} shrink_ary_T;

typedef struct shrink_T {
  shrink_ary_T *ary;
  int length;
} shrink_T;

typedef struct expand_T {
  int *node;
  int length;
} expand_T;

typedef struct stats {
  int expand_count;
  int shrink_count;
  int dualchange_count;
  int dualzero_count;
} stats;

typedef struct graph {
  edge *edgelist;
  node *nodelist;
  int nnodes;
  int nedges;
  int max_nedges;
  int unused;
  int unmatched;
  nodeptr *roots;
  int unmatched_count;
}graph;

typedef struct srkstuff {
  expand_T expands;
  shrink_T shrinks;
  int shrinks_max;
  int expands_max;
  int shrinks_[SHRINKS_SIZE];
} srkstuff;

typedef struct stack {
  int*             ary;
  int              count;
} stack;

#if PRINT_LEVEL
static graph *PG = (graph *) NULL;
#endif

static void
adjust_match (graph *G),
  print_node (graph *G, node *n),
  init_tree (graph *G, node *n),
  clear_tree (graph *G, node *n),
  make_cycle_halves (graph *G, node *node_i, node *node_j, node *node_k,
		     edge *e),
  shrink_tree (graph *G, node *newnode),
  label_penultimate (graph *G, node *n, int label),
  fix_matching (graph *G, node *oldnode, node *x),
  fix_tree (graph *G, node *oldnode, node *x, node *y),
  flip_cycle (graph *G, node *node_i),
  augment_path (graph *G, node *node_i, int fractional),
  vereinigung (graph *G, int x, int y),
  set_vereinigung (graph *G, node *n),
  dual_match_len (graph *G, int fractional, double *val),
  fix_match (graph *G, int blossom),
  mark_match (graph *G, float **dists);

static int 
build_graph (graph *G, int ncount, int ecount, float **dists),
  init (graph *G, srkstuff *srk),
  match_main_frac (graph *G, stats *scount, srkstuff *srk),
  match_main (graph *G, stats *scount, srkstuff *srk, int use_all),
  augment_blossom (graph *G, edge *, int fractional, srkstuff *srk),
  lower_edges (graph *G, node *old),
  expand_blossom (graph *G, node *oldnode, stats *scount),
  lift_edges (graph *G, node *newnode),
  add_to_tree (graph *G, edge *e, node *node_i, node *node_j, int fractional),
  checkout_node (graph *G, node *n, int fractional, srkstuff *srk),
  grow_tree_no_shrink (graph *G, node *n, int fractional, srkstuff *srk),
  grow_tree (graph *G, node *n, int fractional, stats *scount, srkstuff *srk,
	     int *found_aug),
  apply_dual_change (graph *G, node *n, int delta),
  find_parity_sum (graph *G, int n),
  parity_correction (graph *G, stats *scount),
  find_single_dual_change (graph *G, node *n),
  find_dual_change_forest (graph *G, node *n),
  make_dual_change_forest (graph *G, stats *scount),
  match_main_more_in_forest (graph *G, stats *scount, srkstuff *srk),
  match_main_forest (graph *G, stats *scount, srkstuff *srk),
  match_main_tree (graph *G, stats *scount, srkstuff *srk),
  make_match (graph *G),
  write_match (graph *G, float **dists, char* match_file);

static node
*shrink_blossom (graph *G, node *node_i, node *node_j, node *node_k,
		 edge *e, stats *scount),
  *common_parent (graph *G, node *node_i, node *node_j, int *size),
  *find_below (graph *G, node *n, int blossom);

int perfect_match (int ncount, int ecount, float **dists,
                   int just_fractional, int no_fractional,
                   int use_all_trees, int partialprice)
{
  graph G;
  srkstuff srk;
  stats scount;
  double val;
  int finished = 1;
  int sedges = 0;

#if PRINT_LEVEL
  PG = &G;
#endif

  G.edgelist = (edge *) NULL;
  G.nodelist = (node *) NULL;
  G.unused = -1;
  G.unmatched = -1;
  G.roots = (nodeptr *) NULL;

  srk.expands.node = (int *) NULL;
  srk.shrinks.ary = (shrink_ary_T *) NULL;

  if (build_graph (&G, ncount, ecount, dists)) {
    fprintf (stderr, "build_graph failed\n");
    /*
    CC_FREE (*elist, int);
    CC_FREE (*elen, int);
    */
    goto CLEANUP;
  }
  /*CC_FREE (*elist, int);*/
  sedges = G.nedges;
  
  do {
#ifndef NDEBUG
    printf (" Starting Init...");
    fflush (stdout);
#endif
    if (init (&G, &srk)) {
      fprintf (stderr, "init failed\n");
      goto CLEANUP;
    }

    if (no_fractional == 0) {
#ifndef NDEBUG
      printf (" Matching with %i nodes and %i edges ...\n",
	      G.nnodes, G.nedges);
      fflush (stdout);
      printf (" Make a fractional Matching...\n");
      fflush (stdout);
#endif
      scount.expand_count = 0;
      scount.shrink_count = 0;
      scount.dualchange_count = 0;
      scount.dualzero_count = 0;
      if (match_main_frac (&G, &scount, &srk)) {
	fprintf (stderr, "match_main_frac failed\n");
	goto CLEANUP;
      }
    }

    if (just_fractional == 0) {
#ifndef NDEBUG
      printf (" Make an integral Matching ...\n");
      fflush (stdout);
#endif
      if (make_match (&G)) {
	fprintf (stderr, "make_match failed\n");
	fflush (stdout);
      }

#ifndef NDEBUG
      printf (" Starting Integral Matching-Code...\n");
      fflush (stdout);
#endif
      scount.expand_count = 0;
      scount.shrink_count = 0;
      scount.dualchange_count = 0;
      scount.dualzero_count = 0;
      if (match_main (&G, &scount, &srk, use_all_trees)) {
	fprintf (stderr, "match_main failed\n");
	goto CLEANUP;
      }
      dual_match_len (&G, 0, &val);
#ifndef NDEBUG
      printf ("     Dual Matching Length: %.1f\n", val);
      fflush (stdout);
#endif
      CC_IFFREE (G.roots, nodeptr);
    }
    CC_IFFREE (srk.expands.node, int);
    CC_IFFREE (srk.shrinks.ary, shrink_ary_T); 
  } while (!finished);

#ifndef NDEBUG
  printf (" Adjust found Matching...\n");
  fflush (stdout);
#endif
  adjust_match (&G);

#ifndef NDEBUG
  printf ("ZZ Edges:  %d (starting)   %d (total)\n", sedges, G.nedges);
#endif

  mark_match(&G, dists);

 CLEANUP:

  CC_IFFREE (G.nodelist, node);
  CC_IFFREE (G.edgelist, edge);
  CC_IFFREE (G.roots, nodeptr);
  CC_IFFREE (srk.expands.node, int);
  CC_IFFREE (srk.shrinks.ary, shrink_ary_T); 

  return 0;
}

static void print_node (graph *G, node *n)
{
  edge *e;
  int ei;
  int c;
  node *nodelist = G->nodelist;
  edge *edgelist = G->edgelist;

  printf ("Node %d, Pi %d", (int) INODE(n), n->pi);
  if (n->label == PLUS) printf (" label +");
  if (n->label == MINUS) printf (" label -");
  if (n->label == NIX) printf (" label 0");
  if (n->status == MATCHED)
    printf (" MATCHED (%d,%d)\n",
            edgelist[n->matched_edge].nod1,
            edgelist[n->matched_edge].nod2);
  if (n->status == UNMATCHED) printf (" UNMATCHED\n");
  if (n->status == HALVES)
    printf (" HALVES (%i,%i)\n",
            edgelist[n->matched_edge].nod1, edgelist[n->matched_edge].nod2);
  fflush (stdout);
  if (n->label != NIX) {
    printf ("TreeRoot %d ", n->tree_root);
    if (n->parent)
      printf ("TreeParent = %d ", n->parent);
    printf ("Children: ");
    for (c = n->child; c != -1; c = nodelist[c].sibling)
      printf ("%d ", c);
    printf ("\n");
  }
  if (n->blossom_parent != -1) {
    printf ("BlossomNode %i, Parent %i, Next %i\n", (int) INODE(n),
            n->blossom_parent, n->blossom_next);
  }
  if (n->blossom_root != -1) {
    printf ("BlossomParent %i, Root %i\n", (int) INODE(n), n->blossom_root);
  }
  for (ei = n->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
    e = PEDGE(ei/2);
    printf ("| Edge (%i,%i)", e->nod1, e->nod2);
    printf ("| %4i is ", e->nod1);
    if (nodelist[e->nod1].label == PLUS) printf ("+");
    if (nodelist[e->nod1].label == MINUS) printf ("-");
    if (nodelist[e->nod1].label == NIX) printf ("0");
    printf ("| %4i is ", e->nod2);
    if (nodelist[e->nod2].label == PLUS) printf ("+");
    if (nodelist[e->nod2].label == MINUS) printf ("-");
    if (nodelist[e->nod2].label == NIX) printf ("0");
    printf ("| edg Slack %d, x=%d\n", e->slack, (int) e->x);
  }
}

static int init (graph *G, srkstuff *srk)
{
  int i, count;
  node *n;
  edge *e;
  int ei;
  int delta;
  node *nodelist = G->nodelist;
  int   ncount = G->nnodes;
  int   ecount = G->nedges;

  /* set all the pi s */
  for (i = 0; i < ncount; i++)
    nodelist[i].pi = INTEGER_MAX;
  for (i = ncount; i <= 3 * ncount / 2; i++)
    nodelist[i].pi = 0;

  for (i = 0, e = G->edgelist; i < ecount; i++, e++) {
    if (nodelist[e->nod1].pi > e->slack) {
      nodelist[e->nod1].pi = e->slack;
    }
    if (nodelist[e->nod2].pi > e->slack) {
      nodelist[e->nod2].pi = e->slack;
    }
  }

  for (i = 0; i < ncount; i++)
    nodelist[i].pi /= 2;

  /* calculate the slacks with the pi's and set rest to zero/nothing */

  for (i = 0, e = G->edgelist; i < ecount; i++, e++) {
    e->slack = e->slack - nodelist[e->nod1].pi - nodelist[e->nod2].pi;
    e->x = 0;
  }

  for (i = 0, n = nodelist; i <= 3 * ncount / 2; i++, n++) {
    n->child = -1;
    n->sibling = -1;
    n->parent = -1;
    n->label = NIX;
    n->parent_edge = -1;
    n->matched_edge = -1;
    n->status = UNMATCHED;
    n->blossom_root = -1;
    n->blossom_next = -1;
    n->blossom_next_edge = -1;
    n->blossom_parent = -1;
    n->tree_root = -1;
    n->hit = 0;
  }
#if 0
  /* Take all 0-slack edges directly, old version */
  count = 0;
  for (i = 0, e = G->edgelist; i < ecount; i++, e++) {
    if (e->slack == 0) {
      if ((nodelist[e->nod1].status == UNMATCHED) &&
	  (nodelist[e->nod2].status == UNMATCHED)) {
	e->x = 1;
	nodelist[e->nod1].matched_edge = i;
	nodelist[e->nod2].matched_edge = i;
	nodelist[e->nod1].status = MATCHED;
	nodelist[e->nod2].status = MATCHED;
	count += 2;
      }
    }
  }
#else
  /* Take all 0-slack edges directly & make better pi, new version */
  count = 0;
  for (i = 0, n = nodelist; i < ncount; i++, n++) {
    delta = INTEGER_MAX;
    if (n->status == UNMATCHED) {
      for (ei = n->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	e = PEDGE(ei/2); 
	if (e->slack <= delta) {
	  if (e->slack == 0) {
	    if ((nodelist[e->nod1].status == UNMATCHED) && 
		(nodelist[e->nod2].status == UNMATCHED)) {
	      e->x = 1;
	      nodelist[e->nod1].matched_edge = IEDGE(e);
	      nodelist[e->nod2].matched_edge = IEDGE(e);
	      nodelist[e->nod1].status = MATCHED;
	      nodelist[e->nod2].status = MATCHED;
	      count += 2;
	    }
	  }
	  delta = e->slack;
	}
      }
      if (delta != 0) {
	n->pi += delta;
	for (ei = n->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	  e = PEDGE(ei/2);
	  e->slack -= delta;
	}
      }
    }
  }
#endif
#ifndef NDEBUG
  printf (" %i nodes ", count); fflush (stdout);
#endif

  G->unmatched_count = ncount - count;

  /* changed by MW 25/06/2004 */
  //srk->shrinks_max = ncount;  /* was just ncount */
  srk->shrinks_max = ncount / 5 + 10;  /* was just ncount */
  //srk->shrinks_max = ncount / 10 + 10;  /* was just ncount */
  srk->shrinks.ary = CC_SAFE_MALLOC (srk->shrinks_max, shrink_ary_T);
  srk->expands_max = ncount / 10 + 10;
  srk->expands.node = CC_SAFE_MALLOC (srk->expands_max, int);

  if (!srk->shrinks.ary || !srk->expands.node) {
    fprintf (stderr, "out of memory in init\n");
    CC_IFFREE (srk->shrinks.ary, shrink_ary_T);
    CC_IFFREE (srk->expands.node, int);
    return 1;
  }
  return 0;
}

static void init_tree (graph *G, node *n)
{
#if PRINT_LEVEL
  printf ("  init_tree %i\n", (int) (n - PG->nodelist)); fflush (stdout);
#endif
  n->child = -1;
  n->sibling = -1;
  n->parent = -1;
  n->parent_edge = -1;
  n->label = PLUS;
}

static void clear_tree (graph *G, node *n)
{
  node *c = n;
  node *stop = G->nodelist -1;

  while (1) {
    c->mark = 0;
    c->label = NIX;
    c->tree_root = -1;	/* ??? */

    /* go to next c */
    if (c->child!=-1) {
      c=PNODE(c->child);
    }
    else {
      while (PNODE(c->sibling)==stop) {
	if (c==n) return ; /* this if is for an n without childs */
	c=PNODE(c->parent);
	if (c==n) return ;
      }
      c=PNODE(c->sibling);
    }
  }
}

/* make_cycle_halves - edge from i to j, k is the "root" */
static void make_cycle_halves (graph *G, node *node_i, node *node_j,
                               node *node_k, edge *e)
{
  edge *edgelist = G->edgelist;
  node *parent;

#if PRINT_LEVEL
  printf ("    Make cycle halves: root %i ends %i and %i:",
	  (int) (node_k - PG->nodelist), (int) (node_i - PG->nodelist),
          (int) (node_j - PG->nodelist));
  fflush (stdout);
#endif

  /* set matched_edge for node_j */
  node_j->status = HALVES;
  node_j->matched_edge = e - edgelist;
  e->x = HALF;

  /* path from i to k : matched_edges are the parent_edges */
  for (; node_i != node_k; node_i = PNODE(node_i->parent)) {
    edgelist[node_i->parent_edge].x = HALF;
    node_i->matched_edge = node_i->parent_edge;
    node_i->status = HALVES;
  }

  /* path from j to k : matched_edges are the "child_edges" */
  for (; node_j != node_k; node_j = parent) {
    parent = PNODE(node_j->parent);
    edgelist[node_j->parent_edge].x = HALF;
    parent->matched_edge = node_j->parent_edge;
    parent->status = HALVES;
  }
}

static int lift_edges (graph *G, node *newnode)
{
  node *n, *node_k;
  node *nodelist = G->nodelist;
  int nn = INODE(newnode);
  edge *e;
  int ei, einext;

#if PRINT_LEVEL
  printf ("Lift edges %i\n", nn); fflush (stdout);
#endif

  n = node_k = PNODE(newnode->blossom_root);
  do {
    ei = n->edg_list;
    n->edg_list = -1;
    for (; ei != -1; ei = einext) {
      e = PEDGE(ei/2);
      einext = e->ptrs[ei % 2];
      if ((nodelist[e->nod1].blossom_parent == nn) &&
	  (nodelist[e->nod2].blossom_parent == nn)) {
	e->ptrs[ei % 2] = n->edg_list;
	n->edg_list = ei;
      } else {
	if (e->nod1 == INODE(n)) {    /* nod1 is in new blossom */
	  e->nod1 = nn;
	} else {                      /* nod2 is in new blossom */
	  e->nod2 = nn;
	}
	e->ptrs[ei % 2] = newnode->edg_list;
	newnode->edg_list = ei;
      } 
    }
    n = PNODE(n->blossom_next);
  } while (n != node_k);
  return 0;
}

static void shrink_tree (graph *G, node *newnode)
{
  node *n, *node_k, *cnode;
  int c, csibling, p;

#if PRINT_LEVEL
  printf ("   Shrink-Tree %i\n", (int) (newnode - PG->nodelist));
  fflush (stdout);
#endif

  n = node_k = PNODE(newnode->blossom_root);
  do {
    for (c = n->child; c != -1; c = csibling) {
      cnode = PNODE(c);
      csibling = cnode->sibling;
      if (cnode->blossom_parent != INODE(newnode)) {
	cnode->parent = INODE(newnode);
	cnode->sibling = newnode->child;
	newnode->child = c;
      }
    }
    n = PNODE(n->blossom_next);
  } while (n != node_k);

  /* set other items for newnode */
  p = node_k->parent;

  newnode->sibling = -1;    /* Plus node, has no siblings */
  newnode->parent = p;
  newnode->parent_edge = node_k->parent_edge;
  newnode->matched_edge = node_k->matched_edge;
  newnode->label = PLUS;

  if (p != -1)		/* _p has just one child, so this is ok */
    G->nodelist[p].child = INODE(newnode);
}

/* shrink blossom - edge from i to j, k is the "root" */

static node *shrink_blossom (graph *G, node *node_i, node *node_j,
                             node *node_k, edge *e, stats *scount)
{
  node *newnode;
  node *parent;

#if PRINT_LEVEL
  printf ("    Shrink blossom: root %i ends %i and %i:\n",
	  (int) (node_k - PG->nodelist), (int) (node_i - PG->nodelist),
	  (int) (node_j - PG->nodelist));
  fflush (stdout);
#endif

  scount->shrink_count++;

  newnode = PNODE(G->unused);
  G->unused = newnode->mark;
  newnode->mark = 0;
  newnode->edg_list = -1;

  if (node_k->status == UNMATCHED) {
    G->roots[node_k->tree_root].surf = INODE(newnode);
  }

  /* set blossom_next,next_edge,parent for node_j */
  /* blossom_root has to stay the same */

  node_j->blossom_next = INODE(node_i);
  node_j->blossom_next_edge = IEDGE(e);
  node_j->blossom_parent = INODE(newnode);

  /* path from i to k : blossom_next_edges are the parent_edges */

  for (; node_i != node_k; node_i = PNODE(node_i->parent)) {
    node_i->blossom_next_edge = node_i->parent_edge;
    node_i->blossom_next = node_i->parent;
    node_i->blossom_parent = INODE(newnode);
  }

  /* path from j to k : blossom_next_edges are the "child_edges" */

  for (; node_j != node_k; node_j = parent) {
    parent = PNODE(node_j->parent);
    parent->blossom_next_edge = node_j->parent_edge;
    parent->blossom_next = INODE(node_j);
    parent->blossom_parent = INODE(newnode);
  }

  newnode->blossom_root = INODE(node_k);
  newnode->blossom_next = -1;
  newnode->blossom_next_edge = -1;
  newnode->blossom_parent = -1;

  newnode->mark = 0;
  newnode->pi = 0;
  newnode->status = node_k->status;
  newnode->tree_root = node_k->tree_root;

  shrink_tree (G, newnode);
  if (lift_edges (G, newnode)) {   
    fprintf (stderr, "lift_edges failed\n");
    return (node *) NULL;
  } 

  return newnode;
}

static void label_penultimate (graph *G, node *n, int label)
{
  /* all the leafs of the blossom_tree under n get the penultimate label */
  node *n2;

  if (n->blossom_root == -1) {
    n->penultimate = label;
    return;
  }

  n2 = n;
  while (1) {
    if (n2->blossom_root != -1) { 
      n2 = PNODE(n2->blossom_root);
    } else {
      n2->penultimate = label;
      while (n2->blossom_next == 
	     PNODE(n2->blossom_parent)->blossom_root) {
	n2 = PNODE(n2->blossom_parent);
	if (n2 == n) return;
      }
      n2 = PNODE(n2->blossom_next);
    }
  }
}

static int lower_edges (graph *G, node *old)
{
  edge *e;
  int ei, einext;
  node *x;

#if PRINT_LEVEL
  printf ("Lower edges %i !\n", (int) (old - PG->nodelist));
  fflush (stdout);
#endif

  for (ei = old->edg_list; ei != -1; ei = einext) {
    e = PEDGE(ei/2);
    einext = e->ptrs[ei % 2];
    if (e->nod1 == INODE(old)) {
      x = PNODE(e->orig_nod1);
      e->nod1 = x->penultimate;
    } else {
      x = PNODE(e->orig_nod2);
      e->nod2 = x->penultimate;
    }
    e->ptrs[ei % 2] = PNODE(x->penultimate)->edg_list;
    PNODE(x->penultimate)->edg_list = ei;
  }
  return 0;
}

static void fix_matching (graph *G, node *oldnode, node *x)
{
  node *n, *memo;
  node *nodelist = G->nodelist;
  edge *e;

#if PRINT_LEVEL
  printf ("    Fix Matching %i; x = %i\n",
	  (int) (oldnode - PG->nodelist), (int) (x - PG->nodelist));
  fflush (stdout);
#endif

  if (x->blossom_next_edge == x->matched_edge) {
    n = x;
    memo = PNODE(oldnode->blossom_root);
  } else {
    n = PNODE(oldnode->blossom_root);
    memo = x;
  }

  for (; n != memo; n = PNODE(n->blossom_next)) {
    e = PEDGE(n->blossom_next_edge);
    e->x = 1 - e->x;
    if (e->x == 1) {
      nodelist[e->nod1].status = MATCHED;
      nodelist[e->nod2].status = MATCHED;
      nodelist[e->nod1].matched_edge = IEDGE(e);
      nodelist[e->nod2].matched_edge = IEDGE(e);
    }
  }

  oldnode->blossom_root = INODE(x);
  x->matched_edge = oldnode->matched_edge;
  x->status = MATCHED;
}

static void fix_tree (graph *G, node *oldnode, node *x, node *y)
{
  node *p, *n, *cnode;
  int c, label_help;

  p = PNODE(oldnode->parent);

#if PRINT_LEVEL
  printf ("    Fix Tree %i; x = %i; y = %i\n",
	  (int) (oldnode - PG->nodelist), (int) (x - PG->nodelist),
	  (int) (y - PG->nodelist));
  fflush (stdout);
#endif

  /* Restore child-sibling list for p */

  if (p->child == INODE(oldnode)) {
    y->sibling = oldnode->sibling;
    p->child = INODE(y);
  } else {
    for (c = p->child; c != -1; c = cnode->sibling) {
      cnode = PNODE(c);
      if (cnode->sibling == INODE(oldnode)) {
	cnode->sibling = INODE(y);
	y->sibling = oldnode->sibling;
	break;
      }
    }
    assert (c != -1);
  }

  y->parent = INODE(p);
  y->parent_edge = oldnode->parent_edge;
  G->nodelist[oldnode->child].parent = INODE(x);
  x->child = oldnode->child;

  label_help = MINUS;
  if (y->blossom_next_edge == y->matched_edge) {
    /* From y to x */
    for (n = y; n != x; n = PNODE(n->blossom_next)) {
      n->child = n->blossom_next;
      G->nodelist[n->blossom_next].parent = INODE(n);
      G->nodelist[n->blossom_next].parent_edge = n->blossom_next_edge;
      n->label = label_help;
      if (label_help == MINUS) {
	label_help = PLUS;
      } else {
	label_help = MINUS;
      }
    }
    x->label = MINUS;	/* This is not set in for-loop */
  } else {
    for (n = x; n != y; n = PNODE(n->blossom_next)) {
      /* From x to y */
      n->parent = n->blossom_next;
      n->parent_edge = n->blossom_next_edge;
      G->nodelist[n->blossom_next].child = INODE(n);
      n->label = label_help;
      if (label_help == MINUS) {
	label_help = PLUS;
      } else {
	label_help = MINUS;
      }
    }
    y->label = MINUS;	/* This is not set in for-loop */
  }
}

static int expand_blossom (graph *G, node *oldnode, stats *scount)
{
  node *memo, *x, *y, *n;
  int child, parent;

  scount->expand_count++;

#if PRINT_LEVEL
  printf ("    Expand blossom %i \n", (int) (oldnode - PG->nodelist));
  fflush (stdout);
#endif

  child = oldnode->child;
  parent = oldnode->parent;

  n = memo = PNODE(oldnode->blossom_root);
  do {
    label_penultimate (G, n, INODE(n));
    n = PNODE(n->blossom_next);
  } while (n != memo);

  if (G->edgelist[oldnode->matched_edge].nod1 == INODE(oldnode))
    x = PNODE(G->edgelist[oldnode->matched_edge].orig_nod1);
  else
    x = PNODE(G->edgelist[oldnode->matched_edge].orig_nod2);

  if (G->edgelist[oldnode->parent_edge].nod1 == INODE(oldnode))
    y = PNODE(G->edgelist[oldnode->parent_edge].orig_nod1);
  else
    y = PNODE(G->edgelist[oldnode->parent_edge].orig_nod2);

  if (lower_edges (G, oldnode)) {
    fprintf (stderr, "lower_edges failed\n");
    return 1;
  }
  fix_matching (G, oldnode, PNODE(x->penultimate));

  n = memo;
  do {
    /* Clear tree structure of nodes in oldnodes blossom */
    n->blossom_parent = -1;
    n->child = -1;
    n->parent = -1;
    n->sibling = -1;
    n->parent_edge = -1;
    n->label = NIX;

    n = PNODE(n->blossom_next);
    n->hit = n->hit | oldnode->hit;
  } while (n != memo);

  fix_tree (G, oldnode, PNODE(x->penultimate), PNODE(y->penultimate));

  /* update tree roots */

  for (child = G->nodelist[child].parent; child != parent;
       child = G->nodelist[child].parent) {
    G->nodelist[child].tree_root = G->nodelist[parent].tree_root;
  }

  /* clear oldnode */

  oldnode->edg_list = -1;
  oldnode->pi = 0;
  oldnode->status = UNMATCHED;
  oldnode->matched_edge = -1;
  oldnode->child = -1;
  oldnode->sibling = -1;
  oldnode->parent = -1;
  oldnode->parent_edge = -1;
  oldnode->label = NIX;
  oldnode->blossom_root = -1;
  oldnode->blossom_next = -1;
  oldnode->blossom_next_edge = -1;
  oldnode->blossom_parent = -1;
  oldnode->penultimate = -1;
  oldnode->tree_root = -1;
  oldnode->hit = 0;

  /* put oldnode on unused list */

  oldnode->mark = G->unused;
  G->unused = oldnode - G->nodelist;

  return 0;
}

static node *common_parent (graph *G, node *p, node *q, int *size)
{
  int count;
  node *n;
  node *stop = G->nodelist - 1;

#if PRINT_LEVEL
  printf ("     Common-parent from %i and %i is ",
          (int) (p - PG->nodelist), (int) (q - PG->nodelist));
  fflush (stdout);
#endif

  for (count = 0, n = p; n != stop; n = PNODE(n->parent)) {
    count++;
    n->mark = count;
  }

  for (count = 0; q != stop; q = PNODE(q->parent)) {
    if (q->mark) {
      *size = (count + q->mark);
      for (n = p; n != stop; n = PNODE(n->parent)) 
	n->mark = 0;
      return q;
    }
    count++;
  }
  for (n = p; n != stop; n = PNODE(n->parent))
    n->mark = 0;		

  *size = 0;
  return (node *) NULL;
}

static void flip_cycle (graph *G, node *node_i)
{
  int count = 0, ok = 1;
  node *node_j, *node_k;
  edge *e, *memo;

#if PRINT_LEVEL
  printf ("     Flip cycle with root %i: ", (int) (node_i - PG->nodelist));
  fflush (stdout);
#endif

  /* init start (node_j, node_k, e) */

  e = PEDGE(node_i->matched_edge);
  node_j = node_i;
  node_k = OTHEREND (e, node_i, G->nodelist);

  while (ok) {
    /* test if last edge */
    if (node_k == node_i)
      ok = 0;

    e->x = (count % 2);
    memo = PEDGE(node_k->matched_edge);

    if (count % 2) {
      node_j->matched_edge = node_k->matched_edge = IEDGE(e);
      node_j->status = MATCHED;
      node_k->status = MATCHED;
    }

    e = memo;
    node_j = node_k;
    node_k = OTHEREND (e, node_j, G->nodelist);
    count++;
  }
}

static void augment_path (graph *G, node *n, int fractional)
{
  edge *edgelist = G->edgelist;
  node *nodelist = G->nodelist;

#if PRINT_LEVEL
  printf ("      Augment path %i", (int) (n - PG->nodelist));
  fflush (stdout);
#endif

  if (n->parent != -1) {
    for (; n->parent_edge != -1; n = PNODE(n->parent)) {
      edgelist[n->parent_edge].x = 1 - edgelist[n->parent_edge].x;
      if (edgelist[n->parent_edge].x == 1) {
	nodelist[edgelist[n->parent_edge].nod1].matched_edge =
	  n->parent_edge;
	nodelist[edgelist[n->parent_edge].nod2].matched_edge =
	  n->parent_edge;
      }
    }
  }
  n->status = MATCHED;  /* This will be corrected in HALVES case */
  if (!fractional) 
    G->roots[n->tree_root].status = MATCHED;
  clear_tree (G, n);   /* n is the root */
}

static int augment_blossom (graph *G, edge *e, int fractional, srkstuff *srk)
{
  node *node_i = PNODE(e->nod1);
  node *node_j = PNODE(e->nod2);
  node *node_k;
  int size;
  shrink_ary_T *ary = srk->shrinks.ary;

#if PRINT_LEVEL
  printf ("    Augment blossom with ends %i & %i \n", e->nod1, e->nod2);
  fflush (stdout);
#endif

  node_k = common_parent (G, node_i, node_j, &size);
  if (node_k == (node *) NULL) {    /* more then one tree */
    e->x = 1;
    node_i->matched_edge = e -G->edgelist;
    node_j->matched_edge = e -G->edgelist;
    G->unmatched_count -= 2;
    augment_path (G, node_i, fractional);
    augment_path (G, node_j, fractional);
    return 1;
  } else {
    if (fractional) {
      make_cycle_halves (G, node_i, node_j, node_k, e);
      augment_path (G, node_k, fractional);
      node_k->status = HALVES;
      G->unmatched_count--;
      return 1;
    } else {
      if (srk->shrinks.length < srk->shrinks_max) {
	ary[srk->shrinks.length].node_i = node_i - G->nodelist;
	ary[srk->shrinks.length].node_j = node_j - G->nodelist;
	ary[srk->shrinks.length].node_k = node_k - G->nodelist;
	ary[srk->shrinks.length].edge = e - G->edgelist;
	ary[srk->shrinks.length].size = size;
	srk->shrinks.length++;
	if (srk->shrinks.length == srk->shrinks_max) {
	  fprintf (stderr,"   WARNING: shrinks.length==shrinks_max=%i\n",
		  srk->shrinks_max);
	  fflush (stdout);
	}
      }
      return 0;
    }
  }
}

static int add_to_tree (graph *G, edge *e, node *node_i, node *node_j,
                        int fractional)
{
  node *node_k;
  edge *f;

#if PRINT_LEVEL
  printf ("    Add (%i,%i)=(%i,%i) to tree \n", e->nod1,
	  e->nod2, (int) (node_i - PG->nodelist), (int) (node_j - PG->nodelist));
  fflush (stdout);
#endif

  if (node_j->status == UNMATCHED) {
    e->x = 1;
    node_j->status = MATCHED;
    node_i->matched_edge = e - G->edgelist;
    node_j->matched_edge = e - G->edgelist;
    G->unmatched_count -= 2;
    if (!fractional) 
      G->roots[node_j->tree_root].status = MATCHED;
    augment_path (G, node_i, fractional);
    return 1;
  } else if (node_j->status == HALVES) {
    flip_cycle (G, node_j);
    e->x = 1;
    node_i->status = MATCHED;
    node_j->status = MATCHED;
    node_i->matched_edge = e - G->edgelist;
    node_j->matched_edge = e - G->edgelist;
    augment_path (G, node_i, fractional);
    G->unmatched_count--;
    return 1;
  } else {
    /* set node_j in tree */
    node_j->sibling = node_i->child;
    node_j->parent = INODE(node_i);
    node_j->tree_root = node_i->tree_root;
    node_j->parent_edge = IEDGE(e);
    node_j->label = MINUS;

    /* set child for node_i */
    node_i->child = INODE(node_j);

    /* set child for node_j */
    f = PEDGE(node_j->matched_edge);
    node_k = OTHEREND (f, node_j, G->nodelist);
    node_j->child = INODE(node_k);

    /* set node_k in tree */
    node_k->child = -1;
    node_k->sibling = -1;
    node_k->parent = INODE(node_j);
    node_k->parent_edge = IEDGE(f);
    node_k->label = PLUS;
    node_k->tree_root = node_j->tree_root;

    return 0;
  }
}

static int checkout_node (graph *G, node *n, int fractional, srkstuff *srk)
{
  node *m;
  int augmented;
  edge *e;
  int ei;

#if PRINT_LEVEL
  printf ("   Checkout node %i ...\n", (int) (n - PG->nodelist));
  fflush (stdout);
#endif
 
  for (ei = n->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
    e = PEDGE(ei/2);
    if (e->slack == 0) {
      m = OTHEREND (e, n, G->nodelist);
      if (m->label == PLUS) {
	return augment_blossom (G, e, fractional, srk);
      }
      if (m->label == NIX) {
	if ((augmented = add_to_tree (G, e, n, m, fractional)) != 0) {
	  return augmented;
	}
      }
    }
  }

  return 0;
}

static int grow_tree_no_shrink (graph *G, node *n, int fractional,
				srkstuff *srk)
{
  int augmented;
  expand_T *expands = &(srk->expands);
  int expands_max = srk->expands_max;
  node *c=n;
  node *stop = G->nodelist - 1;

#if PRINT_LEVEL
  printf ("   growtree no shrink %i\n", (int) (c - PG->nodelist));
  fflush (stdout);
#endif

  while (1) {
    if (c->label == PLUS) {
      if ((augmented = checkout_node (G, c, fractional, srk)) > 0) {
	return augmented;
      }
    } else {			/* MINUS node */
      if (c->blossom_root != -1 && c->pi == 0) {
	if (expands->length < expands_max) {
	  expands->node[expands->length] = c - G->nodelist;
	  expands->length++;
	  if (expands->length == expands_max) {
	    printf ("   WARNING: expands = expands_max = %i\n",
		    expands_max);
	    fflush (stdout);
	  }
	}
      }
    }
    /* go to next c */
    if (c->child!=-1) {
      c=PNODE(c->child);
    } else {
      while (PNODE(c->sibling) == stop) {
	if (c == n) return 0; /* this if is for childless n */
	c = PNODE(c->parent);
	if (c == n) return 0;
      }
      c = PNODE(c->sibling);
    }
  }
}

static int grow_tree (graph *G, node *n, int fractional, stats *scount,
                      srkstuff *srk, int *found_aug)
{
  int i, size;
  node *node_i, *node_j, *node_k, *newnode;
  edge *e;
  shrink_ary_T *ary = srk->shrinks.ary;
  expand_T *expands = &(srk->expands);
  int *shrinks_ = srk->shrinks_;
  shrink_T *shrinks = &(srk->shrinks);
  int old_shrinks_len;

  *found_aug = 0;

  do {
    FIND_SURFACE (n);
    shrinks->length = 0;
    expands->length = 0;

    if (grow_tree_no_shrink (G, n, fractional, srk)) {
      *found_aug = 1;
      return 0;         
    }

    /* found no new match edge, so do shrinks */

    while (shrinks->length) {
      old_shrinks_len = shrinks->length;
      if (shrinks->length > 4) {   
	/* build buckets for all different sizes */
	for (i = 0; i < SHRINKS_SIZE; i++)
	  shrinks_[i] = -1;
	for (i = 0; i < shrinks->length; i++) {
	  size = ary[i].size;
	  if (size >= SHRINKS_SIZE)
	    size = SHRINKS_SIZE - 1;
	  if (shrinks_[size] == -1) {
	    shrinks_[size] = i;
	    ary[i].next = -1;
	  } else {
	    ary[i].next = shrinks_[size];
	    shrinks_[size] = i;
	  }
	}
	for (size = SHRINKS_SIZE - 1; size >= 3; size -= 2) {
	  for (i = shrinks_[size]; i != -1; i = ary[i].next) {
	    node_i = PNODE(ary[i].node_i);
	    node_j = PNODE(ary[i].node_j);
	    node_k = PNODE(ary[i].node_k);
	    e = PEDGE(ary[i].edge);
	    FIND_SURFACE (node_i);
	    FIND_SURFACE (node_j);
	    FIND_SURFACE (node_k);
	    if (node_i != node_j) {
	      newnode = shrink_blossom (G, node_i, node_j,
					node_k, e, scount);
	      if ((checkout_node(G,newnode,fractional,srk)) >0) {
		*found_aug = 1;
		return 0;
	      }
	    }
	  }
	}
      } else {      /* just do all of the stupid shrinks */
	for (i = 0; i < shrinks->length; i++) {
	  node_i = PNODE(ary[i].node_i);
	  node_j = PNODE(ary[i].node_j);
	  node_k = PNODE(ary[i].node_k);
	  e = PEDGE(ary[i].edge);
	  FIND_SURFACE (node_i);
	  FIND_SURFACE (node_j);
	  FIND_SURFACE (node_k);
	  if (node_i != node_j) {
	    newnode = shrink_blossom (G, node_i, node_j,
				      node_k, e, scount);
	    if ((checkout_node (G, newnode, fractional, srk)) > 0) {
	      *found_aug = 1;
	      return 0;
	    }
	  }
	}
      }
      if (shrinks->length > old_shrinks_len) {
	size = shrinks->length - old_shrinks_len;
	for (i = 0; i < size; i++) {
	  ary[i].node_i = ary[i + old_shrinks_len].node_i;
	  ary[i].node_j = ary[i + old_shrinks_len].node_j;
	  ary[i].node_k = ary[i + old_shrinks_len].node_k;
	  ary[i].edge = ary[i + old_shrinks_len].edge;
	  ary[i].size = ary[i + old_shrinks_len].size;
	}
	shrinks->length = size;
      } else {
	shrinks->length = 0; 
      }
    }
    for (i = 0; i < expands->length; i++) {
      /* blossom_root != -1 is necessary because expand
       * of this node might have happened here before */
      /* blossom_parent == -1 is necessary because
       * this node might be shrunk already */
      if (G->nodelist[expands->node[i]].blossom_root != -1 &&
	  G->nodelist[expands->node[i]].blossom_parent == -1) {
	if (expand_blossom (G, G->nodelist+expands->node[i], scount)) {
	  fprintf (stderr, "expand_blossom failed\n");
	  return 1;
	}
      }
    }
  } while (expands->length > 0);
  return 0;
}

static int apply_dual_change (graph *G, node *n, int delta)
{
  int ei;
  edge *e;
  node *c=n;
  node *stop = G->nodelist - 1;

  if (delta == INTEGER_MAX) {
    fprintf (stderr, "\ndelta=Infinity, node=%i\n", (int) (n-G->nodelist));
    fprintf (stderr, "There seems to be no perfect matching\n");
#ifndef NDEBUG
    print_node (G, n);
#endif
    return 1;
  }

  while (1) {
    if (c->label == PLUS) {
      c->hit = 1;
      c->pi += delta;
      for (ei = c->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	e = PEDGE(ei/2);
	e->slack -= delta;
      }
    } else if (c->label == MINUS) {
      c->pi -= delta;
      for (ei = c->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	e = PEDGE(ei/2);
	e->slack += delta;
      }
    }

    /* go to next c */

    if (c->child != -1) {
      c = PNODE(c->child);
    } else {
      while (PNODE(c->sibling) == stop) {
	if (c == n) return 0; /* this if is for childless n */
	c = PNODE(c->parent);
	if (c == n) return 0;
      }
      c = PNODE(c->sibling);
    }
  }
}

static int find_single_dual_change (graph *G, node *n)
{
  edge *e;
  int ei;
  int lab;
  int delta;
  node *c = n;
  node *stop = G->nodelist - 1;

#if PRINT_LEVEL
  printf (" %i", (int) (n - PG->nodelist));
  fflush (stdout);
#endif

  delta = INTEGER_MAX;

  while (1) {
    if (c->label == PLUS) {
      for (ei = c->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	e = PEDGE(ei/2);
	if (e->slack < 2 * delta) {
	  lab = OTHEREND (e, c, G->nodelist)->label;
	  if (lab == PLUS) {
	    delta = e->slack / 2;
	  } else if (lab == NIX) {
	    if (e->slack < delta) {
	      delta = e->slack;
	    }
	  }
	}
      }
    } else if (c->label == MINUS) {
      if (c->blossom_root != -1 && c->pi < delta) {
	delta = c->pi;
      }
    }

    /* go to next c */

    if (c->child != -1) {
      c = PNODE(c->child);
    } else {
      while (PNODE(c->sibling) == stop) {
	if (c == n) return delta ; /* this if is for childless n */
	c = PNODE(c->parent);
	if (c == n) return delta;
      }
      c = PNODE(c->sibling);
    }
  }
}

static int find_dual_change_forest (graph *G, node *n)
{
  edge *e;
  int ei;
  int delta;
  nodeptr *roots = G->roots;
  node *c = n, *o;
  node *stop = G->nodelist - 1;

#if PRINT_LEVEL
  printf (" %i", (int) (c - PG->nodelist));
  fflush (stdout);
#endif

  delta = INTEGER_MAX;

  while (1) {
    if (c->label == PLUS) {
      for (ei = c->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
	e = PEDGE(ei/2);
	o = OTHEREND (e, c, G->nodelist);
	switch (o->label) {
#ifdef GREEDY_DUAL_CHANGE
	case PLUS:
	  if (roots[c->tree_root].dad != roots[o->tree_root].dad ) {
	    /* two trees, not connected with vereinigung */
	    if (roots[o->tree_root].tree_processed) {
	      if (delta > e->slack -
		  roots[roots[o->tree_root].dad].delta) {
		delta = e->slack -
		  roots[roots[o->tree_root].dad].delta;
	      }
	    } else {
	      if (e->slack < delta) {
		delta = e->slack;
	      }
	    }
	  } else {
	    if (e->slack < 2 * delta) {
	      delta = e->slack / 2;
	    }
	  }
	  break;
	case NIX:
	  if (e->slack < delta) {
	    delta = e->slack;
	  }
	  break;
	case MINUS:
	  if (roots[c->tree_root].dad != roots[o->tree_root].dad) {
	    /* two trees, not connected with vereinigung */
	    if (roots[o->tree_root].tree_processed) { 
	      if (delta > e->slack +
		  roots[roots[o->tree_root].dad].delta) {
		delta = e->slack +
		  roots[roots[o->tree_root].dad].delta;
	      }
	    } else {
	      if (e->slack < delta ) {
		delta = e->slack;
	      }
	    }
	  }
	}
#else /* GREEDY_DUAL_CHANGE */
      case PLUS:
	if (e->slack < 2 * delta) {
	  delta = e->slack / 2;
	}
	break;
      case NIX:
	if (e->slack < delta) {
	  delta = e->slack;
	}
	break;
      case MINUS:
	if (roots[c->tree_root].dad != roots[o->tree_root].dad) { 
	  if (e->slack < delta ) {
	    delta = e->slack;
	  }
	}
      }
#endif /* GREEDY_DUAL_CHANGE */
    }
  } else if (c->label == MINUS) {
    if (c->blossom_root != -1 && c->pi <= delta) {
      delta = c->pi;
    }
  }
	
    /* go to next c */

    if (c->child != -1) {
      c = PNODE(c->child);
    } else {
      while (PNODE(c->sibling) == stop) {
	if (c == n) return delta ; /* this is for childless n */
	c = PNODE(c->parent);
	if (c == n) return delta;
      }
      c = PNODE(c->sibling);
    }
}
}

static void vereinigung (graph *G, int x, int y)
{
int xroot, yroot, dad;
nodeptr *roots = G->roots;

#if PRINT_LEVEL
printf (" vereinigung (%i,%i)", x, y); fflush (stdout);
#endif

xroot = x;
while (roots[xroot].dad >= 0)
     xroot = roots[xroot].dad;
     while (roots[x].dad >= 0) {
dad = roots[x].dad;
roots[x].dad = xroot;
x = dad;
}

yroot = y;
while (roots[yroot].dad >= 0)
     yroot = roots[yroot].dad;
     while (roots[y].dad >= 0) {
dad = roots[y].dad;
roots[y].dad = yroot;
y = dad;
}

if (xroot != yroot) {
if (roots[yroot].dad < roots[xroot].dad) {
roots[yroot].dad += roots[xroot].dad; 
roots[xroot].dad = yroot;
} else {
roots[xroot].dad += roots[yroot].dad; 
roots[yroot].dad = xroot;
}
}
}

static void set_vereinigung (graph *G, node *n)
{
  edge *e;
  int ei;
  node *nod = G->nodelist;
  node *c = n;
  node *stop = G->nodelist - 1;

#if PRINT_LEVEL
  printf (" set_vereinigung (%i)", (int) (n - PG->nodelist));
  fflush (stdout);
#endif

  /* this is called only for PLUS nodes */

  while (1) {
    for (ei = c->edg_list; ei != -1; ei = e->ptrs[ei % 2]) {
      e = PEDGE(ei/2);
      if (nod[e->nod1].tree_root != nod[e->nod2].tree_root) {  
	if (e->slack == 0) {	
	  if (OTHEREND (e, c, nod)->label == MINUS) {
	    vereinigung (G, nod[e->nod1].tree_root,
			 nod[e->nod2].tree_root);
	  }
	}
      }
    }
	
    /* now check the chilren of c's children (these are PLUS nodes) */

    if (c->child != -1) {
      c = PNODE(PNODE(c->child)->child);
    } else {
      if (c == n) return; /* this if is for childless n */
      c = PNODE(c->parent);
      while (PNODE(c->sibling) == stop) {
	c = PNODE(c->parent);
	if (c == n) return;
	c = PNODE(c->parent);
      }
      c = PNODE(PNODE(c->sibling)->child);
    }
  }
}

static int find_parity_sum (graph *G, int n)
{
  node *v;
  int sum, ei;
  edge *e;

  ei = PNODE(n)->edg_list;
  e = PEDGE(ei/2);

  if (e->nod1 == n) {
    v = PNODE(e->orig_nod1); 
  } else {
    v = PNODE(e->orig_nod2);
  }

  sum = v->pi;
  while (v->blossom_parent != -1) {
    v = PNODE(v->blossom_parent);
    sum += v->pi;
  }

  return sum;
}

static int parity_correction (graph *G, stats *scount)
{
  nodeptr *np, *npprev = (nodeptr *) NULL;
  int hitme = 0;
  nodeptr *roots = G->roots;
  nodeptr *rstop = G->roots - 1;

#if PRINT_LEVEL
  printf ("   parity_correction ...\n"); 
  fflush (stdout);
#endif

  np = roots + G->unmatched;
  while (np->status != UNMATCHED)  
    np = roots + np->next;
  G->unmatched = np - roots;

  while (np != rstop) {
    if (np->status != UNMATCHED) {
      npprev->next = np->next;
    } else {
      if (np->sum % 2) {
	hitme = 1;
	scount->dualchange_count++;
	if (apply_dual_change (G, PNODE(np->surf), 1)) {
	  fprintf (stderr, "apply_dual_change failed\n");
	  return -1;
	}
	np->sum += 1; 
      }
      npprev = np;
    }
    np = roots + np->next;
  }
  return hitme;
}

static int make_dual_change_forest (graph *G, stats *scount)
{
  nodeptr *np;
  int t;
  int delta;
  nodeptr *roots = G->roots;
  nodeptr *rstop = G->roots - 1;

#if PRINT_LEVEL
  printf ("\n   make_dual_change_forest ...\n");
  fflush (stdout);
#endif

  if (parity_correction (G, scount) == -1) {
    fprintf (stderr, "parity_correction failed\n");
    return -1;
  }

  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) 
    np->dad = -1;     /*  Set dad for all trees to -1 */

  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) 
    set_vereinigung (G, PNODE(np->surf));

  /* Set roots to point to theirselves */
  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
    if (np->dad < 0) {
      np->dad = np - roots;
      np->delta = INTEGER_MAX; 
    }
  }

  /* Set others to point to root */
  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
    t = np - roots;
    while (roots[t].dad != t)
      t = roots[t].dad;
    np->dad = t;
  }

#ifdef GREEDY_DUAL_CHANGE
  {
    int npnext;
    int glist = -1;

    /* order the unmatched nodes so that each union appears consecutively */

    for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
      np->tree_processed = 0;
      if (np->dad == np - roots) {
	np->gnext = glist;
	glist = np - roots;
      }
    }

    for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
      if (np->dad != np - roots) {
	np->gnext = roots[np->dad].gnext;
	roots[np->dad].gnext = np - roots;
      }
    }

    for (np = roots + G->unmatched; np != rstop; np = roots + npnext) {
      npnext = np->next;
      np->next = np->gnext;
    }
    G->unmatched = glist;
  }
#endif

  /* Set delta for all vereinigung trees */
  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
    delta = find_dual_change_forest (G, PNODE(np->surf));
    if (delta < roots[np->dad].delta)
      roots[np->dad].delta = delta;
#ifdef GREEDY_DUAL_CHANGE
    np->tree_processed = 1;
#endif
  } 

  /* Apply_dual_change for all vereinigung trees */
  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
    delta = roots[np->dad].delta;
    scount->dualchange_count++;
    if (delta == 0) {
      scount->dualzero_count++;
    } else {
      if (apply_dual_change (G, PNODE(np->surf), delta)) {
	fprintf (stderr, "apply_dual_change failed\n");
	return -1;
      }
    }
    np->sum += delta; 
  }

  return 0;
}

static int match_main_frac (graph *G, stats *scount, srkstuff *srk)
{
  node *n;
  int i, delta;
  int found_aug = 0;

  /* Just work with one tree   */

  for (i = 0, n = G->nodelist; i < G->nnodes; i++, n++) {
    if (n->status == UNMATCHED) {
      init_tree (G, n);
      if (grow_tree (G, n, 1, scount, srk, &found_aug)) {
	fprintf (stderr, "grow_tree failed\n");
	return 1;
      }
      while (!found_aug) {
	scount->dualchange_count++;
	delta = find_single_dual_change (G, n);
	if (delta != 0) {
	  if (apply_dual_change (G, n, delta)) {
	    fprintf (stderr, "apply_dual_change failed\n");
	    return 1;
	  }
	} else {
	  scount->dualzero_count++;
	}
	if (grow_tree (G, n, 1, scount, srk, &found_aug)) {
	  fprintf (stderr, "grow_tree failed\n");
	  return 1;
	}
      }
#if PRINT_LEVEL
      printf ("."); fflush (stdout);
#endif
    }
  }
#ifndef NDEBUG
  printf (" %i Dual Changes, %i with delta=0 ",
	  scount->dualchange_count, scount->dualzero_count);
  fflush (stdout);
#endif
  return 0;
}

static int match_main_more_in_forest (graph *G, stats *scount, srkstuff *srk)
{
  nodeptr *np;
  node *n;
  int i, ni, ninext;
  int grow_status_forest;
  int found_aug = 0;
  nodeptr *roots, *rstop;

  /* The unmatched surface nodes should be linked using the mark  */
  /* fields, with G->unmatched pointing to the start of the list. */
  /* G->umatched is assumed to be accurate.                       */

  if (G->unmatched_count == 0) 
    return 0;

  /* ********************************* */
  /* First init all forests            */
  /* ********************************* */

  G->roots = CC_SAFE_MALLOC (G->unmatched_count, nodeptr);
  if (!G->roots) {
    fprintf (stderr, "out of memory in match_main_more_in_forest\n");
    return 1;
  }
  roots = G->roots;
  rstop = G->roots - 1;

  for (i = 0, ni = G->unmatched; ni != -1; ni = ninext, i++) {
    n = PNODE(ni);
    ninext = n->mark;
    n->mark = 0;
    roots[i].status = UNMATCHED;
    roots[i].surf = ni;
    roots[i].sum = find_parity_sum (G, ni);
    roots[i].next = i+1;
    init_tree (G, n);
    n->tree_root = i;
  }
  roots[i - 1].next = -1;
  G->unmatched = 0;

#ifndef NDEBUG
  printf ("\nTry to grow the trees in the forest directly\n");
  fflush (stdout);
#endif

  /* ************************************ */
  /* Try to grow all forests directly     */
  /* ************************************ */

  for (np = roots + G->unmatched; np != rstop; np = roots + np->next) {
    n = PNODE(np->surf);
    if (np->status == UNMATCHED) {
      if (grow_tree (G, n, 0, scount, srk, &found_aug)) {
	fprintf (stderr, "grow_tree failed\n");
	return 1;
      }
      if (found_aug) {	 
#if PRINT_LEVEL
      printf ("."); fflush (stdout);
#endif
	if (G->unmatched_count == 0) {
	  return 0;
	}
      }
    }
  }

  /* *************************** */
  /* Work with tree-vereinigung  */
  /* *************************** */

#ifndef NDEBUG
  printf ("\nNow work on tree-vereinigung in a forest (%i points)\n",
	  G->unmatched_count / 2);
  fflush (stdout);
#endif

  /* first make sure that all trees have the correct parity */

  if (parity_correction (G, scount) == -1) {
    fprintf (stderr, "parity_correction failed\n");
    return 1;
  }

  while (G->unmatched_count > 0) {
    if (make_dual_change_forest (G, scount) == -1) {
      fprintf (stderr, "make_dual_change_forest failed\n");
      return 1;
    }
    do {
      grow_status_forest = 0;
      for (np = roots + G->unmatched; np != rstop; np = roots+np->next) {
	n = PNODE(np->surf);
	if (np->status == UNMATCHED) {
	  if (grow_tree (G, n, 0, scount, srk, &found_aug)) {
	    fprintf (stderr, "grow_tree failed\n");
	    return 1;
	  }
	  if (found_aug) { 
#if PRINT_LEVEL
	    printf ("."); fflush (stdout);
#endif
	    if (G->unmatched_count == 0) {
	      goto DONE;
	    }
	    grow_status_forest = 1;
	  }
	}
      }
    } while (grow_status_forest);
  }

 DONE:

  G->unmatched = -1;
  CC_IFFREE (G->roots, nodeptr);

#ifndef NDEBUG
  printf ("\n %i Dual Changes, %i with delta=0 ",
	  scount->dualchange_count, scount->dualzero_count);
  printf ("| %i Expands ", scount->expand_count);
  printf ("| %i Shrinks ", scount->shrink_count);
  fflush (stdout);
#endif

  return 0;
}

static int match_main_forest (graph *G, stats *scount, srkstuff *srk)
{
  nodeptr *np;
  node *n;
  int i, ni, ninext;
  int grow_status_forest;
  int found_aug = 0;
  int delta, delta2;
  int do_parity = 1;
  nodeptr *roots, *rstop;

#ifndef NDEBUG
  printf ("match_main_forest (%d points)\n", G->unmatched_count / 2);
  fflush (stdout);
#endif

  if (G->unmatched_count == 0) {
    return 0;
  }

  /* ********************************* */
  /* First init all forests            */
  /* ********************************* */

  G->roots = CC_SAFE_MALLOC (G->unmatched_count, nodeptr);
  if (!G->roots) {
    fprintf (stderr, "out of memory in match_main_more_in_forest\n");
    return 1;
  }
  roots = G->roots;
  rstop = G->roots - 1;

  for (i = 0, ni = G->unmatched; ni != -1; ni = ninext, i++) {
    n = PNODE(ni);
    ninext = n->mark;
    n->mark = 0;
    roots[i].status = UNMATCHED;
    roots[i].surf = ni;
    roots[i].sum = find_parity_sum (G, ni);
    roots[i].next = i+1;
    init_tree (G, n);
    n->tree_root = i;
  }
  roots[i - 1].next = -1;
  G->unmatched = 0;

  /* ********************************* */
  /* Work with a forest                */
  /* ********************************* */

  while (G->unmatched_count > 0) {
    do {
      grow_status_forest = 0;
      for (np = roots + G->unmatched; np != rstop; np = roots+np->next) {
	n = PNODE(np->surf);
	if (np->status == UNMATCHED) {
	  if (grow_tree (G, n, 0, scount, srk, &found_aug)) {
	    fprintf (stderr, "grow_tree failed\n");
	    return 1;
	  }
	  if (found_aug) {
#if PRINT_LEVEL
	    printf ("."); fflush (stdout);
#endif
	    if (G->unmatched_count == 0) {
	      goto DONE;
	    }
	    grow_status_forest = 1;
	  }
	}
      }
    } while (grow_status_forest);

    if (G->unmatched_count > 0 && do_parity) {
      if (parity_correction (G, scount) == -1) {
	fprintf (stderr, "parity_correction failed\n");
	return 1;
      }
      do_parity = 0;
    }
    
    if (G->unmatched_count > 0 && !do_parity) {
      delta = INTEGER_MAX;
      scount->dualchange_count++;
      for (np = roots + G->unmatched; np != rstop; np = roots+np->next) {
	n = PNODE(np->surf);
	if (np->status == UNMATCHED) {
	  delta2 = find_single_dual_change (G, n);
	  if (delta2 < delta) {
	    delta = delta2;
	  }
	}
      } 

      if (delta == 0) {
	scount->dualzero_count++;
      } else {
	for (np = roots+G->unmatched; np!=rstop; np = roots+np->next) {
	  n = PNODE(np->surf);
	  if (np->status == UNMATCHED) {
	    if (apply_dual_change (G, n, delta)) {
	      fprintf (stderr, "apply_dual_change failed\n");
	      return 1;
	    }
	    np->sum += delta; 
	  }
	}
      }
    }
  }

 DONE:

  G->unmatched = -1;
  return 0;
}

static int match_main_tree (graph *G, stats *scount, srkstuff *srk)
{
  nodeptr *np;
  node *n;
  int i, ni, ninext;
  int found_aug = 0;
  int delta;
  nodeptr *roots, *rstop;

#ifndef NDEBUG
  printf ("match_main_tree (%d points)\n", G->unmatched_count / 2);
  fflush (stdout);
#endif

  if (G->unmatched_count == 0) {
    return 0;
  }

  /* ********************************* */
  /* First init all forests            */
  /* ********************************* */

  G->roots = CC_SAFE_MALLOC (G->unmatched_count, nodeptr);
  if (!G->roots) {
    fprintf (stderr, "out of memory in match_main_more_in_forest\n");
    return 1;
  }
  roots = G->roots;
  rstop = G->roots - 1;

  for (i = 0, ni = G->unmatched; ni != -1; ni = ninext, i++) {
    n = PNODE(ni);
    ninext = n->mark;
    n->mark = 0;
    roots[i].status = UNMATCHED;
    roots[i].surf = ni;
    roots[i].sum = find_parity_sum (G, ni);
    roots[i].next = i+1;
    /* init_tree (G, n); */
    n->tree_root = i;
  }
  roots[i - 1].next = -1;
  G->unmatched = 0;

  /* ********************************* */
  /* Work with a tree                  */
  /* ********************************* */

  while (G->unmatched_count > 0) {
    np = roots + G->unmatched; 
    while (np->status != UNMATCHED) {
      np = roots + np->next;
    }
    G->unmatched = np - roots;
    np = roots + G->unmatched; 
    n = PNODE(np->surf);
    init_tree (G, n);
    do {
      if (grow_tree (G, n, 0, scount, srk, &found_aug)) {
	fprintf (stderr, "grow_tree failed\n");
	return 1;
      }
      if (found_aug) {
#if PRINT_LEVEL
	printf ("."); fflush (stdout);
#endif
	if (G->unmatched_count == 0) {
	  goto DONE;
	}
      } else {
	scount->dualchange_count++;
	n = PNODE(np->surf);
	delta = find_single_dual_change (G, n);

	if (delta == 0) {
	  scount->dualzero_count++;
	} else {
	  if (apply_dual_change (G, n, delta)) {
	    fprintf (stderr, "apply_dual_change failed\n");
	    return 1;
	  }
	}
      }
    } while (!found_aug);
  }

 DONE:

  G->unmatched = -1;
  return 0;
}

static int match_main (graph *G, stats *scount, srkstuff *srk, int use_all)
{
#ifndef NDEBUG
  printf (" Need to find %i edges ...\n", G->unmatched_count / 2);
  fflush (stdout);
#endif

  if (!G->unmatched_count) 
    return 0;

  if (use_all == 2) {
    if (match_main_tree (G, scount, srk)) {
      fprintf (stderr, "match_main_tree failed\n");
      return 1;
    }
  } else if (use_all == 1) {
    if (match_main_forest (G, scount, srk)) {
      fprintf (stderr, "match_main_forest failed\n");
      return 1;
    }
  } else {
    if (match_main_more_in_forest (G, scount, srk)) {
      fprintf (stderr, "match_main_more_in_forest failed\n");
      return 1;
    }
  }
#ifndef NDEBUG
  printf ("\n %i Dual Changes, %i with delta=0 \n",
	  scount->dualchange_count, scount->dualzero_count);
  printf (" %i Expands \n", scount->expand_count);
  printf (" %i Shrinks \n", scount->shrink_count);
#endif

  return 0;
}

static node *find_below (graph *G, node *n, int blossom)
{
  if (!n)
    return (node *) NULL;

  for (; n->blossom_parent != blossom; n = PNODE(n->blossom_parent)) {
    if (n->blossom_parent == -1)
      return (node *) NULL;
  }
  return n;
}

static void fix_match (graph *G, int blossom)
{
  node *x, *n;
  node *b = blossom + G->nodelist;
  
#if PRINT_LEVEL
  printf ("   Fixing blossom %i ", blossom);
  fflush (stdout);
#endif

  x = find_below (G, PNODE(G->edgelist[b->matched_edge].orig_nod1),
		  blossom);
  if (x == (node *) NULL) {
    x = find_below (G, PNODE(G->edgelist[b->matched_edge].orig_nod2),
		    blossom);
  }
  /* x is now something like the match_edge->nod(1 or 2)->penultimate */
  fix_matching (G, b, x);
  b->mark = 1;

  n = PNODE(b->blossom_root);
  do {
    if (n->blossom_root != -1)          /* if n is a real blossom */
      fix_match (G, INODE(n));        /* call fix match again */
    n = PNODE(n->blossom_next);
  } while (n != PNODE(b->blossom_root));
}

static void adjust_match (graph *G)
{
  node *n, *b;
  int i;
  int ncount = G->nnodes;
  node *nodelist = G->nodelist;

  for (i = 0; i < ncount; i++) {
    n = PNODE(i);
    if (n->blossom_parent != -1 && nodelist[n->blossom_parent].mark == 0) {
      for (b = n; b->blossom_parent != -1 &&
	     nodelist[b->blossom_parent].mark == 0;
	   b = PNODE(b->blossom_parent));
      fix_match (G, INODE(b));
    }
  }
}

static void dual_match_len (graph *G, int fractional, double *val) 
{
  int i;
  double b = 0.0;
  int ncount = G->nnodes;
  node *nodelist = G->nodelist;

  if (fractional) {
    for (i = 0; i < ncount; i++) {
      b += (double) nodelist[i].pi;
    }
  } else {
    for (i = 0; i <= 3*ncount/2; i++)
      b += (double) nodelist[i].pi;
  }
  *val = b / 2.0;
}

static int make_match (graph *G) 
{
  int i, i2, j, k, start, counter=0;
  int ncount = G->nnodes;
  edge *e;
  node *nodelist = G->nodelist;
  edge *edgelist = G->edgelist;

  for (i = 0; i < ncount; i++) {
    if (nodelist[i].status == MATCHED)
      counter++;
    if (nodelist[i].status == HALVES) {
      k = 1;
      start = i;
      i2 = i;
      j = OTHEREND_INT (PEDGE(nodelist[i2].matched_edge), i2);
      do {
	e = PEDGE(nodelist[j].matched_edge);
	if (k % 2) {
	  nodelist[i2].status = MATCHED;
	  nodelist[j].status = MATCHED;
	  edgelist[nodelist[j].matched_edge].x = 0;
	  nodelist[j].matched_edge = nodelist[i2].matched_edge;
	  edgelist[nodelist[j].matched_edge].x = 1;
	}
	i2 = j;
	j = OTHEREND_INT (e, i2);
	k++;
      } while (j!=start);
      nodelist[i2].status = UNMATCHED;
      edgelist[nodelist[i2].matched_edge].x = 0;
      nodelist[i2].matched_edge = -1;
      counter++;
    }
    nodelist[i].label = NIX;
    nodelist[i].child = -1;
    nodelist[i].sibling = -1;
    nodelist[i].parent = -1;
    nodelist[i].parent_edge = -1;
    nodelist[i].mark = 0;
    nodelist[i].tree_root = -1;
  }

  counter=0;
  for (i = 0; i < ncount; i++)
    if (nodelist[i].status == MATCHED)
      counter++;
#ifndef NDEBUG
  printf(" %i nodes are matched now !\n",counter);
#endif

  /* build list of unmatched nodes */

  G->unmatched_count = 0;
  for (i = 0; i < ncount; i++) {
    if (nodelist[i].status == UNMATCHED) {
      G->unmatched_count++;
    }
  }

  j = -1;
  for (i = 0; i < ncount; i++) {
    if (nodelist[i].status == UNMATCHED) {
      if (j != -1) 
	nodelist[j].mark = i;
      j = i;
    }
  }
  if (j != -1)
    nodelist[j].mark = -1;


  if (G->unmatched_count == 0) {
    G->unmatched = -1;
  } else {
    for (i = 0; nodelist[i].status != UNMATCHED; i++);
    G->unmatched = i;
  }
#ifndef NDEBUG
  printf(" unmatched_count = %i\n", G->unmatched_count);
  fflush (stdout);
#endif

  /* hit is used to identify the blossoms that receive a label +  during */
  /* the matching run; this will help the pricing code                   */

  for (i = 0; i < 3 * ncount / 2; i++)  
    nodelist[i].hit = 0;

  return 0;
}

static int build_graph (graph *G, int ncount, int ecount, float **dists)
{
  int i, j;
  edge *e;

  if (ncount % 2 != 0) {
    fprintf (stderr, "problem has an odd number of nodes\n");
    return 1;
  }

  G->nnodes = ncount;
  G->nodelist = CC_SAFE_MALLOC ((3*ncount/2+1), node);
  if (!G->nodelist) {
    fprintf (stderr, "out of memory in build_graph\n");
    return 1;
  }
  for (i = 0; i < 3*ncount/2 + 1; i++) {
    G->nodelist[i].blossom_parent = -1;
    G->nodelist[i].matched_edge = -1;
    G->nodelist[i].edg_list = -1;
  }

  /* init unused */

  j = 3*ncount/2;
  for (i = ncount; i < j; i++) {
    G->nodelist[i].mark = i + 1;
  }
  G->nodelist[j].mark = -1;
  G->unused = ncount;

  G->nedges = ecount;
  G->max_nedges = ecount + (1.5 * ncount);

  G->edgelist = CC_SAFE_MALLOC (G->max_nedges, edge);
  if (!G->edgelist) {
    fprintf (stderr, "out of memory in build_graph\n");
    CC_FREE (G->nodelist, node);
    return 1;
  }

  for (i = 0; i < ncount; i++) {
    G->nodelist[i].mark = 0;
  }

  unsigned int edge_num = 0;

  for(i = 0, e = G->edgelist; i < ncount; ++i) {
    for(j = 0; j < i; ++j) {
      if(dists[i][j] == MAXFLOAT) continue;
      e->nod1 = i;
      e->nod2 = j;
      e->orig_nod1 = e->nod1;
      e->orig_nod2 = e->nod2;
      e->slack = 2*dists[i][j];  /* to work with ints */
      e->mark = 0;
      e->ptrs[0] = G->nodelist[e->nod1].edg_list;
      G->nodelist[e->nod1].edg_list = 2*edge_num;
      e->ptrs[1] = G->nodelist[e->nod2].edg_list;
      G->nodelist[e->nod2].edg_list = 2*edge_num+1;
      ++e;
      ++edge_num;
    }
  }

  /*
    for (i = 0, e = G->edgelist; i < ecount; i++, e++) {
    e->nod1 = elist[2*i];
    e->nod2 = elist[2*i+1];
    e->orig_nod1 = e->nod1;
    e->orig_nod2 = e->nod2;
    e->slack =  2 * elen[i];
    e->mark = 0;
    e->ptrs[0] = G->nodelist[e->nod1].edg_list;
    G->nodelist[e->nod1].edg_list = 2*i;
    e->ptrs[1] = G->nodelist[e->nod2].edg_list;
    G->nodelist[e->nod2].edg_list = 2*i + 1;
    }
  */

  for (; i < G->max_nedges; i++)
    G->edgelist[i].mark = 0;

  return 0;
}

static int write_match (graph *G, float **dists, char* fname) 
{
  int i, c;
  float len;
  edge *e;
  FILE* fp;

  printf(" Writing %s ...", fname);
  fflush(stdout);

  if ((fp = fopen (fname, "w")) ==  (FILE *) NULL) {
    fprintf(stderr," Error: Can't open file %s\n",fname);
    return 1;
  }

  c = 0;
  for (i = 0; i < G->nnodes; i++) {
    if (i == G->edgelist[G->nodelist[i].matched_edge].orig_nod2 ||
	G->edgelist[G->nodelist[i].matched_edge].x == HALF) {
      c++;
    }
  }
  printf(" %i nodes, %i edges ", G->nnodes, c);
  fprintf (fp, "%i %i\n", G->nnodes, c);
  for (i = 0; i < G->nnodes; i++) {
    e = PEDGE(G->nodelist[i].matched_edge);
    if (i == e->orig_nod2 || e->x == HALF) {
      len = dists[G->edgelist[G->nodelist[i].matched_edge].orig_nod1][G->edgelist[G->nodelist[i].matched_edge].orig_nod2];
      fprintf (fp,"%i %i %d\n",
	       G->edgelist[G->nodelist[i].matched_edge].orig_nod1,
	       G->edgelist[G->nodelist[i].matched_edge].orig_nod2, (int)len);
    }
  }
  fclose(fp);
  printf("\n");

  return 0;
}

static void mark_match (graph *G, float **dists) 
{
  int i;
  float len;
  edge *e;

  for (i = 0; i < G->nnodes; i++) {
    e = PEDGE(G->nodelist[i].matched_edge);
    if (i == e->orig_nod2 || e->x == HALF) {
      len =  dists[G->edgelist[G->nodelist[i].matched_edge].orig_nod1][G->edgelist[G->nodelist[i].matched_edge].orig_nod2];
      dists[G->edgelist[G->nodelist[i].matched_edge].orig_nod1][G->edgelist[G->nodelist[i].matched_edge].orig_nod2] =
	(len) ? -len : -1;
    }
  }
}
