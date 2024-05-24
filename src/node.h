/////////////////////////////////////////////////////////////////////
// All rights reserved
/////////////////////////////////////////////////////////////////////

#ifndef MAZEH
#define MAZEH

#include "include.h"
#define MAXOBJ 5

#define MAXGCL 80000
#define MAXGOP 100


struct nodeSmall;
typedef struct nodeSmall nodeSmall;

struct nodeOp;
typedef struct nodeOp nodeOp;

struct gnode;
typedef struct gnode gnode;

struct snode;
typedef struct snode snode;

struct nodeCl;
typedef struct nodeCl nodeCl;

struct nodeCl
{
  int g[MAXGCL][MAXOBJ];
  int cont;
};

struct nodeOp
{
  snode *Gop[MAXGOP];
  int front;
  int rear;
};

struct nodeSmall
{
  nodeSmall *searchtree;
  int state;
};


struct gnode // stores info needed for each graph node
{
  long long int id;
  unsigned h[MAXOBJ];
  unsigned h1;
  unsigned h2;
  unsigned long long int key;
  unsigned gmin;
  unsigned long heapindex;
  nodeOp *Op;
  nodeCl *Cl;
};


struct snode // MOA*'s search nodes
{
  int state;
  unsigned g[MAXOBJ];
  unsigned g1;
  unsigned g2;
  double key;
  unsigned long heapindex;
  snode *searchtree;
};


#endif
