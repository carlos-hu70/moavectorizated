/////////////////////////////////////////////////////////////////////
// Carlos Hernandez
// All rights reserved
/////////////////////////////////////////////////////////////////////

#include "heap.h"
#include "node.h"
#include "include.h"
#include "moastarv.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include "immintrin.h"
#include <limits.h>

gnode* graph_node;
unsigned num_gnodes;
unsigned adjacent_table[MAXNODES][MAXNEIGH];
unsigned pred_adjacent_table[MAXNODES][MAXNEIGH];
unsigned goal, start, nobjs;
char filename[128];
gnode* start_state;
gnode* goal_state;
snode* start_node;

unsigned long long int stat_expansions = 0;
unsigned long long int stat_generated = 0;
unsigned long long int stat_gopoperations = 0;
unsigned long long int stat_gcloperations = 0;
unsigned long long int stat_gclupdate = 0;
unsigned long long int minf_solution = LARGE;
unsigned long long int stat_check1 = 0;

unsigned solutions[MAX_SOLUTIONS][2];
snode* recycled_nodes[MAX_RECYCLE];
int next_recycled = 0;

unsigned nsolutions = 0;
unsigned stat_pruned = 0;
unsigned stat_created = 0;

float runtimeUp = 0,runtimeCh = 0;
struct timeval tstart1, tend1;

void initialize_parameters() {
    start_state = &graph_node[start];
    goal_state = &graph_node[goal];
    stat_percolations = 0;
}

int backward_dijkstra(int dim) {
    for (int i = 0; i < num_gnodes; ++i)
        graph_node[i].key = LARGE;
    emptyheap_dij();
    goal_state->key = 0;
    insertheap_dij(goal_state);

    while (topheap_dij() != NULL) {
        gnode* n;
        gnode* pred;
        short d;
        n = popheap_dij();
        n->h[dim-1] = n->key;
        ++stat_expansions;
        for (d = 1; d < pred_adjacent_table[n->id][0] * (nobjs+1); d += nobjs+1) {
            pred = &graph_node[pred_adjacent_table[n->id][d]];
            int new_weight = n->key + pred_adjacent_table[n->id][d + dim];
            if (pred->key > new_weight) {
                pred->key = new_weight;
                insertheap_dij(pred);
            }
        }
    }
    return 1;
}

snode* new_node() {
	snode* state;
	
	if (next_recycled > 0) { //to reuse pruned nodes in memory
		state = recycled_nodes[--next_recycled];
     }
     else{
		state = (snode*)malloc(sizeof(snode));
        ++stat_created;
    }
    state->heapindex = 0;
		return state;
}

int isDominated(unsigned *v1, unsigned *v2){
	int i;
	for (i = 1; i < nobjs;i++){
		if(v1[i] < v2[i])
			return 0;
	}
	return 1;
}



void printfCl(gnode *n){
	if (n->Cl == NULL)
		return;
	int i;
	int l = n->Cl->cont;
	printf("[%d]\n",l);
	for (i = 0; i < l; i++){
		int k;
		printf("%d ",i);
		for (k = 0; k < nobjs;k++)
			printf("(%d)",n->Cl->g[i][k]);
		printf("\n");
	}
	
}

//*************** Begin Vectorized Functions ***********************************//
#define NPV 5000 //num_packed_vectors

__m512i goalv1[NPV];
__m512i goalv2[NPV];
__m512i goalv3[NPV];

int goalnp = 0;

void print_m512i(__m512i v) {
    int32_t result[16];
    _mm512_store_si512(result, v);

    printf("Vector: [");
    for (int i = 0; i < 16; ++i) {
        printf("%d", result[i]);
        if (i < 15) {
            printf(", ");
        }
    }
    printf("]\n");
}

void copyVector(__m512i* destination, const __m512i* source) {
    // Copy the contents from source to destination
    _mm512_store_si512(destination, _mm512_load_si512(source));
}


int PacketAddAll_Check_Dominance3(gnode *n, unsigned *v) 
{ 
	int i;
	int num_elements = n->Cl->cont;
	int num_packed_vectors = num_elements / 16;
    int num_unpacked_vectors = num_elements - num_packed_vectors * 16;
	if(num_elements > 15){
	__m512i v1[num_packed_vectors];
	__m512i v2[num_packed_vectors];

	__m512i x1 = _mm512_set1_epi32(v[1]); //Broadcast v[1] to all elements of v1
	__m512i x2 = _mm512_set1_epi32(v[2]); //Broadcast v[2] to all elements of v2
	
	for (i =0; i < num_packed_vectors; i++) {
		int base = i * 16;		
        //init dim 1
		v1[i] = _mm512_set_epi32(
		n->Cl->g[base + 0][1],
		n->Cl->g[base + 1][1],
		n->Cl->g[base + 2][1],
		n->Cl->g[base + 3][1],
		n->Cl->g[base + 4][1],
		n->Cl->g[base + 5][1],
		n->Cl->g[base + 6][1],
		n->Cl->g[base + 7][1],
		n->Cl->g[base + 8][1],
		n->Cl->g[base + 9][1],
		n->Cl->g[base + 10][1],
		n->Cl->g[base + 11][1],
		n->Cl->g[base + 12][1],
		n->Cl->g[base + 13][1],
		n->Cl->g[base + 14][1],
		n->Cl->g[base + 15][1]);
        //init dim 2
		v2[i] = _mm512_set_epi32(
		n->Cl->g[base + 0][2],
		n->Cl->g[base + 1][2],
		n->Cl->g[base + 2][2],
		n->Cl->g[base + 3][2],
		n->Cl->g[base + 4][2],
		n->Cl->g[base + 5][2],
		n->Cl->g[base + 6][2],
		n->Cl->g[base + 7][2],
		n->Cl->g[base + 8][2],
		n->Cl->g[base + 9][2],
		n->Cl->g[base + 10][2],
		n->Cl->g[base + 11][2],
		n->Cl->g[base + 12][2],
		n->Cl->g[base + 13][2],
		n->Cl->g[base + 14][2],
		n->Cl->g[base + 15][2]);
        __mmask32 mask1 = _mm512_cmp_epi32_mask(v1[i], x1, _MM_CMPINT_LE);
        __mmask32 mask2 = _mm512_cmp_epi32_mask(v2[i], x2, _MM_CMPINT_LE);

		//now perform *and* over results 
		__mmask32 and_res = mask1 & mask2;
		unsigned int result = _cvtmask32_u32(and_res);
		if (result != 0){    //test if and_res contains a positive vale
			return 1;
		}

	}
	}
	for (i = num_packed_vectors * 16; i< n->Cl->cont;i++){
		stat_gcloperations++;		
		if (isDominated(v, n->Cl->g[i])){
			return 1;
		}
	}
	return 0;
}

int PacketAddAll_Check_Dominance4(gnode *n, unsigned *v) 
{ 
	int i;
	int num_elements = n->Cl->cont;
	int num_packed_vectors = num_elements / 16;
    int num_unpacked_vectors = num_elements - num_packed_vectors * 16;
	if(num_elements > 15){
	__m512i v1[num_packed_vectors];
	__m512i v2[num_packed_vectors];
	__m512i v3[num_packed_vectors];

	__m512i x1 = _mm512_set1_epi32(v[1]); //Broadcast v[1] to all elements of v1
	__m512i x2 = _mm512_set1_epi32(v[2]); //Broadcast v[2] to all elements of v2
	__m512i x3 = _mm512_set1_epi32(v[3]); //Broadcast v[2] to all elements of v2
	
	
	for (i =0; i < num_packed_vectors; i++) {
		int base = i * 16;		
        //init dim 1
		v1[i] = _mm512_set_epi32(
		n->Cl->g[base + 0][1],
		n->Cl->g[base + 1][1],
		n->Cl->g[base + 2][1],
		n->Cl->g[base + 3][1],
		n->Cl->g[base + 4][1],
		n->Cl->g[base + 5][1],
		n->Cl->g[base + 6][1],
		n->Cl->g[base + 7][1],
		n->Cl->g[base + 8][1],
		n->Cl->g[base + 9][1],
		n->Cl->g[base + 10][1],
		n->Cl->g[base + 11][1],
		n->Cl->g[base + 12][1],
		n->Cl->g[base + 13][1],
		n->Cl->g[base + 14][1],
		n->Cl->g[base + 15][1]);
        //init dim 2
		v2[i] = _mm512_set_epi32(
		n->Cl->g[base + 0][2],
		n->Cl->g[base + 1][2],
		n->Cl->g[base + 2][2],
		n->Cl->g[base + 3][2],
		n->Cl->g[base + 4][2],
		n->Cl->g[base + 5][2],
		n->Cl->g[base + 6][2],
		n->Cl->g[base + 7][2],
		n->Cl->g[base + 8][2],
		n->Cl->g[base + 9][2],
		n->Cl->g[base + 10][2],
		n->Cl->g[base + 11][2],
		n->Cl->g[base + 12][2],
		n->Cl->g[base + 13][2],
		n->Cl->g[base + 14][2],
		n->Cl->g[base + 15][2]);
		v3[i] = _mm512_set_epi32(
		n->Cl->g[base + 0][3],
		n->Cl->g[base + 1][3],
		n->Cl->g[base + 2][3],
		n->Cl->g[base + 3][3],
		n->Cl->g[base + 4][3],
		n->Cl->g[base + 5][3],
		n->Cl->g[base + 6][3],
		n->Cl->g[base + 7][3],
		n->Cl->g[base + 8][3],
		n->Cl->g[base + 9][3],
		n->Cl->g[base + 10][3],
		n->Cl->g[base + 11][3],
		n->Cl->g[base + 12][3],
		n->Cl->g[base + 13][3],
		n->Cl->g[base + 14][3],
		n->Cl->g[base + 15][3]);

        __mmask32 mask1 = _mm512_cmp_epi32_mask(v1[i], x1, _MM_CMPINT_LE);
        __mmask32 mask2 = _mm512_cmp_epi32_mask(v2[i], x2, _MM_CMPINT_LE);
        __mmask32 mask3 = _mm512_cmp_epi32_mask(v3[i], x3, _MM_CMPINT_LE);

		//now perform *and* over results 
		__mmask32 and_res = mask1 & mask2 & mask3;
		unsigned int result = _cvtmask32_u32(and_res);
		if (result != 0){    //test if and_res contains a positive vale
			return 1;
		}

	}
	}
	
	for (i = num_packed_vectors * 16; i< n->Cl->cont;i++){
		stat_gcloperations++;		
		if (isDominated(v, n->Cl->g[i])){
			return 1;
		}
	}
	return 0;
}


void PacketAdd3(gnode *n) 
{
	goalv1[goalnp] = _mm512_set_epi32(
	n->Cl->g[0][1],
	n->Cl->g[1][1],
	n->Cl->g[2][1],
	n->Cl->g[3][1],
	n->Cl->g[4][1],
	n->Cl->g[5][1],
	n->Cl->g[6][1],
	n->Cl->g[7][1],
	n->Cl->g[8][1],
	n->Cl->g[9][1],
	n->Cl->g[10][1],
	n->Cl->g[11][1],
	n->Cl->g[12][1],
	n->Cl->g[13][1],
	n->Cl->g[14][1],
	n->Cl->g[15][1]);
	goalv2[goalnp] = _mm512_set_epi32(
	n->Cl->g[0][2],
	n->Cl->g[1][2],
	n->Cl->g[2][2],
	n->Cl->g[3][2],
	n->Cl->g[4][2],
	n->Cl->g[5][2],
	n->Cl->g[6][2],
	n->Cl->g[7][2],
	n->Cl->g[8][2],
	n->Cl->g[9][2],
	n->Cl->g[10][2],
	n->Cl->g[11][2],
	n->Cl->g[12][2],
	n->Cl->g[13][2],
	n->Cl->g[14][2],
	n->Cl->g[15][2]);
	n->Cl->cont = 0;
	goalnp++;
}

void PacketAdd4(gnode *n) 
{

	goalv1[goalnp] = _mm512_set_epi32(
	n->Cl->g[0][1],
	n->Cl->g[1][1],
	n->Cl->g[2][1],
	n->Cl->g[3][1],
	n->Cl->g[4][1],
	n->Cl->g[5][1],
	n->Cl->g[6][1],
	n->Cl->g[7][1],
	n->Cl->g[8][1],
	n->Cl->g[9][1],
	n->Cl->g[10][1],
	n->Cl->g[11][1],
	n->Cl->g[12][1],
	n->Cl->g[13][1],
	n->Cl->g[14][1],
	n->Cl->g[15][1]);
	goalv2[goalnp] = _mm512_set_epi32(
	n->Cl->g[0][2],
	n->Cl->g[1][2],
	n->Cl->g[2][2],
	n->Cl->g[3][2],
	n->Cl->g[4][2],
	n->Cl->g[5][2],
	n->Cl->g[6][2],
	n->Cl->g[7][2],
	n->Cl->g[8][2],
	n->Cl->g[9][2],
	n->Cl->g[10][2],
	n->Cl->g[11][2],
	n->Cl->g[12][2],
	n->Cl->g[13][2],
	n->Cl->g[14][2],
	n->Cl->g[15][2]);
	goalv3[goalnp] = _mm512_set_epi32(
	n->Cl->g[0][3],
	n->Cl->g[1][3],
	n->Cl->g[2][3],
	n->Cl->g[3][3],
	n->Cl->g[4][3],
	n->Cl->g[5][3],
	n->Cl->g[6][3],
	n->Cl->g[7][3],
	n->Cl->g[8][3],
	n->Cl->g[9][3],
	n->Cl->g[10][3],
	n->Cl->g[11][3],
	n->Cl->g[12][3],
	n->Cl->g[13][3],
	n->Cl->g[14][3],
	n->Cl->g[15][3]);
	n->Cl->cont = 0;
	goalnp++;
}


int does_dominate_goal3(gnode *n, unsigned *v ){
	__m512i x1 = _mm512_set1_epi32(v[1]); //Broadcast v[1] to all elements of v1
	__m512i x2 = _mm512_set1_epi32(v[2]); //Broadcast v[2] to all elements of v2
	
	for (int i = 0; i < goalnp; i++) {
		//compare each dimension            
        __mmask32 mask1 = _mm512_cmp_epi32_mask(goalv1[i], x1, _MM_CMPINT_LE);
        __mmask32 mask2 = _mm512_cmp_epi32_mask(goalv2[i], x2, _MM_CMPINT_LE);

		//now perform *and* over three results 
		__mmask32 and_res = /*mask0 &*/ mask1 & mask2;
		unsigned int result = _cvtmask32_u32(and_res);
		if (result != 0){    //test if and_res contains a positive vale
			return 1;
		}
	}
	return 0;
}

int does_dominate_goal4(gnode *n, unsigned *v ){

	__m512i x1 = _mm512_set1_epi32(v[1]); //Broadcast v[1] to all elements of v1
	__m512i x2 = _mm512_set1_epi32(v[2]); //Broadcast v[2] to all elements of v2
	__m512i x3 = _mm512_set1_epi32(v[3]); //Broadcast v[2] to all elements of v2
	
	for (int i = 0; i < goalnp; i++) {
		//compare each dimension            
        __mmask32 mask1 = _mm512_cmp_epi32_mask(goalv1[i], x1, _MM_CMPINT_LE);
        __mmask32 mask2 = _mm512_cmp_epi32_mask(goalv2[i], x2, _MM_CMPINT_LE);
        __mmask32 mask3 = _mm512_cmp_epi32_mask(goalv3[i], x3, _MM_CMPINT_LE);

		//now perform *and* over three results 
		__mmask32 and_res = mask1 & mask2 & mask3;
		unsigned int result = _cvtmask32_u32(and_res);
		if (result != 0){    //test if and_res contains a positive vale
			return 1;
		}
	}
	return 0;
}


int does_dominate_goal(gnode *n, unsigned *v ){
	int dom;
	if (nobjs == 3)
		dom = does_dominate_goal3(n, v);
	else
		dom = does_dominate_goal4(n, v);
	return dom;
}

void PacketAddGoal(gnode *n){
	if (nobjs == 3)
		PacketAdd3(n);
	else
		PacketAdd4(n);
}


//*************** End Vectorized Functions ***********************************//


int GoalCheck(gnode *n, unsigned *v)
{
	int i;

	if(does_dominate_goal(n, v))
		return 1;
	for (i = 0; i< n->Cl->cont;i++){
		stat_gcloperations++;		
		if (isDominated(v, n->Cl->g[i])){
			return 1;
		}
	}
    return 0;
}

int NodesCheck(gnode *n, unsigned *v){
	int dom;
	if (nobjs == 3)
		dom = PacketAddAll_Check_Dominance3(n, v);
	else
		dom = PacketAddAll_Check_Dominance4(n, v);
	return dom;
} 

void simpleAdd(gnode *n, unsigned *v, int s)
{
	int k;
	for (k = 0; k < nobjs;k++)
		n->Cl->g[n->Cl->cont][k] = v[k]; 
	n->Cl->cont++;
	if (n->id == goal && n->Cl->cont == 16)
		PacketAddGoal(n);
	stat_gclupdate++;	
}


int moastarv() {
    nsolutions = 0;
    stat_pruned = 0;

    emptyheap();

    start_node = new_node();
    ++stat_created;
    start_node->state = start;
    int i;
    for (i = 0; i < nobjs; i++)
		start_node->g[i] = 0;
    start_node->searchtree = NULL;
    insertheap(start_node);

    stat_expansions = 0;
    
	unsigned f[MAXOBJ-1],g[MAXOBJ-1];
    while (topheap() != NULL) {
        snode* n = popheap(); //best node in open
        short d;  
 		for (i = 0; i< nobjs;i++){
			f[i] = n->g[i]+graph_node[n->state].h[i];
			g[i] = n->g[i];
		}
		if (GoalCheck((&graph_node[goal]),f) || NodesCheck((&graph_node[n->state]),g)){
            stat_pruned++;
            if (next_recycled < MAX_RECYCLE) {
                recycled_nodes[next_recycled++] = n;
            }
            continue;
        }
		simpleAdd((&graph_node[n->state]),g,n->state);
		
        if (n->state == goal){
		 //  for (i = 0; i< nobjs;i++)
		//		printf("%d ",n->g[i]);
          //  printf("nsolutions:%d expanded:%llu generated:%llu heapsize:%d pruned:%d perc:%llu\n"
           // , nsolutions, stat_expansions, stat_generated, sizeheap(), stat_pruned,stat_percolations);
           
            nsolutions++;
            if (nsolutions > MAX_SOLUTIONS) {
                printf("Maximum number of solutions reached, increase MAX_SOLUTIONS!\n");
                exit(1);
            }
            if (next_recycled < MAX_RECYCLE) {
                recycled_nodes[next_recycled++] = n;
            }
            continue;
		}
        ++stat_expansions;
        for (d = 1; d < adjacent_table[n->state][0] * (nobjs+1); d += (nobjs+1)) {
            snode* succ;
            int i;
            unsigned nsucc = adjacent_table[n->state][d];
            succ = new_node();
            succ->state = nsucc;
            succ->searchtree = n;
            stat_generated++;
            for (i=0;i<nobjs;i++){
				succ->g[i] = n->g[i] + adjacent_table[n->state][d + (i+1)];
				g[i] = succ->g[i];
				f[i] = succ->g[i] + graph_node[succ->state].h[i]; 
			}
			if (GoalCheck((&graph_node[goal]),f) || NodesCheck((&graph_node[nsucc]),g)){
				stat_pruned++;
				if (next_recycled < MAX_RECYCLE) {
					recycled_nodes[next_recycled++] = succ;
				}
            continue;
			}
			insertheap(succ);
        }
    }
    return nsolutions > 0;
}

/* ------------------------------------------------------------------------------*/
void call_moastarv() {
    float runtime;
    struct timeval tstart, tend;

    initialize_parameters();

    gettimeofday(&tstart, NULL);

    //Dijkstra    
    int i;
    for (i = 1; i<= nobjs;i++){ 
		backward_dijkstra(i);
	}

    //MOA*v
    moastarv();

    gettimeofday(&tend, NULL);
    runtime = 1.0 * (tend.tv_sec - tstart.tv_sec) + 1.0 * (tend.tv_usec - tstart.tv_usec) / 1000000.0;
    printf("%lld;%lld;%d;%f:%llu;%llu;%d;%llu;%llu;%d;%llu;%d;%s\n",
        start_state->id + 1,
        goal_state->id + 1,
        nsolutions,
        runtime * 1000,
        stat_generated,
        stat_expansions,
        stat_created,
        stat_gcloperations,
        stat_gclupdate,
        stat_pruned,
        stat_percolations,
        graph_node[start].h[0],
        filename);
}
