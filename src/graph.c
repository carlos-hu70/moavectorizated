#include "include.h"
#include "moastarv.h"
#include "graph.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>


void read_adjacent_table(const char* filename) {
	FILE* f;
	int i,j, ori, dest;
	f = fopen(filename, "r");
	int num_arcs = 0;
	if (f == NULL) 	{
		printf("Cannot open file %s.\n", filename);
		exit(1);
	}
	fscanf(f, "%d %d", &num_gnodes, &num_arcs);
	fscanf(f, "\n");

	for (i = 0; i < num_gnodes; i++)
		adjacent_table[i][0] = 0;

	for (i = 0; i < num_arcs; i++) {
		fscanf(f, "%d %d", &ori, &dest);
		adjacent_table[ori - 1][0]++;
		adjacent_table[ori - 1][adjacent_table[ori - 1][0] * (nobjs+1) - nobjs] = dest - 1;
		
		pred_adjacent_table[dest - 1][0]++;
		pred_adjacent_table[dest - 1][pred_adjacent_table[dest - 1][0] * (nobjs+1) - nobjs] = ori - 1;

		for (j = (nobjs-1); j >=0; j--){
			int c;
			fscanf(f, "%d", &c);
			adjacent_table[ori - 1][adjacent_table[ori - 1][0] * (nobjs+1) - j] = c;
			pred_adjacent_table[dest - 1][pred_adjacent_table[dest - 1][0] * (nobjs+1) - j] = c;
		}
		fscanf(f, "\n");
	}
	
	
	fclose(f);
}

void new_graph() {
	int y,i;

	if (graph_node == NULL) {
		graph_node = (gnode*) calloc(num_gnodes, sizeof(gnode));
		for (y = 0; y < num_gnodes; ++y) 		{
			nodeCl *CLlist = NULL;
			CLlist = (nodeCl *)malloc(sizeof(nodeCl));
			graph_node[y].id = y;
			graph_node[y].gmin = LARGE;
			graph_node[y].h1 = LARGE;
			graph_node[y].h2 = LARGE;
			for (i = 0; i < nobjs; i++){
				graph_node[y].h[i] = LARGE;
			}
			CLlist->cont = 0;
			graph_node[y].Cl = CLlist;
		}
	}
}
