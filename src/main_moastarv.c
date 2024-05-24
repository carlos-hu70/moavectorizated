/////////////////////////////////////////////////////////////////////
// Carlos Hernandez
// All rights reserved
/////////////////////////////////////////////////////////////////////

#include "include.h"
#include "moastarv.h"
#include "graph.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>

/*----------------------------------------------------------------------------------*/
int main(int argc, char** argv) {
	if (argc != 5) 	{
		printf("Usage: %s [graph_file] [start_node] [goal_node] [n_objectives]\n", argv[0]);
		exit(1);
	}
	strcpy(filename, argv[1]);
	start = atoi(argv[2]) - 1;
	goal = atoi(argv[3]) - 1;
	nobjs = atoi(argv[4]);
	read_adjacent_table(filename);
	new_graph();
	call_moastarv();
}
