# moavectorizated
#Multi-Objective A* with vectorized operations

This is the C implementation of LT-MOA* with vectorized operations(LT-MOA*v). LT-MOA*v computes a Pareto frontier containing paths from a given initial node to a given goal node on a given graph. The implementation is for 3 and 4 objectives. More details in the following paper:

* [Carlos Hernandez Ulloa, Han Zhang, Sven Koenig, Ariel Felner, Oren Salzman:
Efficient Set Dominance Checks in Multi-Objective Shortest-Path Algorithms via
Vectorized Operations. SOCS 2024]

The implementation contained in this package assumes the graph is explicitly given.

# Compilation

Compile using `make`.

# Running LT-MOA*v

To run a single problem use  `./moastarv [graph_file] [start_node] [goal_node] [n_objectives]`. For example: `./moastarv Maps/NY-road-m-t-d.txt 28493 5425 3` or `./moastarv Maps/NY-road-l-d-t-m.txt 36255 28938 4`. 