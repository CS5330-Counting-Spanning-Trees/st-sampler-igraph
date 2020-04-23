#ifndef GRAPH_GENERATOR_H
#define GRAPH_GENERATOR_H

#include <igraph/igraph.h>
#include <limits.h>

void generate_small_test_graph(igraph_t* g)
{
	igraph_vector_t v1;

	igraph_real_t edges[] = {	0, 1, 0, 3, 
								1, 2, 1, 3,
								2, 3, 2, 5,
								3, 4,
								4, 5
							};

	igraph_vector_view(&v1, edges, sizeof(edges)/sizeof(double));
	igraph_create(g, &v1, 0, IGRAPH_UNDIRECTED);
	
}

// assume uninitialised graph pointer
void generate_random_graph(igraph_t* g, igraph_integer_t n, double density, igraph_integer_t max_degree, igraph_integer_t min_degree)
{
	igraph_empty(g, n, IGRAPH_UNDIRECTED);
	// Using density to add edges randomly

	igraph_adjlist_t adjlist;
	igraph_adjlist_init(g, &adjlist, IGRAPH_OUT);

	for(igraph_integer_t i = 0; i < n; i++)
	{
		for(igraph_integer_t j = i + 1; j < n; j++) // no self-loop allowed
		{
			double p = rand() / (double)RAND_MAX;

			long degi = igraph_vector_int_size(&adjlist.adjs[i]);
			long degj = igraph_vector_int_size(&adjlist.adjs[j]);

			if (degi >= max_degree || degj >= max_degree)
				continue;
			
			if (degi < min_degree || p <= density)
			{
				igraph_vector_int_push_back(&adjlist.adjs[i], j);
				igraph_vector_int_push_back(&adjlist.adjs[j], i);
			}

		}
	}
	igraph_adjlist(g, &adjlist, IGRAPH_ALL, 1);
}

void generate_random_connected_graph(igraph_t* g, igraph_integer_t n, double density, igraph_integer_t max_degree, igraph_integer_t min_degree)
{

	while(1)
	{
		generate_random_graph(g, n, density, max_degree, min_degree);
		igraph_bool_t is_connected;
		igraph_is_connected(g, &is_connected, IGRAPH_STRONG);
		if(is_connected)
			break;
		// clean up for the next round
		igraph_destroy(g);
		printf("Generation of connected graph failed, retry...\n");
	}

}

#endif /* GRAPH_GENERATOR_H */
