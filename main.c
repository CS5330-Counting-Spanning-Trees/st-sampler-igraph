#include <igraph/igraph.h>

#include <time.h>
#include <limits.h>
// https://www.techiedelight.com/find-execution-time-c-program/

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
void get_random_graph(igraph_t* g, igraph_integer_t n, double density, igraph_integer_t max_degree, igraph_integer_t min_degree)
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
			

			// THIS IMPLEMENTATION IS TOO SLOW, DUE TO CALLING igraph_add_edge INDIVIDUALLY
			// // store the degrees for i and j
			// igraph_vector_t deg;
			// igraph_vector_init(&deg, 2);

			// // obtain degrees for vertices i and j
			// igraph_vs_t vs;
			// igraph_vs_vector_small(&vs, i, j, -1);
			// igraph_degree(g, &deg, vs, IGRAPH_ALL, IGRAPH_LOOPS);

			// igraph_integer_t degi, degj;
			// degi = VECTOR(deg)[0];
			// degj = VECTOR(deg)[1];

			// igraph_vector_destroy(&deg);
			
			// if (degi >= max_degree || degj >= max_degree)
			// 	continue;
			
			// if (degi < min_degree || p <= density)
			// {
			// 	igraph_add_edge(g, i, j);
			// }

		}
	}
	igraph_adjlist(g, &adjlist, IGRAPH_ALL, 1);
}

void generate_random_connected_graph(igraph_t* g, igraph_integer_t n, double density, igraph_integer_t max_degree, igraph_integer_t min_degree)
{

	while(1)
	{
		get_random_graph(g, n, density, max_degree, min_degree);
		igraph_bool_t is_connected;
		igraph_is_connected(g, &is_connected, IGRAPH_STRONG);
		if(is_connected)
			break;
		// clean up for the next round
		igraph_destroy(g);
		printf("Generation of connected graph failed, retry...\n");
	}

}

int main(void)
{
	srand(time(NULL)); // Initialization, should only be called once.

	igraph_t g;

	// generate_small_test_graph(&g);
	generate_random_connected_graph(&g, 10000, 0.1, 5, 0);
	// igraph_full(&g, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
	igraph_ring(&g,10000,IGRAPH_UNDIRECTED, 0, 1);

	printf("Created graph with %d vertices and %d edges\n", igraph_vcount(&g), igraph_ecount(&g));

	igraph_vector_t res;
	igraph_vector_init(&res,0);

	int i = 0;
	long begin = clock();
	for (; i < 1000; i++)
	{
		igraph_random_spanning_tree(&g, &res, 0);

		if (i%100 == 0)
			printf("Generated sample %d\n", i+1);
	}
		


	long end = clock();

	printf("%d samples taken, with per sample time taking %.2lf ms\n", i, (double)(end - begin) / i / CLOCKS_PER_SEC * 1e3);
	printf("Total time spent %ld seconds\n\n", (end - begin) / CLOCKS_PER_SEC);
	// printf("returned vector size = %ld\n", igraph_vector_size(&res));

	igraph_vector_destroy(&res);
	igraph_destroy(&g);
	return 0;
}