

#include <time.h>
// https://www.techiedelight.com/find-execution-time-c-program/
#include <graph_generator.h>
#include <approx_count_st.h>

int main(void)
{
	srand(time(NULL)); // Initialization, should only be called once.

	igraph_t g;

	// generate_small_test_graph(&g);
	
	generate_random_connected_graph(&g, 1000, 0.1, INT_MAX, 0);
	
	// igraph_full(&g, 10000, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

	// igraph_ring(&g,1000,IGRAPH_UNDIRECTED, 0, 1);

	// igraph_vector_t dimvector;
	// igraph_vector_init(&dimvector, 3);
	// VECTOR(dimvector)[0]=33;
	// VECTOR(dimvector)[1]=33;
	// VECTOR(dimvector)[2]=33;
	// igraph_lattice(&g, &dimvector, 0, IGRAPH_UNDIRECTED, 0,1);

	printf("Created graph with %d vertices and %d edges\n", igraph_vcount(&g), igraph_ecount(&g));

	// igraph_vector_t res;
	// igraph_vector_init(&res,0);

	// int i = 0;
	// long begin = clock();
	// for (; i < 1000; i++)
	// {
	// 	igraph_random_spanning_tree(&g, &res, 0);

	// 	if (i%100 == 0)
	// 		printf("Generated sample %d\n", i+1);
	// }



	// long end = clock();

	// printf("%d samples taken, with per sample time taking %.2lf ms\n", i, (double)(end - begin) / i / CLOCKS_PER_SEC * 1e3);
	// printf("Total time spent %ld seconds\n\n", (end - begin) / CLOCKS_PER_SEC);
	// printf("returned vector size = %ld\n", igraph_vector_size(&res));
	// igraph_vector_destroy(&res);

	approx_count_result_t res = approx_count_st(&g);

	
	igraph_destroy(&g);
	return 0;
}