

#include <time.h>
// https://www.techiedelight.com/find-execution-time-c-program/
#include <graph_generator.h>
#include <approx_count_st.hpp>

#include <graph_lite.hpp>

#include <stdio.h>


void log2file(FILE* fp, const char *__restrict __format, ...)
{
	va_list args;
    va_start(args, __format);
	// printf(__format, ...);
	vprintf(__format, args);
	va_end(args);
	va_start(args, __format);
	vfprintf(fp,__format, args);
	va_end(args);
}

int main(void)
{
	FILE *fp;
	fp = fopen("output.txt", "a");

	srand(123);
	// srand(time(NULL)); // Initialization, should only be called once.

	// Use igraph library for generation purpose only
	
	igraph_t g;

	// generate_small_test_graph(&g);
	
	generate_random_connected_graph(&g, 1000, 0.1, 5);
	
	// igraph_full(&g, 500, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

	// igraph_ring(&g,200,IGRAPH_UNDIRECTED, 0, 1);


	// igraph_vector_t dimvector;
	// igraph_vector_init(&dimvector, 3);
	// VECTOR(dimvector)[0]=10;
	// VECTOR(dimvector)[1]=10;
	// VECTOR(dimvector)[2]=10;
	// igraph_lattice(&g, &dimvector, 0, IGRAPH_UNDIRECTED, 0,1);


	GraphLite gl(&g);

	// gl.print();

	igraph_destroy(&g);


	time_t current_time;
    char* c_time_string;
	/* Convert to local time format. */
	current_time = time(NULL);
    c_time_string = ctime(&current_time);
	log2file(fp, "%s", c_time_string);

	log2file(fp, "Created graph with %d vertices and %d edges\n", gl.vertex_count_all() , gl.edge_count_all());

	fflush(fp);

	long begin = clock();

	//// Actual Counting work Starts
	ApproxCountST::result_t res;
	ApproxCountST* ast = new ApproxCountST(&gl); // putting on heap is needed, otherwise double free error
	res = ast->approx_count_st();

	//// Actual Counting work Ends

	long end = clock();

	log2file(fp,"%lld samples taken, with per sample time taking %.3lf ms\n", res.actual_samples, (double)(end - begin) / res.actual_samples / CLOCKS_PER_SEC * 1e3);
	log2file(fp,"Total time spent %ld seconds\n\n", (end - begin) / CLOCKS_PER_SEC);

	log2file(fp,"FINAL result = %.4e (e^%.4e) with %lld effective samples, avg %d samples per edge. \n\n", res.count, res.count_log, res.effective_samples, res.effective_samples / gl.edge_count_all());

	fclose(fp);

	// ast.print_all();
	
	// igraph_destroy(&g);
	return 0;
}