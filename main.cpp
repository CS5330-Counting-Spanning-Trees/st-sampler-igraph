

#include <time.h>
// https://www.techiedelight.com/find-execution-time-c-program/
#include <graph_generator.h>
#include <approx_count_st.hpp>

#include <graph_lite.hpp>

#include <stdio.h>

#include <mtt.hpp>


ApproxCountST::convergence_mode_t ApproxCountST::convergence_mode = CONVERGENCE_MODE;
double ApproxCountST::convergence_ratio_threshold = RATIO_THRESHOLD_DEFAULT;
double ApproxCountST::convergence_variance_threshold = VARIANCE_THRESHOLD_DEFAULT;
int ApproxCountST::convergence_constant_threshold = CONSTANT_THRESHOLD_DEFAULT;

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

int main(int argc, char* argv[])
{

	if (argc < 4) {
        // Tell the user how to run the program
       printf("Usage: %s NUM_OF_VERTICES NUM_OF_LOOPS THRESHOLD\n", argv[0] );
        /* "Usage messages" are a conventional way of telling the user
         * how to run a program if they enter the command incorrectly.
         */
        return 1;
    }

	const int N = atoi(argv[1]);

	int L = atoi(argv[2]);

	const double threshold = atof(argv[3]);

	ApproxCountST::convergence_constant_threshold =
	ApproxCountST::convergence_ratio_threshold =
	ApproxCountST::convergence_variance_threshold = threshold;


	FILE *fp;
	char filename[200];
	strcpy(filename, argv[0]);
	fp = fopen(strcat(filename,".txt"), "a");

	log2file(fp, "-------------------------------------------------------------\n");

	srand(123);
	// srand(time(NULL)); // Initialization, should only be called once.

	// Use igraph library for generation purpose only
	
	igraph_t g;

	// generate_small_test_graph(&g);
	
#ifdef SPARSE_GRAPH
	generate_random_connected_graph(&g, N, 0.1, 5);
#endif

#ifdef FULL_GRAPH
	igraph_full(&g, N, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
#endif

#ifdef RING_GRAPH
	igraph_ring(&g,N,IGRAPH_UNDIRECTED, 0, 1);
#endif

	// igraph_vector_t dimvector;
	// igraph_vector_init(&dimvector, 3);
	// VECTOR(dimvector)[0]=10;
	// VECTOR(dimvector)[1]=10;
	// VECTOR(dimvector)[2]=10;
	// igraph_lattice(&g, &dimvector, 0, IGRAPH_UNDIRECTED, 0,1);

	GraphLite gl(&g);

	time_t current_time;
	struct tm * timeinfo;
    char* c_time_string;
	/* Convert to local time format. */
	current_time = time(NULL);
	timeinfo = localtime(&current_time);
    c_time_string = ctime(&current_time);
	log2file(fp, "%s", c_time_string);
	log2file(fp, "Created graph with %d vertices and %d edges\n", gl.vertex_count_all() , gl.edge_count_all());
	fflush(fp);


	// for writing csv file
	FILE *fp_csv;
	char csv_filename[200];
	strftime(filename,30,"%m-%d-%H-%M-%S",timeinfo);
	sprintf(csv_filename, "%s-%s.csv" ,argv[0], filename);
	fp_csv = fopen(csv_filename, "w");

	long begin = clock();

	printf("Calculating Laplacian Matrix...\n"); 
	igraph_matrix_t laplacian_matrix; // igraph uses row major, but doesnt matter for laplacian
	igraph_matrix_init(&laplacian_matrix, N, N);
	igraph_laplacian(&g, &laplacian_matrix, nullptr, 0, nullptr);

	printf("Calculating MTT count...\n");

	Eigen::Map<Eigen::MatrixXd> g_eigen(laplacian_matrix.data.stor_begin,N,N);

	// double det_value = g_eigen.topRightCorner(N-1,N-1).determinant();
	// log2file(fp, "\nMTT Result = e^%.4e\n", det_value);

	double logdet_value = logdet(g_eigen.topLeftCorner(N-1,N-1), true);

	log2file(fp, "MTT Result = (e^%.4e), ", logdet_value);

	long end = clock();

	log2file(fp,"Total time spent for MTT: %ld seconds\n\n", (end - begin) / CLOCKS_PER_SEC);

	igraph_destroy(&g);

	//////////////////////////////////////////////////////////////////////////////////////////

	//// Logging Parameters

	char params[1024];
	sprintf(params,"N=%d, M=%d, presample_size=%d, buffer_size=%d, convergence_mode=%d, threshold=%lf, initial_batch_size=%d\n", N, gl.edge_count_all(), PRESAMPLE_SIZE_REQUIRED, PIVOT_BUFFER_SIZE, CONVERGENCE_MODE, threshold, INITIAL_REQUESTED_BATCH_SIZE);

	log2file(fp, "%s", params);
	fprintf(fp_csv,"%s", params);

	ApproxCountST::result_t res;
	
	fprintf(fp_csv,"trial, count_log, mtt_log, time_spent, error rate, actual samples\n");
	log2file(fp,"stats writting to file %s\n", csv_filename);

	fflush(fp_csv);
	fflush(fp);

	for (int l = 0; l < L; l++){

		GraphLite gl_temp = gl;

		log2file(fp,"<<<<<<<<<<<ROUND %d<<<<<<<<<<<<<<\n", l+1);
		begin = clock();
		ApproxCountST* ast = new ApproxCountST(&gl_temp); // putting on heap is needed, otherwise double free error
		res = ast->approx_count_st();

		end = clock();

		log2file(fp,"%lld actual samples taken, with per sample time taking %.3lf ms\n", res.actual_samples, (double)(end - begin) / res.actual_samples / CLOCKS_PER_SEC * 1e3);
		

		log2file(fp,"ROUND %d FINAL result = %.4e (e^%.4e) with %lld effective samples, avg %d samples per edge. \n", 
			l+1, res.count, res.count_log, res.effective_samples, res.effective_samples / gl.edge_count_all());

		log2file(fp,"error percentage %.2lf%%", 100.0 * (std::exp(res.count_log - logdet_value) - 1.0) );
		log2file(fp,", time spent for randomised algo: %ld seconds\n\n", (end - begin) / CLOCKS_PER_SEC);

		fprintf(fp_csv, "%d, %.4e, %.4e, %ld, %.4lf, %lld\n", l+1, res.count_log, logdet_value, (end - begin) / CLOCKS_PER_SEC, (std::exp(res.count_log - logdet_value) - 1.0), res.actual_samples);
		fflush(fp_csv);
		fflush(fp);
		delete ast;
		ast = nullptr;
	}

	log2file(fp,"stats written to file %s.csv\n", csv_filename);
	
	fclose(fp);
	fclose(fp_csv);
	return 0;
}