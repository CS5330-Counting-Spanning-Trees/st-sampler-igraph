#ifndef APPROX_COUNT_ST_H
#define APPROX_COUNT_ST_H

#include <igraph/igraph.h>
#include <limits.h>

#include <assert.h>

inline double MAX(double a, double b) { return((a) > (b) ? a : b); }
inline double MIN(double a, double b) { return((a) < (b) ? a : b); }

typedef struct approx_count_result{
	int count;
	double epsilon; // probabilistic multiplicative error bound
	double delta; // probabilistic confidence
}approx_count_result_t;


inline void sample_mini_batch(igraph_t* g, int K, igraph_integer_t pivot, int* pivot_present, int* pivot_absent)
{
	igraph_vector_t res;
	igraph_vector_init(&res,0);

	for (int k = 0 ; k < K ; k++)
	{
		igraph_random_spanning_tree(g, &res, 0);
		if (igraph_vector_contains(&res, pivot)){
			(*pivot_present)++;
		}else{
			(*pivot_absent)++;
		}
	}
	
}


approx_count_result_t approx_count_st(igraph_t* g)
{
	// to store the carry-forward count from parent recursion
	int cf_pivot_present = 0;
	int cf_pivot_absent = 0;

	// the pivot used for forming self-reduction
	igraph_integer_t pivot;

	const igraph_integer_t N = igraph_vcount(g);
	const igraph_integer_t M = igraph_ecount(g);

	// Determine the number of times the ratio to be calculated, K
	assert(M >= N-1);
	const int K = M - (N-1);

	// Setup a circular buffer to store the last three
	const size_t PIVOT_BUFFER_SIZE = 4;
	double pivot_buffer[PIVOT_BUFFER_SIZE];

	// Initialise before the first round
	double accum_ratio = 1;
	pivot = rand() / (double) RAND_MAX * M;

	// This loop will calculate ratio of G_M / G_{M-1} ..... G_{N} / G_{N-1}, total of M-1 - (N-1) + 1 terms
	for(int k = 0; k < K; k++)
	{
		int pivot_present = 0, pivot_absent = 0;
		const size_t INITIAL_BATCH_SIZE = 20;
		size_t mini_K = INITIAL_BATCH_SIZE;

		int mini_iter = 1;
		double buffer_min = 1, buffer_max = 0;
		for(;;mini_iter++)
		{
			int batch_pivot_present = 0;
			int batch_pivot_absent = 0;
			sample_mini_batch(g, mini_K, pivot, &batch_pivot_present, &batch_pivot_absent);

			printf("mini batch %d %d\n", batch_pivot_present, batch_pivot_absent);

			pivot_present += batch_pivot_present;
			pivot_absent += batch_pivot_absent;

			double ratio = pivot_present / (double)(pivot_present + pivot_absent);
			pivot_buffer[k%PIVOT_BUFFER_SIZE] = ratio;

			printf("ratio = %.6lf\n", ratio);

			buffer_min = MIN(ratio, buffer_min);
			buffer_max = MAX(ratio, buffer_max);


			// buffer is filled for a complete round!
			if (mini_iter%PIVOT_BUFFER_SIZE == 0)
			{
				// now compare the samples in the bins of pivot_buffer
				
				if (0.5 * (buffer_max + buffer_min) < 0.5);
				if (buffer_max / buffer_min < 1.05)
					break; // CONVERGENCE!

				// Otherwise, continue with sample doubling
				mini_K = (pivot_present + pivot_absent) / PIVOT_BUFFER_SIZE;
				buffer_min = 1;
				buffer_max = 0;

			}

		}
		printf("Pivot %d/%d converged after %d iterations\n", k+1, K, mini_iter);
	}

	// Return the results
	approx_count_result_t ret;

	ret.count = accum_ratio;

	return ret;
}

#endif /* APPROX_COUNT_ST_H */
