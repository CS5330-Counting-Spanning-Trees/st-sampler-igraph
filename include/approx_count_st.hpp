#pragma once

#include <igraph/igraph.h>

#include <cassert>
#include <vector>

#include <array>

#include <algorithm>

const int PRESAMPLE_SIZE_REQUIRED = 50;
const int PIVOT_BUFFER_SIZE = 4;
const double FLUCTURATION_THRESHOLD = 0.003;

class RandomSpanningTrees;

class ApproxCountST{

public:

    enum count_mode_t{
        UNSPECIFIED = 0,
        PRESENCE,
        ABSENCE
    };

    enum random_walk_mode_t{
        NORMAL = 0, // visit one means visit all, but the walk within the jointed vertices are permitted
        AS_ONE_VERTEX // think of the set vertices jointed by present edges as 1 single vertex, while doing the walk
    };

    typedef struct count_result{
        double count = 1.0;
        long long effective_samples = 0;
        long long actual_samples = 0;
        double epsilon; // probabilistic multiplicative error bound
        double delta; // probabilistic confidence
    }result_t;

    
    typedef struct pivot_stats{
        igraph_integer_t eid = -1;
        int total = 0;
        int rippled_total = 0;
        int present = 0;
        int absent = 0;
        count_mode_t count_mode = UNSPECIFIED;
        double ratio = -1;

        int update_count = 0;
        int requested_batch_size = 50;
        std::array<double, PIVOT_BUFFER_SIZE> ratio_buffer;

        void update(int new_present = 0, int new_absent = 0){
            present += new_present;
            absent += new_absent;
            total = present + absent;

            switch (count_mode){
                case PRESENCE:
                    ratio = present / (double)total;
                    break;
                case ABSENCE:
                    ratio = absent / (double)total;
                    break;
                default:
                    return;
            }

            
            ratio_buffer[update_count % PIVOT_BUFFER_SIZE] = ratio;
            update_count++;
        }

        bool converged(){
            if(count_mode == UNSPECIFIED)
                return false;
            
            // only test convergence when the buffer is fully filled
            if(update_count == 0 || update_count % PIVOT_BUFFER_SIZE)
                return false;

            double rmax = *std::max_element(ratio_buffer.begin(), ratio_buffer.end());
            double rmin = *std::min_element(ratio_buffer.begin(), ratio_buffer.end());
            
            assert(rmax > 0.2);
            // printf("%lf\n", rmin);
            assert(rmin > 0.1);

            // printf("rmax / rmin = %.6lf \t rmax %.3lf, rmin %.3lf \n", rmax / rmin , rmax, rmin);
            if (rmax / rmin < 1 + FLUCTURATION_THRESHOLD)
            {
                printf("\nFinal Ratio = %.3lf, batch size = %d\n\n", 0.5 * (rmax + rmin), requested_batch_size);
                return true; // CONVERGENCE!
            }

            printf("->%.3lf", rmax / rmin);

            // increase the sampling size, capped for 20% increase
            requested_batch_size = std::max(requested_batch_size, (total - rippled_total) / PIVOT_BUFFER_SIZE / 5);
            return false;
        }


        // This function should be called everytime new samples are collected, but the count_mode is still UNSPECIFIED
        bool try_set_count_mode(){
            if (total < PRESAMPLE_SIZE_REQUIRED)
                return false;

            if (present / (double)total > 0.5)
                count_mode = PRESENCE;
            else
                count_mode = ABSENCE;

            return true;
        }

        void print(){
            printf("total %d (%d|%d), ratio %.3lf, mode %d, update_count %d\n", total, present, absent, ratio, count_mode, update_count);
        }
    }pivot_stats_t;

    ApproxCountST() = delete;
    // Initialise constants and g_contracted
    ApproxCountST(igraph_t* g, random_walk_mode_t rw_mode);


    // result stored in ps_vec
    result_t approx_count_st();

    void print_all();

    std::vector<pivot_stats_t> ps_vec;

private:

    // This routine will make n ST samples, rooted at vertex vid. It also updates ps_vec until the first unspecified entries
    void sample_mini_batch_with_updates(RandomSpanningTrees* rst, int k_start);

    random_walk_mode_t random_walk_mode;
    igraph_t* g;
    

    const igraph_integer_t N;
	const igraph_integer_t M;
    int K;

};