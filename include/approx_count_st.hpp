#pragma once

#include "graph_lite.hpp"

#include <cassert>
#include <vector>

#include <array>

#include <algorithm>

/*
const int PRESAMPLE_SIZE_REQUIRED = 50;
const int PIVOT_BUFFER_SIZE = 8;
const double RATIO_FLUCTURATION_THRESHOLD = 0.002;
const double VARIANCE_THRESHOLD = 0.001;
const int INITIAL_REQUESTED_BATCH_SIZE = 500;
*/

class RandomSpanningTrees;

class ApproxCountST{

public:

    enum convergence_mode_t{
        RATIO = 0,
        VARIANCE,
        CONSTANT
    };

    enum count_mode_t{
        UNSPECIFIED = 0,
        PRESENCE,
        ABSENCE
    };

    
    typedef struct count_result{
        double count = 1.0;
        double count_log = 0.0;
        long long effective_samples = 0;
        long long actual_samples = 0;
        double epsilon; // probabilistic multiplicative error bound
        double delta; // probabilistic confidence
    }result_t;
    
    typedef struct pivot_stats{
        eid_t eid = -1;
        int total = 0;
        int rippled_total = 0; // record keeping of how many samples we gotten from rippling
        int present = 0;
        int absent = 0;
        count_mode_t count_mode = UNSPECIFIED;
        double ratio = -1;

        int update_count = 0;
        int requested_batch_size = INITIAL_REQUESTED_BATCH_SIZE;
        std::array<double, PIVOT_BUFFER_SIZE> inverse_ratio_buffer;


        void update(){

            // not enough samples to create a new ratio data in the buffer
            if (present + absent < total + requested_batch_size)
                return;
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

            
            inverse_ratio_buffer[update_count % PIVOT_BUFFER_SIZE] = 1.0 / ratio;
            update_count++;
        }

        bool test_constant(){
            if (total > convergence_constant_threshold)
                return true;
            else
                return false;
        }

        bool test_converged_ratio(){
            double rmax = *std::max_element(inverse_ratio_buffer.begin(), inverse_ratio_buffer.end());
            double rmin = *std::min_element(inverse_ratio_buffer.begin(), inverse_ratio_buffer.end());
            
            assert(rmax > 0.2);
            // printf("%lf\n", rmin);
            assert(rmin > 0.1);

           if (rmax / rmin < 1 + convergence_ratio_threshold){
                printf("\nFinal Ratio = %.3lf, batch size = %d\n\n", 0.5 * (rmax + rmin), requested_batch_size);
                ratio = 1.0 / (0.5 * (rmax + rmin));
                return true;
            }
            else
                return false;
        }

        bool test_converged_variance(){
            double sum = std::accumulate(inverse_ratio_buffer.begin(), inverse_ratio_buffer.end(), 0.0);
            double mean = sum / inverse_ratio_buffer.size();

            double sq_sum = std::inner_product(inverse_ratio_buffer.begin(), inverse_ratio_buffer.end(), inverse_ratio_buffer.begin(), 0.0);
            double stddev = std::sqrt(sq_sum / inverse_ratio_buffer.size() - mean * mean);

            if(stddev / mean < convergence_variance_threshold){
                ratio = 1.0 / mean;
                printf("\nFinal Ratio = %.3lf, batch size = %d\n\n", mean, requested_batch_size);
                return true;
            }
            else
                return false;
        }

        bool converged(){
            if(count_mode == UNSPECIFIED)
                return false;
            
            // only test convergence when the buffer is fully filled
            if(update_count < PIVOT_BUFFER_SIZE)
                return false;

            switch (convergence_mode){ // CONVERGENCE!
                case RATIO:
                    if(test_converged_ratio()) return true;
                    break;
                case VARIANCE:
                    if(test_converged_variance()) return true;
                    break;
                case CONSTANT:
                    if(test_constant()) return true;
                    break;
                default:
                    assert(0);
            };

            // printf("->%.3lf", rmax / rmin);

            // increase the sampling size, capped for 20% increase
            requested_batch_size = std::max(requested_batch_size, total / PIVOT_BUFFER_SIZE / 10);
            return false;
        }


        // This function should be called everytime new samples are collected, but the count_mode is still UNSPECIFIED
        bool try_set_count_mode(){
            if(count_mode != UNSPECIFIED)
                return true;
            
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
    ApproxCountST(GraphLite* g);


    // result stored in ps_vec
    result_t approx_count_st();

    void print_all();

    std::vector<pivot_stats_t> ps_vec;

    static convergence_mode_t convergence_mode;
    static double convergence_ratio_threshold;
    static double convergence_variance_threshold;
    static int convergence_constant_threshold;
    
private:
    inline bool check_convergence(eid_t* k);

    typedef struct sampling_struct{
        std::vector<eid_t> path;
        std::vector<eid_t> next;
        std::vector<bool> in_tree;
    }sampling_struct_t;

    // This routine will make n ST samples, rooted at vertex vid. It also updates ps_vec until the first unspecified entries
    void sample_mini_batch_with_updates(RandomSpanningTrees* rst, int k_start, sampling_struct_t* sampling_struct);

    GraphLite* gl;
    

    const vid_t N_initial;
	const eid_t M_initial;
    const int K;

    std::vector<eid_t> e_shuffle;

    
};
