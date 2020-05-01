#pragma once

#include "graph_lite.hpp"

#include <cassert>
#include <vector>

#include <array>

#include <algorithm>

const int PRESAMPLE_SIZE_REQUIRED = 50;
const int PIVOT_BUFFER_SIZE = 8;
const double RATIO_FLUCTURATION_THRESHOLD = 0.002;
const double VARIANCE_THRESHOLD = 0.001;

class RandomSpanningTrees;

class ApproxCountST{

public:

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
        int requested_batch_size = 50;
        std::array<double, PIVOT_BUFFER_SIZE> ratio_buffer;


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

            
            ratio_buffer[update_count % PIVOT_BUFFER_SIZE] = ratio;
            update_count++;
        }

        bool test_converged_ratio(){
            double rmax = *std::max_element(ratio_buffer.begin(), ratio_buffer.end());
            double rmin = *std::min_element(ratio_buffer.begin(), ratio_buffer.end());
            
            assert(rmax > 0.2);
            // printf("%lf\n", rmin);
            assert(rmin > 0.1);

            // printf("%.3lf %.3lf %.3lf %.3lf\n", ratio_buffer[0], ratio_buffer[1] ,ratio_buffer[2] ,ratio_buffer[3]);

            // printf("rmax / rmin = %.6lf \t rmax %.3lf, rmin %.3lf total %d \n", rmax / rmin , rmax, rmin, total);
            if (rmax / rmin < 1 + RATIO_FLUCTURATION_THRESHOLD){
                printf("\nFinal Ratio = %.3lf, batch size = %d\n\n", 0.5 * (rmax + rmin), requested_batch_size);
                return true;
            }
            else
                return false;
        }

        bool test_converged_variance(){
            double sum = std::accumulate(ratio_buffer.begin(), ratio_buffer.end(), 0.0);
            double mean = sum / ratio_buffer.size();

            double sq_sum = std::inner_product(ratio_buffer.begin(), ratio_buffer.end(), ratio_buffer.begin(), 0.0);
            double stddev = std::sqrt(sq_sum / ratio_buffer.size() - mean * mean);

            if(stddev / mean < VARIANCE_THRESHOLD){
                ratio = mean;
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

            if(test_converged_variance())
                return true; // CONVERGENCE!

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