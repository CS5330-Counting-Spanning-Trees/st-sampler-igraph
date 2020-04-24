#include <approx_count_st.hpp>
#include <random_spanning_trees.hpp>

#include <random>
#include <algorithm>

ApproxCountST::ApproxCountST(igraph_t* g, random_walk_mode_t rw_mode) : random_walk_mode(rw_mode), g(g), N(igraph_vcount(g)), M(igraph_ecount(g))
{
    assert(M >= N-1);
    // K = M - (N-1);
    K = M; // the last round should always yield ratio of 1, for presence

    ps_vec.resize(K);

    printf("Approximate Count ST initialised with a graph of %d vertices and %d edges\n", N, M);
    printf("Iterations of ratio estimators to run: %d rounds\n", K);
}



ApproxCountST::result_t ApproxCountST::approx_count_st()
{
    printf("approx_count_st...\n");
    // Obtain a shuffled edge sequence, or just iterate as per original sequence
    // {
    //     std::vector<igraph_integer_t> e(M);
    //     for(int i = 0 ; i < N; i++)
    //         e[i] = i;
        
    //     std::random_device rd;
    //     std::mt19937 g(rd());

    //     std::shuffle(e.begin(), e.end(), g);

    //     for(int i = 0 ; i < M; i++)
    //         ps_vec[i].eid = e[i];

    //     // so far all ps_vec element's mode should be UNSPECIFIED
    // }

    for(int i = 0 ; i < K; i++)
        ps_vec[i].eid = i;

    printf("ps_vec[i].eid done...\n");
    
    // Initialise Random Spanning Tree Sampler
    RandomSpanningTrees rst(g);

    printf("rst initialised...\n");

    // for each k-th loop, we keep drawing batch samples for the k-th ratio, until convergence
    // This loop will calculate ratio of G_M / G_{M-1} ..... G_{N} / G_{N-1}, total of M-1 - (N-1) + 1 terms
    for(int k = 0; ; )
    {
        sample_mini_batch_with_updates(&rst, k); // affected by random_walk_mode
        if (ps_vec[k].converged()){

            // make graph changes
            switch(ps_vec[k].count_mode){
                case PRESENCE:
                    // to contract the edge
                    rst.contract_edge(ps_vec[k].eid);
                    break;
                case ABSENCE:
                    // to delete the edge in the incident list
                    rst.remove_edge(ps_vec[k].eid);
                    break;
                default:
                    assert(0); // should never come here
            }

            printf("%d-th of %d ratio converged to %.3lf\n", k, K, ps_vec[k].ratio);

            if (++k == K)
                break;

            printf("  moving to new ratio %.3lf with existing %d rippled samples \n ", ps_vec[k].ratio, ps_vec[k].total);
            ps_vec[k].rippled_total = ps_vec[k].total;
        }

        
    }

    // Prepare final result
    assert( ps_vec[K-1].ratio == 1.0);

    result_t res;
    res.count = 1.0;

    for (int k = 0; k < K; k++)
    {
        res.count *= 1/ps_vec[k].ratio;
        res.effective_samples += ps_vec[k].total;
        res.actual_samples += ps_vec[k].total - ps_vec[k].rippled_total;
    }

    
    return res;
}

// Draw new samples starting from index k
void ApproxCountST::sample_mini_batch_with_updates(RandomSpanningTrees* rst, int k_start)
{
    
    const int BATCH_SIZE = ps_vec[k_start].requested_batch_size;
    // printf("mini_batch at [%d] for %d samples\n", k_start, BATCH_SIZE);

    // store the sample obtained from the sampler
    igraph_vector_t res;
	igraph_vector_init(&res,0);

    // perform the batch sampling
    for (int i = 0 ; i < BATCH_SIZE ; i++)
	{
        igraph_vector_clear(&res);
		rst->sample(&res, 0);

        // NOTE: change K to k_start + 1, to disable ripple feature
        for(int k = k_start; k < K ;k++)
        {
            if (igraph_vector_contains(&res, ps_vec[k].eid)){
                ps_vec[k].present++;
                if(ps_vec[k].count_mode != PRESENCE)
                    break; // We should not continue to update the downstreams in this case
            }else{
                ps_vec[k].absent++;
                if(ps_vec[k].count_mode != ABSENCE)
                    break; // We should not continue to update the downstreams in this case
            }
        }
	}


    // ripple update
    for(int k = k_start; k < K ;k++)
    {
        // VECTOR(g->from)[(long int)ps_vec[k].eid] = 0; // for undirected, from always greater than to
        // igraph_integer_t eid = ps_vec[k].eid;
        
        ps_vec[k].update();

        // printf("%d -> ", k);
        // ps_vec[k].print();

        // check if the count mode can be specified now
        if(ps_vec[k].count_mode == UNSPECIFIED)
        {
            // we must ensure all upstream chains has been set properly before we update
            if (k == 0 || ps_vec[k-1].count_mode != UNSPECIFIED)
                ps_vec[k].try_set_count_mode();
            break; // no point update downstream ratios
        }
        
    }
    
}

void ApproxCountST::print_all()
{
    for(int i=0;i<K;i++)
    {
        printf("%d: %.3lf(%d)\t", i, 1/ps_vec[i].ratio, ps_vec[i].total);
    }

    printf("\n");
}