#include <approx_count_st.hpp>
#include <random_spanning_trees.hpp>

#include <random>
#include <algorithm>

ApproxCountST::ApproxCountST(GraphLite* gl) : gl(gl), N_initial(gl->vertex_count_all()), M_initial(gl->edge_count_all())
{
    assert(M_initial >= N_initial-1);
    // K = M_initial - (N-1);
    K = M_initial; // the last round should always yield ratio of 1, for presence
    N = N_initial;

    ps_vec.resize(K);

    printf("Approximate Count ST initialised with a graph of %d vertices and %d edges\n", N, M_initial);
    printf("Iterations of ratio estimators to run: %d rounds\n", K);
}



ApproxCountST::result_t ApproxCountST::approx_count_st()
{
    printf("approx_count_st...\n");
    // Obtain a shuffled edge sequence, or just iterate as per original sequence
    bool do_shuffle_edges = true;
    if(do_shuffle_edges)
    {
        std::vector<int> e(M_initial);
        for(int i = 0 ; i < M_initial; i++)
            e[i] = i;
        
        std::random_device rd;
        std::mt19937 gen(rd());

        std::shuffle(e.begin(), e.end(), gen);

        for(int i = 0 ; i < M_initial; i++)
            ps_vec[i].eid = e[i];

        // so far all ps_vec element's mode should be UNSPECIFIED
    }else{
        for(int i = 0 ; i < K; i++)
            ps_vec[i].eid = i;
    }

    

    printf("ps_vec[i].eid done...\n");
    
    // Initialise Random Spanning Tree Sampler
    RandomSpanningTrees rst(gl);

    printf("rst initialised...\n");

    // for each k-th loop, we keep drawing batch samples for the k-th ratio, until convergence
    // This loop will calculate ratio of G_M / G_{M-1} ..... G_{2} / G_{1} * G{1}, total of M-1 - (N-1) + 1 terms
    for(int k = 0; ; )
    {
        const auto& e = gl->edge(ps_vec[k].eid);

        // TODO: make it more streamlined
        // if the edge is invalid, means some other present edge has contracted this one, so the ratio automatically should be 1
        if(!gl->is_edge_valid(ps_vec[k].eid)){
            assert(ps_vec[k].ratio == 1.0);
            printf("short circuiting edge %d, as the edge has been contracted by others before...\n", ps_vec[k].eid);
            if (++k == K)
                break;
        }
        
        // Logging the ripple sample count if required
        if(!ps_vec[k].total)
        {
            printf("\nEdge %d (%d->%d)\n",ps_vec[k].eid, e.from, e.to);
            ps_vec[k].rippled_total = -1;
        }
        else if(!ps_vec[k].rippled_total)
        {
            printf("\n[%.1lf%%] Edge %d (%d->%d) with existing %d rippled samples, mode %d\n ", 
            100.0 * (k+1)/ K, ps_vec[k].eid, e.from, e.to,  ps_vec[k].total, ps_vec[k].count_mode);

            ps_vec[k].rippled_total = ps_vec[k].total;
            // Good! Enough samples were obtained to output past stats
            if(ps_vec[k].update_count >= PIVOT_BUFFER_SIZE){
                printf("Past ratio buffer: %.3lf %.3lf %.3lf %.3lf\n", ps_vec[k].ratio_buffer[0], ps_vec[k].ratio_buffer[1] ,ps_vec[k].ratio_buffer[2] ,ps_vec[k].ratio_buffer[3]);
            }
        }
    
        sample_mini_batch_with_updates(&rst, k); // affected by random_walk_mode
        if (ps_vec[k].converged()){

            // make graph changes
            switch(ps_vec[k].count_mode){
                case PRESENCE:
                    // to contract the edge
                    gl->contract_edge(ps_vec[k].eid);
                    break;
                case ABSENCE:
                    // to delete the edge in the incident list
                    gl->remove_edge(ps_vec[k].eid);
                    break;
                default:
                    assert(0); // should never come here
            }

            printf("%d-th of %d ratio converged to %.3lf\n", k, K, ps_vec[k].ratio);

            if (++k == K)
                break;

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
    // TODO: put this outside the function to make initialisation faster
    std::vector<eid_t> path;
    std::vector<eid_t> next(N);
    std::vector<bool> in_tree(N);

    const vid_t root = gl->first_connected_vertex();
    // const vid_t root = gl->random_connected_vertex();

    // perform the batch sampling
    for (int i = 0 ; i < BATCH_SIZE ; i++)
	{
        rst->wilsons_get_st(&path, root, &next, &in_tree); // NOTE: for now, always sample from the node 0

        // NOTE: change K to k_start + 1, to disable ripple feature
        for(int k = k_start; k < K ;k++)
        {
            if ( path.end() != std::find(path.begin(), path.end(), ps_vec[k].eid) ){
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