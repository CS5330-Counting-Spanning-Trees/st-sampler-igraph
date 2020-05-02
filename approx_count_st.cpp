#include <approx_count_st.hpp>
#include <random_spanning_trees.hpp>

#include <random>
#include <algorithm>

#include "graph_lite.hpp"

ApproxCountST::ApproxCountST(GraphLite* gl) : ps_vec(gl->edge_count_all()), gl(gl), N_initial(gl->vertex_count_all()), M_initial(gl->edge_count_all()), K(M_initial)
{
    assert(M_initial >= N_initial-1);

    printf("Approximate Count ST initialised with a graph of %d vertices and %d edges\n", N_initial, M_initial);
    printf("Iterations of ratio estimators to run: %d rounds\n", K);
}



ApproxCountST::result_t ApproxCountST::approx_count_st()
{
    printf("approx_count_st...\n");
    // Obtain a shuffled edge sequence, or just iterate as per original sequence

    e_shuffle.resize(M_initial);
    for(int i = 0 ; i < M_initial; i++)
        e_shuffle[i] = i;

    const bool do_shuffle_edges = false;
    if(do_shuffle_edges)
    {
        std::random_device rd;
        std::mt19937 gen(rd());

        std::shuffle(e_shuffle.begin(), e_shuffle.end(), gen);

        // so far all ps_vec element's mode should be UNSPECIFIED
    }

    for(int i = 0 ; i < M_initial; i++)
        ps_vec[e_shuffle[i]].eid = i;

    printf("ps_vec[i].eid done...\n");
    
    // Initialise Random Spanning Tree Sampler
    RandomSpanningTrees rst(gl);

    printf("rst initialised...\n");

    // for each k-th loop, we keep drawing batch samples for the k-th ratio, until convergence
    // This loop will calculate ratio of G_M / G_{M-1} ..... G_{2} / G_{1} * G{1}, total of M-1 - (N-1) + 1 terms

    // store the sample obtained from the sampler
        
    sampling_struct_t sampling_struct;
    sampling_struct.next.resize(N_initial);
    sampling_struct.in_tree.resize(N_initial);

    for(int k = 0; k < K ; )
    {
        // TODO: make it more streamlined
        // if the edge is invalid, means some other present edge has contracted this one, so the ratio automatically should be 1
        if(!gl->is_edge_valid(ps_vec[k].eid)){
            
            printf("short circuiting edge %d, as the edge has been contracted by others before...\n", ps_vec[k].eid);
            
            // possible to have ratio = -1, uninitialised, as we contract it before it got even attempted.
            if (ps_vec[k].ratio == -1.0)
                ps_vec[k].ratio = 1.0;
            assert(ps_vec[k].ratio == 1.0);
            k++;
            continue;
        }
        
        const auto e = gl->edge(ps_vec[k].eid);

        //// Logging the ripple sample count if required
        if(!ps_vec[k].total)
        {
            printf("\nEdge %d (%d->%d)\n",ps_vec[k].eid, e.from, e.to);
            ps_vec[k].rippled_total = -1;

            if(k != 0){
                ps_vec[k].print();
                abort();
            }
        }
        else if(!ps_vec[k].rippled_total)
        {
            printf("\n[%.1lf%%] Edge %d (%d->%d) with existing %d rippled samples (count = %d), mode %d\n ", 
            100.0 * (k+1)/ K, ps_vec[k].eid, e.from, e.to,  ps_vec[k].total, ps_vec[k].update_count, ps_vec[k].count_mode);

            ps_vec[k].rippled_total = ps_vec[k].total;
            // Good! Enough samples were obtained to output past stats
            if(ps_vec[k].update_count >= PIVOT_BUFFER_SIZE){
                printf("Past ratio buffer:");
                for(auto e : ps_vec[k].inverse_ratio_buffer)
                    printf("%.3lf ", e);
                printf("\n");
            }
        }
        //// End logging

        // k will increment if convergence is checked
        // This check is before the first ever sampling is done
        if (check_convergence(&k)){
            k++;
            continue;
        }
    
        sample_mini_batch_with_updates(&rst, k, &sampling_struct); // affected by random_walk_mode
        printf(".");
        fflush(stdout);
    }

    // Prepare final result
    assert( ps_vec[K-1].ratio == 1.0);

    result_t res;
    res.count = 1.0;

    for (int k = 0; k < K; k++)
    {
        assert(ps_vec[k].ratio > 0.1);
        res.count *= 1/ps_vec[k].ratio;
        res.count_log += std::log(1/ps_vec[k].ratio);
        res.effective_samples += ps_vec[k].total;
        res.actual_samples += ps_vec[k].total - ps_vec[k].rippled_total;
    }
    printf("approx_count_st COMPLETED...\n");
    return res;
}

inline bool ApproxCountST::check_convergence(eid_t* pk)
{
    eid_t& k = *pk;
    assert(k >= 0 && k < K);
    if (ps_vec[k].converged()){

        printf("%d-th of %d ratio converged to %.3lf\n", k + 1, K, ps_vec[k].ratio);

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

        return true;
    }

    return false;
}

// Draw new samples starting from index k
void ApproxCountST::sample_mini_batch_with_updates(RandomSpanningTrees* rst, int k_start, sampling_struct_t* sampling_struct)
{
    assert(k_start >=0 && k_start < K);
    const int BATCH_SIZE = ps_vec[k_start].requested_batch_size;
    // printf("mini_batch at [%d] for %d samples\n", k_start, BATCH_SIZE);

    // const vid_t root = gl->first_connected_vertex();
    const vid_t root = gl->random_connected_vertex();

    // perform the batch sampling
    for (int i = 0 ; i < BATCH_SIZE ; i++)
	{
        rst->wilsons_get_st(&(sampling_struct->path), root, &(sampling_struct->next), &(sampling_struct->in_tree)); // NOTE: for now, always sample from the node 0

        // NOTE: change K to k_start + 1, to disable ripple feature
        for(int k = k_start; k < K ;k++)
        {
            // If the edge is contracted away by edges before it, no need to update further!
            if (!gl->is_edge_valid(ps_vec[k].eid))
            {
                assert(k!=k_start);
                continue;
            }
            if ( sampling_struct->path.end() != std::find(sampling_struct->path.begin(), sampling_struct->path.end(), ps_vec[k].eid) ){
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

        // If the edge is contracted away by edges before it, no need to update further!
        if (!gl->is_edge_valid(ps_vec[k].eid))
            continue;
        
        ps_vec[k].update();

        if(ps_vec[k].count_mode == UNSPECIFIED && !ps_vec[k].try_set_count_mode()){
            break;
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