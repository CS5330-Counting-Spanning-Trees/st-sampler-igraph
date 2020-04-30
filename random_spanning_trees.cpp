#include <random_spanning_trees.hpp>
#include <assert.h>

int RandomSpanningTrees::wilsons_get_st(std::vector<eid_t> *path, vid_t root, std::vector<eid_t>* next, std::vector<bool>* in_tree)
{   

    const vid_t N = gl->vertex_count_all(); // The count may change every time

    assert(gl->inclist()[root].size()); // assert reachability of the root

    // in_tree should be initialised to 0
    std::fill(in_tree->begin(), in_tree->end(), false);
    (*in_tree)[root] = true;

    // next vector need not to be cleared, as it will be overwritted properly. it should have size N
    
    path->clear(); // capacity unchanged
    path->reserve(gl->vertex_count()-1);

    RNG_BEGIN();

    for (vid_t i = 0; i < N; i++)
    {
        if (!gl->inclist()[i].size())
            continue;
        
        vid_t u = i;

        while(!(*in_tree)[u])
        {
            // generate random successor
            auto& edges = gl->inclist()[u]; // O(1)
            eid_t edge = edges[ RNG_INTEGER(0, edges.size() - 1) ];
            
            (*next)[u] = edge;
            u = gl->edge_other_end(edge, u); // O(1)
            
        }

        // collecting the walked path
        u = i;
        while(!(*in_tree)[u])
        {
            (*in_tree)[u] = true;

            eid_t edge = (*next)[u];
            path->push_back(edge);
            u = gl->edge_other_end(edge, u); // O(1)
        }
    }

    RNG_END();


    return IGRAPH_SUCCESS;
}
