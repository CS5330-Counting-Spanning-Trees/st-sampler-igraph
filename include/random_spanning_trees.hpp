#pragma once

#include <igraph/igraph.h>


class RandomSpanningTrees{

public:
    RandomSpanningTrees(igraph_t* g) : g(g), N(igraph_vcount(g)), M(igraph_ecount(g)){
        igraph_inclist_init(g, &il, IGRAPH_ALL); // Complexity O(|V|+|E|)
        igraph_empty(&g_contracted, N, IGRAPH_UNDIRECTED);
    };

    // obtain a random spanning tree sample
    void sample(igraph_vector_t *res, igraph_integer_t vid);

    void contract_edge(igraph_integer_t eid);

    // This is achieved by manipulating the inclist
    void remove_edge(igraph_integer_t eid);

    ~RandomSpanningTrees(){
        igraph_inclist_destroy(&il);
    }


private:
    igraph_t* g = nullptr;
    igraph_t g_contracted; // Store all the contracted edges, using the same vertices ordering as g
    igraph_inclist_t il;

    const igraph_integer_t N;
	const igraph_integer_t M;
    // igraph_integer_t component_size;

    inline int igraph_i_lerw_contraction_routine(igraph_vector_t *res, igraph_integer_t start, igraph_vector_bool_t *visited,
        igraph_integer_t* visited_count, igraph_vector_t* reachable_vs, igraph_vector_t* reachable_es);

    int igraph_i_lerw_modified(igraph_vector_t *res, igraph_integer_t start,
                         igraph_integer_t comp_size, igraph_vector_bool_t *visited);
};

// int igraph_random_spanning_trees(const igraph_t *graph, const igraph_inclist_t* il, igraph_vector_t *res, igraph_integer_t vid);

