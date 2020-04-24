#include <random_spanning_trees.hpp>

// modified from https://github.com/igraph/igraph/blob/master/src/spanning_trees.c
// to accommodate efficient sampling of multiple spanning trees.

#include "igraph_structural.h"
#include "igraph_dqueue.h"
#include "igraph_interface.h"
// #include "igraph_interrupt_internal.h"
#include "igraph_memory.h"
#include "igraph_adjlist.h"
#include "igraph_random.h"
#include "igraph_components.h"
#include "igraph_progress.h"
// #include "igraph_types_internal.h"

#include <assert.h>


int igraph_incident_modified(const igraph_t *graph, igraph_vector_t *eids,
                    igraph_integer_t pnode, igraph_neimode_t mode) {

    long int length = 0, idx = 0;
    long int i, j;

    long int node = pnode;

    if (node < 0 || node > igraph_vcount(graph) - 1) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVVID);
    }
    if (mode != IGRAPH_OUT && mode != IGRAPH_IN &&
        mode != IGRAPH_ALL) {
        IGRAPH_ERROR("cannot get neighbors", IGRAPH_EINVMODE);
    }

    // if (! graph->directed) {
    //     mode = IGRAPH_ALL;
    // }

    /* Calculate needed space first & allocate it*/

    if (mode & IGRAPH_OUT) {
        length += (VECTOR(graph->os)[node + 1] - VECTOR(graph->os)[node]);
    }
    if (mode & IGRAPH_IN) {
        length += (VECTOR(graph->is)[node + 1] - VECTOR(graph->is)[node]);
    }

    IGRAPH_CHECK(igraph_vector_resize(eids, length));

    if (mode & IGRAPH_OUT) {
        j = (long int) VECTOR(graph->os)[node + 1];
        for (i = (long int) VECTOR(graph->os)[node]; i < j; i++) {
            VECTOR(*eids)[idx++] = VECTOR(graph->oi)[i];
        }
    }
    if (mode & IGRAPH_IN) {
        j = (long int) VECTOR(graph->is)[node + 1];
        for (i = (long int) VECTOR(graph->is)[node]; i < j; i++) {
            VECTOR(*eids)[idx++] = VECTOR(graph->ii)[i];
        }
    }

    return 0;
}

/* igraph_random_spanning_tree */

/* Loop-erased random walk (LERW) implementation.
 * res must be an initialized vector. The edge IDs of the spanning tree
 * will be added to the end of it. res will not be cleared before doing this.
 *
 * The walk is started from vertex start. comp_size must be the size of the connected
 * component containing start.
 */

inline int RandomSpanningTrees::igraph_i_lerw_contraction_routine(igraph_vector_t *res, igraph_integer_t start, igraph_vector_bool_t *visited,
     igraph_integer_t* visited_count, igraph_vector_t* reachable_vs, igraph_vector_t* reachable_es)
{
    // now we have to check if this new start is linked to a connected contracted edge
    igraph_subcomponent(&g_contracted,  reachable_vs, start, IGRAPH_ALL);
    // assert(igraph_vector_contains(reachable_vs, start));

    const int N = igraph_vector_size(reachable_vs);


    int edge_count = 0;
    for(int i = 0; i < N ; i++)
    {
        igraph_integer_t vid = (int)VECTOR(*reachable_vs)[i];

        assert(!VECTOR(*visited)[vid]);

        // Find the directed edge, at most 1 per vertex
        igraph_incident_modified(&g_contracted, reachable_es, vid, IGRAPH_OUT);
        // assert(igraph_vector_size(&reachable_es) <= 1);
        for(int j = 0; j < igraph_vector_size(reachable_es); j++)
        {
            IGRAPH_CHECK(igraph_vector_push_back(res, VECTOR(*reachable_es)[j]));
            edge_count++;   
        }

        // igraph_vector_print(reachable_es);
        VECTOR(*visited)[vid] = 1;
        (*visited_count)++;

    }
    // igraph_vector_print(reachable_vs);
    

    // printf("N = %d, edge_count = %d\n", N, edge_count);

    assert(edge_count == N - 1);

    return IGRAPH_SUCCESS;
}

int RandomSpanningTrees::igraph_i_lerw_modified(igraph_vector_t *res, igraph_integer_t start,
                         igraph_integer_t comp_size, igraph_vector_bool_t *visited) {
    igraph_integer_t visited_count = 0;

    IGRAPH_CHECK(igraph_vector_reserve(res, igraph_vector_size(res) + comp_size - 1));

    RNG_BEGIN();

    igraph_vector_t reachable_vs, reachable_es;
    igraph_vector_init(&reachable_vs, 0);
    igraph_vector_init(&reachable_es, 0);

    // VECTOR(*visited)[start] = 1;
    // visited_count = 1;
    igraph_i_lerw_contraction_routine(res, start, visited, &visited_count, &reachable_vs, &reachable_es);


    while (visited_count < comp_size) {
        long degree, edge;
        igraph_vector_int_t *edges;

        edges = igraph_inclist_get(&il, start);

        /* choose a random edge */
        degree = igraph_vector_int_size(edges);
        edge = VECTOR(*edges)[ RNG_INTEGER(0, degree - 1) ];

        /* set 'start' to the next vertex */
        start = IGRAPH_OTHER(g, edge, start);

        /* if the next vertex hasn't been visited yet, register the edge we just traversed */
        if (! VECTOR(*visited)[start]) {
            IGRAPH_CHECK(igraph_vector_push_back(res, edge));
            // VECTOR(*visited)[start] = 1;
            // visited_count++;

            igraph_i_lerw_contraction_routine(res, start, visited, &visited_count, &reachable_vs, &reachable_es);

        }

    }
    // printf("igraph_vector_size(res) = %d\n", igraph_vector_size(res));
    assert(igraph_vector_size(res) == N - 1);

    igraph_vector_destroy(&reachable_vs);
    igraph_vector_destroy(&reachable_es);
    RNG_END();

    return IGRAPH_SUCCESS;
}


void RandomSpanningTrees::sample(igraph_vector_t *res, igraph_integer_t vid)
{
    igraph_vector_bool_t visited;
    igraph_vector_bool_init(&visited, N);

    igraph_i_lerw_modified(res, vid, N, &visited);

    igraph_vector_bool_destroy(&visited);
}

void RandomSpanningTrees::contract_edge(igraph_integer_t eid)
{
    igraph_integer_t from, to;
    igraph_edge(g, eid, &from, &to);
    igraph_add_edge(&g_contracted, from, to);

    printf("CONTRACTED edge %d \n", eid);

    // When the sampler does random walk, it should check if the newly visited vertex is connected in g_contracted
    // If yes, all the connected edges should be labelled visited
}

void RandomSpanningTrees::remove_edge(igraph_integer_t eid)
{
    igraph_integer_t from, to;
    igraph_edge(g, eid, &from, &to);
    
    // for (int i =0 ; i < il.length ; i++)
    // {
    //     igraph_vector_int_print(igraph_inclist_get(&il, i));
    // }

    // printf("\n");

    igraph_vector_int_t *edges;

    assert(from < to); // quirks of igraph_edge
    long pos;

    // remove eid from the first vertex
    edges = igraph_inclist_get(&il, from);
    // igraph_vector_int_print(edges);
    // printf("eid=%d\n", eid);
    assert(igraph_vector_int_contains(edges,eid));
    

    igraph_vector_int_search(edges, 0, eid, &pos);
    igraph_vector_int_remove(edges, pos);

    // remove eid from the second vertex
    edges = igraph_inclist_get(&il, to);
    // igraph_vector_int_print(edges);
    assert(igraph_vector_int_contains(edges,eid));
    igraph_vector_int_search(edges, 0, eid, &pos);
    igraph_vector_int_remove(edges, pos);


    // printf("REMOVED edge %d \n", eid);

    // for (int i =0 ; i < il.length ; i++)
    // {
    //     igraph_vector_int_print(igraph_inclist_get(&il, i));
    // }

    // printf("\n");

    // exit(-1);

    // After this operation, the edge is still in the graph, but the algorithm will not random walk from this edge anymore
}