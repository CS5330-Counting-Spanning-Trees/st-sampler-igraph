#pragma once

#include <igraph/igraph.h>

#include <vector>

#include <cassert>

#include <algorithm>

#include <random>

#define eid_t int
#define vid_t int

// Simple and Memory Efficient Undirected Graph, to Allow Contraction and Edge Removal
class GraphLite{
public:
    typedef struct edge{
        vid_t from;
        vid_t to;

        edge(vid_t from, vid_t to) : from(from), to(to){}
    }edge_t;

    inline GraphLite(const igraph_t* g);

    void clear(){edge_list_.clear(); inclist_.clear();e_removed_count = 0; v_removed_count = 0;}

    // vertices and edges could be marked removed, hence requires more care when counting
    vid_t vertex_count(){return inclist_.size() - v_removed_count;}
    vid_t vertex_count_all(){return inclist_.size();}
    eid_t edge_count(){return edge_list_.size() - e_removed_count;}
    eid_t edge_count_all(){return edge_list_.size();}

    // Obtain element
    const auto& edge(eid_t e){return edge_list_[e];}
    inline vid_t edge_other_end(eid_t e, vid_t v1);
    inline vid_t first_connected_vertex();
    inline vid_t random_connected_vertex();
    void invalidate_edge(eid_t e){edge_list_[e].from =  edge_list_[e].to = -1;}
    bool is_edge_valid(eid_t e){return edge_list_[e].from >= 0 && edge_list_[e].to >= 0;}

    // Obtain whole structures
    const auto& edge_list(){return edge_list_;}
    const auto& inclist(){return inclist_;}

    // Modify the Graph
    inline const std::vector<eid_t> contract_edge(eid_t e, bool return_other_edges_removed = true); // return number of contracted edges other than the current one
    inline void remove_edge(eid_t e);

    inline void sanity_check(bool check_connected = false);
    inline void print();

private:
    
    // Memory Complexity = 2M * vid_t + 2M * eid_t + N * std::vector
    std::vector<edge_t> edge_list_;
    std::vector<std::vector<eid_t>> inclist_; // incident edges for each vertex
    eid_t e_removed_count;
    vid_t v_removed_count;


};

GraphLite::GraphLite(const igraph_t* g)
{
    clear();

    printf("GraphLite: Building Graph from igraph structure...\n");
    inclist_.resize(igraph_vcount(g));
    edge_list_.reserve(igraph_ecount(g));

    printf("V = %d, E = %d\n", igraph_vcount(g), igraph_ecount(g));


    const eid_t M = igraph_ecount(g);

    for(eid_t i = 0; i < M ; i++){

        // add edge
        auto& e = edge_list_.emplace_back((vid_t)VECTOR(g->from)[i],(vid_t)VECTOR(g->to)[i]);

        // add incident edge to the two vertices
        inclist_[e.from].push_back(edge_list_.size()-1);
        inclist_[e.to].push_back(edge_list_.size()-1);
    }

    sanity_check(true);

    printf("Done\n");
}

inline vid_t GraphLite::edge_other_end(eid_t e, vid_t v1)
{
    if (edge_list_[e].from == v1)
        return edge_list_[e].to;
    else if (edge_list_[e].to == v1)
        return edge_list_[e].from;
    else{
        printf("Wrong edge %d for vertex %d\n", e, v1);
        print();
        abort();
    }
        
}

vid_t GraphLite::first_connected_vertex()
{
    vid_t v = 0;

    while(!inclist_[v].size()){
        v++;
        assert(v < vertex_count_all());
    }
    return v;
}

vid_t GraphLite::random_connected_vertex()
{
    vid_t v;

    std::mt19937 mt_rand(time(0));

    while(true){
        v = mt_rand() %  vertex_count_all();
        if (inclist_[v].size())
            break;
    }
    return v;
}


const std::vector<eid_t> GraphLite::contract_edge(eid_t e_in, bool return_other_edges_removed)
{
    std::vector<eid_t> ret;
    assert(e_in < edge_count_all());

    vid_t from = edge_list_[e_in].from;
    vid_t to = edge_list_[e_in].to;

    //// decide which vertex to remove

    if (inclist_[from].size() > inclist_[to].size() )
        std::swap(from,to);

    //// redirect all edges, remove all edges connect between the two

    vid_t edge_removed = 0;

    for (size_t i = 0; i < inclist_[to].size();){
        
        eid_t e = inclist_[to][i];

        // assert no self loops, or wrong edge pointed
        assert( (edge_list_[e].from == to) != (edge_list_[e].to == to) );

        if(edge_list_[e].from == from || edge_list_[e].to == from){
            remove_edge(e);
            edge_removed++;
            if(return_other_edges_removed && e != e_in)
                ret.push_back(e);
            continue;
        }

        if(edge_list_[e].from == to){
            edge_list_[e].from = from;
            inclist_[from].push_back(e);
        }
        else if (edge_list_[e].to == to){
            edge_list_[e].to = from;
            inclist_[from].push_back(e);
        }
        else{
            printf("Impossible case\n");
            abort();
        }
        i++;
    }

    //// remove one of the vertex
    inclist_[to].clear();
    v_removed_count++;

    printf("Contracted edge %d and removed vertex %d and additional %d edges \n", e_in, to, edge_removed-1);
    return ret;
}

// remove edge only remove from the inclist, not the edge_list, to conserve the edge id
void GraphLite::remove_edge(eid_t e)
{
    assert(e < edge_count_all());

    vid_t from = edge_list_[e].from;
    vid_t to = edge_list_[e].to;

    // invalidate
    invalidate_edge(e);


    auto iter = std::find(inclist_[from].begin(),inclist_[from].end(), e);
    assert(iter != inclist_[from].end());
    inclist_[from].erase(iter);

    iter = std::find(inclist_[to].begin(),inclist_[to].end(), e);
    assert(iter != inclist_[to].end());
    inclist_[to].erase(iter);

    e_removed_count++;
    printf("Removed edge %d\n", e);
}


void GraphLite::sanity_check(bool check_connected)
{

    eid_t ecount_dir = 0;
    vid_t v_removed = 0;
    for(auto& e : inclist_){
        if(check_connected)
            assert(e.size()); // should always have incident edge(s)
        else if(!e.size())
            v_removed++;

        ecount_dir += e.size();
    }

    // must have pairs of incident entries
    assert(ecount_dir % 2 == 0);

    // check consistency of edge count in both lists
    assert(ecount_dir / 2 == edge_count());

    assert(v_removed == v_removed_count);
}

void GraphLite::print()
{
    printf("incident list:\n");
    for(size_t i = 0; i < inclist_.size(); i++){
        printf("%3d: ", (int)i);
        for(auto e : inclist_[i])
            printf("%3d ", e);
        printf("\n");
    }

    printf("edge list:\n");

    for(size_t i = 0; i < edge_list_.size(); i++)
        printf("%3d: %3d, %3d\n", int(i), edge_list_[i].from, edge_list_[i].to);
    printf("\n");
}