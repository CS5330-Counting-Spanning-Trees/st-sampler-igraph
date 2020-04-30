#pragma once

#include "graph_lite.hpp"

class RandomSpanningTrees{

public:



    RandomSpanningTrees(GraphLite* gl) : gl(gl){}

    
    // Wilson's Algorithm Implementation
    int wilsons_get_st(std::vector<eid_t> *path, vid_t root, std::vector<eid_t>* next, std::vector<bool>* in_tree);

private:

    GraphLite* gl = nullptr;
    
};

