#pragma once


#ifndef ARC_MATCH_CONTROL
    # define ARC_MATCH_CONTROL 14

    #if ARC_MATCH_CONTROL == 1
    	#define MAMA_1
    	#define SOLVER_0

    #elif ARC_MATCH_CONTROL == 2
    	#define MAMA_1
    	#define SOLVER_ED

    #elif ARC_MATCH_CONTROL == 3
    	#define MAMA_1
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 4
    	#define MAMA_NS
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 5
    	#define REDUCE_EDGES
    	#define PATH_LENGTH 6
    	#define MAMA_NS
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 6
    	#define MAMA_NSL
    	#define SOLVER_LF

    #elif ARC_MATCH_CONTROL == 7
    	#define REDUCE_EDGES
    	#define PATH_LENGTH 6
    	#define MAMA_NSL
    	#define SOLVER_LF

    #elif ARC_MATCH_CONTROL == 8
    	#define NODE_D_CONV
    	#define MAMA_1
    	#define SOLVER_0

    #elif ARC_MATCH_CONTROL == 9
    	#define NODE_D_CONV
    	#define MAMA_1
    	#define SOLVER_ED

    #elif ARC_MATCH_CONTROL == 10
    	#define NODE_D_CONV
    	#define MAMA_1
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 11
    	#define NODE_D_CONV
    	#define MAMA_NS
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 12
    	#define NODE_D_CONV
    	#define REDUCE_EDGES
    	#define PATH_LENGTH 6
    	#define MAMA_1
    	#define SOLVER_0

    #elif ARC_MATCH_CONTROL == 13
    	#define NODE_D_CONV
    	#define MAMA_NSL
    	#define SOLVER_LF

    #elif ARC_MATCH_CONTROL == 14
    	#define NODE_D_CONV
    	#define REDUCE_EDGES
    	#define PATH_LENGTH 6
    	#define MAMA_NSL
    	#define SOLVER_LF

    #elif ARC_MATCH_CONTROL == 15
    	#define NODE_D_CONV
    	#define EDGE_D_CONV
    	#define MAMA_NS
    	#define SOLVER_DP

    #elif ARC_MATCH_CONTROL == 16
    	#define NODE_D_CONV
    	#define EDGE_D_CONV
    	#define MAMA_NSL
    	#define SOLVER_LF

    #endif
#endif




