/*
 * MatchingMachine.h
 *
 */
/*
Copyright (c) 2023

This library contains portions of other open source products covered by separate
licenses. Please see the corresponding source files for specific terms.

ArcMatch is provided under the terms of The MIT License (MIT):

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#pragma once

#ifndef MATCHINGMACHINE_H_
#define MATCHINGMACHINE_H_

#include "Graph.h"

namespace rilib
{
class FMatchingMachineEdge
{
public:
    int Source;
    int Target;
    void* EdgeAttribute;
    int Id;

    FMatchingMachineEdge(int _source, int _target, void* _attr, int _id)
    {
        Source = _source;
        Target = _target;
        EdgeAttribute = _attr;
        Id = _id;
    }

    FMatchingMachineEdge()
    {
        Source = -1;
        Target = -1;
        EdgeAttribute = NULL;
        Id = -1;
    }
};

enum EParentType { PARENTTYPE_IN, PARENTTYPE_OUT, PARENTTYPE_NULL };

class FRiMatchingMachine
{
public:
    int NumOfQueryVertex;

    void** NodesAttributes; // indexed by state_id
    int* EdgeSizes; // indexed by state_id
    int* OutEdgeSizes; // indexed by state_id
    int* InEdgeSizes; // indexed by state_id
    FMatchingMachineEdge** Edges; // indexed by state_id, map on states  (0,1) = (state0, state1)

    int* QueryVertexToState; // indexed by node_id
    int* StateToQueryVertex; // indexed by state_id

    int* ParentState; // indexed by state_id
    EParentType* ParentType; // indexed by state id

    int nof_leafs;

    FRiMatchingMachine(FRiGraph& query)
    {
#ifdef MDEBUG
        std::cout << "mama constructor (" << query.NumOfVertex << ")...\n";
#endif
        NumOfQueryVertex = query.NumOfVertex;
        NodesAttributes = new void*[NumOfQueryVertex];
        EdgeSizes = (int*)calloc(NumOfQueryVertex, sizeof(int));
        OutEdgeSizes = (int*)calloc(NumOfQueryVertex, sizeof(int));
        InEdgeSizes = (int*)calloc(NumOfQueryVertex, sizeof(int));
        Edges = new FMatchingMachineEdge*[NumOfQueryVertex];

        QueryVertexToState = (int*)calloc(NumOfQueryVertex, sizeof(int));
        StateToQueryVertex = (int*)calloc(NumOfQueryVertex, sizeof(int));
        ParentState = (int*)calloc(NumOfQueryVertex, sizeof(int));
        ParentType = new EParentType[NumOfQueryVertex];

        nof_leafs = 0; // only used by MaMaxxxLeafs
#ifdef MDEBUG
        std::cout << "done...\n";
#endif
    }

    virtual ~FRiMatchingMachine()
    {
        delete[] NodesAttributes;
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            delete[] Edges[i];
        }
        delete[] Edges;
        free(EdgeSizes);
        free(OutEdgeSizes);
        free(InEdgeSizes);
        free(QueryVertexToState);
        free(StateToQueryVertex);
        free(ParentState);
        delete[] ParentType;
    }

    void fix_eids(FRiGraph& query)
    {
        int source, target, eid;
        for (int si = 0; si < NumOfQueryVertex; si++)
        {
            for (int ei = 0; ei < EdgeSizes[si]; ei++)
            {
                source = StateToQueryVertex[Edges[si][ei].Source];
                target = StateToQueryVertex[Edges[si][ei].Target];

                eid = 0;
                for (int i = 0; i < source; i++)
                {
                    eid += query.OutAdjSizes[i];
                }
                for (int i = 0; i < query.OutAdjSizes[source]; i++)
                {
                    if (query.OutAdjList[source][i] == target)
                    {
                        Edges[si][ei].Id = eid;
                        break;
                    }
                    eid++;
                }
            }
        }
    }

    void print()
    {
        std::cout << "| MatchingMachine:  nof sates " << NumOfQueryVertex << "\n";
        std::cout << "| \tmap state_to_node(";
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            std::cout << "[" << i << ":" << StateToQueryVertex[i] << "]";
        }
        std::cout << ")\n";
        std::cout << "| \tmap node_to_state(";
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            std::cout << "[" << i << ":" << QueryVertexToState[i] << "]";
        }
        std::cout << ")\n";
        std::cout << "| \tstates (node)(parent state, parent type)\n";
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            std::cout << "| \t\t[" << i << "] (" << StateToQueryVertex[i] << ") (" << ParentState[i] << ", ";
            switch (ParentType[i])
            {
                case PARENTTYPE_IN:
                    std::cout << "IN";
                    break;
                case PARENTTYPE_OUT:
                    std::cout << "OUT";
                    break;
                case PARENTTYPE_NULL:
                    std::cout << "NULL";
                    break;
            }
            std::cout << ")\n";
            std::cout << "| \t\t\tchecking[" << EdgeSizes[i] << "] ";
            for (int j = 0; j < EdgeSizes[i]; j++)
            {
                std::cout << "{s(" << Edges[i][j].Id << ")(" << Edges[i][j].Source << "," << Edges[i][j].Target << "):";
                std::cout << "n(" << StateToQueryVertex[Edges[i][j].Source] << "," << StateToQueryVertex[Edges[i][j].Target] << ")}";
            }
            std::cout << "\n";
        }
    }

public:
    virtual void Build(FRiGraph& ssg) = 0;
};
} // namespace rilib

#endif /* MATCHINGMACHINE_H_ */
