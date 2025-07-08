/*
 * MaMaConstrFirst.h
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

#ifndef MAMACONSTRFIRSTNSCC_H_
#define MAMACONSTRFIRSTNSCC_H_

#include "AttributeComparator.h"
#include "Graph.h"
#include "MatchingMachine.h"
#include "sbitset.h"
#include <algorithm>
#include <set>
#include <vector>

namespace rilib
{
class FAmMaMaConstrFirstNSCC : public FRiMatchingMachine
{
    enum NodeFlag { NS_CORE, NS_CNEIGH, NS_UNV };

    FAmsbitset* domains;
    int* domains_size;

    FAmAttributeComparator& nodeComparator;
    FAmAttributeComparator& edgeComparator;

public:
    FAmMaMaConstrFirstNSCC(FRiGraph& query, FAmsbitset* _domains, int* _domains_size, FAmAttributeComparator& _nodeComparator, FAmAttributeComparator& _edgeComparator)
        : FRiMatchingMachine(query)
        , domains(_domains)
        , domains_size(_domains_size)
        , nodeComparator(_nodeComparator)
        , edgeComparator(_edgeComparator)
    {
    }

    void Build(FRiGraph& ssg) override
    {
        auto node_flags = new NodeFlag[NumOfQueryVertex]; // indexed by node_id
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            node_flags[i] = NS_UNV;

            // used for recognizing core compatibility parents
            ParentState[i] = -1;
        }

        int si = 0;
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            if (domains_size[i] == 1)
            {
                std::cout << "ssi[" << si << "] = " << i << "\n";
                push_node_to_core(i, si, node_flags, ssg, StateToQueryVertex, QueryVertexToState);
                si++;
            }
        }

        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            std::cout << i << "[" << node_flags[i] << "] ";
        }
        std::cout << "\n";

        for (; si < NumOfQueryVertex; si++)
        {
            std::cout << "SI[" << si << "]\n";
            int best_nid = -1;
            int best_nid_score[] = {0, 0, 0, 0, 0};
            int current_nid_score[] = {0, 0, 0, 0, 0};

            for (int nid = 0; nid < NumOfQueryVertex; nid++)
            {
                std::cout << nid << "(" << node_flags[nid] << ") ";
                get_scores(nid, current_nid_score, node_flags, ssg);
                std::cout << "[";
                for (int i = 0; i < 5; i++)
                {
                    std::cout << current_nid_score[i] << " ";
                }
                std::cout << "]\n";
            }

            for (int nid = 0; nid < NumOfQueryVertex; nid++)
            {
                if (node_flags[nid] == NS_CNEIGH)
                {
                    if (best_nid == -1)
                    {
                        best_nid = nid;
                        get_scores(nid, best_nid_score, node_flags, ssg);
                    }
                    else
                    {
                        get_scores(nid, current_nid_score, node_flags, ssg);
                        if (compare_scores(current_nid_score, best_nid_score, nid, best_nid) > 0)
                        {
                            best_nid = nid;
                            for (int i = 0; i < 5; i++)
                            {
                                best_nid_score[i] = current_nid_score[i];
                            }
                        }
                    }
                }
            }

            if (best_nid == -1)
            {
                // firs node without singletons or disconnected query
                for (int nid = 0; nid < NumOfQueryVertex; nid++)
                {
                    if (node_flags[nid] == NS_UNV)
                    {
                        if (best_nid == -1)
                        {
                            best_nid = nid;
                            get_scores(nid, best_nid_score, node_flags, ssg);
                        }
                        else
                        {
                            get_scores(nid, current_nid_score, node_flags, ssg);
                            if (compare_scores(current_nid_score, best_nid_score, nid, best_nid) > 0)
                            {
                                best_nid = nid;
                                for (int i = 0; i < 5; i++)
                                {
                                    best_nid_score[i] = current_nid_score[i];
                                }
                            }
                        }
                    }
                }
            }
            std::cout << "si[" << si << "] = " << best_nid << "\n";

            std::set<int> cascade;
            if (best_nid_score[0] > 0)
            {
                for (int i = 0; i < NumOfQueryVertex; i++)
                {
                    if ((i != best_nid) && (node_flags[i] == NS_CNEIGH))
                    {
                        if (are_core_compatible(best_nid, i, node_flags, ssg))
                        {
                            cascade.insert(i);
                        }
                    }
                }
            }

            push_node_to_core(best_nid, si, node_flags, ssg, StateToQueryVertex, QueryVertexToState);

            std::cout << "core compatible: ";
            for (auto& i : cascade)
            {
                std::cout << i << " ";
            }
            std::cout << "\n";
            int osi = si;
            for (auto& i : cascade)
            {
                si++;
                push_node_to_core(i, si, node_flags, ssg, StateToQueryVertex, QueryVertexToState);
                ParentState[si] = osi;
            }

            for (int i = 0; i < NumOfQueryVertex; i++)
            {
                std::cout << i << "[" << node_flags[i] << "] ";
            }
            std::cout << "\n";
        }

        int e_count, o_e_count, i_e_count, n, nn;
        for (int si = 0; si < NumOfQueryVertex; si++)
        {

            n = StateToQueryVertex[si];
            e_count = 0;
            o_e_count = 0;
            for (int i = 0; i < ssg.OutAdjSizes[n]; i++)
            {
                if (QueryVertexToState[ssg.OutAdjList[n][i]] < si)
                {
                    e_count++;
                    o_e_count++;
                }
            }
            i_e_count = 0;
            for (int i = 0; i < ssg.InAdjSizes[n]; i++)
            {
                if (QueryVertexToState[ssg.InAdjList[n][i]] < si)
                {
                    e_count++;
                    i_e_count++;
                }
            }

            EdgeSizes[si] = e_count;
            OutEdgeSizes[si] = o_e_count;
            InEdgeSizes[si] = i_e_count;

            Edges[si] = new FMatchingMachineEdge[e_count];

            if (e_count > 0)
            {

                e_count = 0;

                for (int i = 0; i < ssg.OutAdjSizes[n]; i++)
                {
                    if (QueryVertexToState[ssg.OutAdjList[n][i]] < si)
                    {
                        Edges[si][e_count].Source = QueryVertexToState[n];
                        Edges[si][e_count].Target = QueryVertexToState[ssg.OutAdjList[n][i]];
                        e_count++;
                    }
                }
                for (int i = 0; i < ssg.InAdjSizes[n]; i++)
                {
                    if (QueryVertexToState[ssg.InAdjList[n][i]] < si)
                    {
                        Edges[si][e_count].Target = QueryVertexToState[n];
                        Edges[si][e_count].Source = QueryVertexToState[ssg.InAdjList[n][i]];
                        e_count++;
                    }
                }
            }
        }
    }

private:
    void push_node_to_core(int nid, int si, NodeFlag* node_flags, FRiGraph& qg, int* map_state_to_node, int* map_node_to_state)
    {
        node_flags[nid] = NS_CORE;

        for (int i = 0; i < qg.OutAdjSizes[nid]; i++)
        {
            if (node_flags[qg.OutAdjList[nid][i]] == NS_UNV)
            {
                node_flags[qg.OutAdjList[nid][i]] = NS_CNEIGH;
            }
        }
        for (int i = 0; i < qg.InAdjSizes[nid]; i++)
        {
            if (node_flags[qg.InAdjList[nid][i]] == NS_UNV)
            {
                node_flags[qg.InAdjList[nid][i]] = NS_CNEIGH;
            }
        }
        map_state_to_node[si] = nid;
        map_node_to_state[nid] = si;
    }

    void get_scores(int nid, int* scores, NodeFlag* node_flags, FRiGraph& qg)
    {
        std::set<int> all;

        std::set<int> cores;
        std::set<int> neighs;
        std::set<int> unvs;

        for (int i = 0; i < qg.OutAdjSizes[nid]; i++)
        {
            if (node_flags[qg.OutAdjList[nid][i]] == NS_CORE)
            {
                cores.insert(qg.OutAdjList[nid][i]);
            }
            else if (node_flags[qg.OutAdjList[nid][i]] == NS_CNEIGH)
            {
                neighs.insert(qg.OutAdjList[nid][i]);
            }
            else
            {
                unvs.insert(qg.OutAdjList[nid][i]);
            }
            all.insert(qg.OutAdjList[nid][i]);
        }

        for (int i = 0; i < qg.InAdjSizes[nid]; i++)
        {
            if (node_flags[qg.InAdjList[nid][i]] == NS_CORE)
            {
                cores.insert(qg.InAdjList[nid][i]);
            }
            else if (node_flags[qg.InAdjList[nid][i]] == NS_CNEIGH)
            {
                neighs.insert(qg.InAdjList[nid][i]);
            }
            else
            {
                unvs.insert(qg.InAdjList[nid][i]);
            }
            all.insert(qg.InAdjList[nid][i]);
        }

        scores[0] = cores.size();
        scores[1] = neighs.size();
        scores[2] = unvs.size();
        scores[3] = all.size();
        scores[4] = domains_size[nid];
    }

    int compare_scores(int* s1, int* s2, int n1, int n2)
    {
        for (int i = 0; i < 4; i++)
        {
            if (s1[i] != s2[i])
            {
                return s1[i] - s2[i];
            }
        }
        if (s1[4] != s2[4])
        {
            return s2[4] - s1[4];
        }

        return n2 - n1;
    }

    bool are_core_compatible(int nid1, int nid2, NodeFlag* node_flags, FRiGraph& qg)
    {
        if (nodeComparator.compare(qg.VertexAttributes[nid1], qg.VertexAttributes[nid2]))
        {
            bool found;

            for (int i = 0; i < qg.OutAdjSizes[nid1]; i++)
            {
                if (node_flags[qg.OutAdjList[nid1][i]] == NS_CORE)
                {
                    found = false;
                    for (int j = 0; j < qg.OutAdjSizes[nid2]; j++)
                    {
                        if (qg.OutAdjList[nid2][j] == qg.OutAdjList[nid1][i])
                        {
                            if (edgeComparator.compare(qg.OutAdjAttributes[nid2][j], qg.OutAdjAttributes[nid1][i]))
                            {
                                found = true;
                                break;
                            }
                        }
                    }
                    if (!found)
                    {
                        return false;
                    }
                }
            }

            for (int c = 0; c < qg.NumOfVertex; c++)
            {
                if (node_flags[c] == NS_CORE)
                {
                    for (int i = 0; i < qg.OutAdjSizes[c]; i++)
                    {
                        if (qg.OutAdjList[c][i] == nid1)
                        {

                            found = false;
                            for (int j = 0; j < qg.OutAdjSizes[c]; j++)
                            {
                                if (qg.OutAdjList[c][j] == nid2)
                                {
                                    if (edgeComparator.compare(qg.OutAdjAttributes[c][j], qg.OutAdjAttributes[c][i]))
                                    {
                                        found = true;
                                        break;
                                    }
                                }
                            }
                            if (!found)
                            {
                                return false;
                            }
                        }
                    }
                }
            }

            return true;
        }
        return false;
    }
};
} // namespace rilib

#endif /* MAMACONSTRFIRSTNSCC_H_ */
