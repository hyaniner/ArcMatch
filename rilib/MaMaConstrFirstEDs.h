/*
 * MaMaConstrFirst.h
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

#ifndef MAMACONSTRFIRSTEDS_H_
#define MAMACONSTRFIRSTEDS_H_

#include "Domains.h"
#include "Graph.h"
#include "MatchingMachine.h"
#include "sbitset.h"

namespace rilib
{
class FAmMaMaConstrFirstEDs : public FRiMatchingMachine
{
    FAmsbitset* domains;
    int* domains_size;
    FAmEdgeDomains& edge_domains;

public:
    FAmMaMaConstrFirstEDs(FRiGraph& query, FAmsbitset* _domains, int* _domains_size, FAmEdgeDomains& _edomains)
        : FRiMatchingMachine(query)
        , domains(_domains)
        , domains_size(_domains_size)
        , edge_domains(_edomains)
    {
    }

    virtual void Build(FRiGraph& ssg)
    {

#ifdef MDEBUG
        std::cout << "init mama...\n";
#endif

        enum NodeFlag { NS_CORE, NS_CNEIGH, NS_UNV };
        NodeFlag* node_flags = new NodeFlag[NumOfQueryVertex]; // indexed by node_id
        int** weights = new int*[NumOfQueryVertex]; // indexed by node_id
        int* t_parent_node = (int*)calloc(NumOfQueryVertex, sizeof(int)); // indexed by node_id
        EParentType* t_parent_type = new EParentType[NumOfQueryVertex]; // indexed by node id

        int nof_single_domains = 0;

        double** o_query_e_weights = new double*[NumOfQueryVertex];
        for (int i = 0; i < ssg.NumOfVertex; i++)
        {
            o_query_e_weights[i] = new double[ssg.OutAdjSizes[i]];
            for (int j = 0; j < ssg.OutAdjSizes[i]; j++)
            {
                o_query_e_weights[i][j] = 1;
            }
        }
        double** i_query_e_weights = new double*[NumOfQueryVertex];
        for (int i = 0; i < ssg.NumOfVertex; i++)
        {
            i_query_e_weights[i] = new double[ssg.InAdjSizes[i]];
            for (int j = 0; j < ssg.InAdjSizes[i]; j++)
            {
                i_query_e_weights[i][j] = 1;
            }
        }

        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            node_flags[i] = NS_UNV;
            weights[i] = new int[3];
            weights[i][0] = 0;
            weights[i][1] = 0;

            weights[i][2] = 0.0;
            for (int j = 0; j < ssg.OutAdjSizes[i]; j++)
            {
                weights[i][2] += o_query_e_weights[i][j];
            }
            for (int j = 0; j < ssg.InAdjSizes[i]; j++)
            {
                weights[i][2] += i_query_e_weights[i][j];
            }

            t_parent_node[i] = -1;
            t_parent_type[i] = PARENTTYPE_NULL;

            if (domains_size[i] == 1)
                nof_single_domains++;
        }

        int si = 0;
        int n;
        int nIT;
        int ni;
        int nnIT;
        int nni;
        int nqueueL = 0, nqueueR = 0;
        int maxi, maxv;
        int tmp;

#ifdef MDEBUG
        std::cout << "single domains [" << nof_single_domains << "]...\n";
#endif

        if (nof_single_domains != 0)
        {

            nqueueR = nof_single_domains;

            for (int n = 0; n < NumOfQueryVertex; n++)
            {
                if (domains_size[n] == 1)
                {

#ifdef MDEBUG
                    std::cout << "nqueueR(" << nqueueR << ") node[" << n << "] si[" << si << "]\n";
#endif

                    StateToQueryVertex[si] = n;
                    QueryVertexToState[n] = si;
                    t_parent_type[n] = PARENTTYPE_NULL;
                    t_parent_node[n] = -1;

                    // update nodes' flags & weights
                    node_flags[n] = NS_CORE;
                    nIT = 0;
                    while (nIT < ssg.OutAdjSizes[n])
                    {
                        ni = ssg.OutAdjList[n][nIT];
                        if (ni != n && domains_size[ni] > 1)
                        {
                            weights[ni][0] += o_query_e_weights[n][nIT];
                            weights[ni][1] -= o_query_e_weights[n][nIT];

                            if (node_flags[ni] == NS_UNV)
                            {
                                node_flags[ni] = NS_CNEIGH;
                                t_parent_node[ni] = n;
                                t_parent_type[ni] = PARENTTYPE_OUT;
                                // add to queue
                                StateToQueryVertex[nqueueR] = ni;
                                QueryVertexToState[ni] = nqueueR;
                                nqueueR++;

                                nnIT = 0;
                                while (nnIT < ssg.OutAdjSizes[ni])
                                {
                                    nni = ssg.OutAdjList[ni][nnIT];

                                    weights[nni][1] += o_query_e_weights[n][nnIT];

                                    nnIT++;
                                }
                            }
                        }
                        nIT++;
                    }

                    nIT = 0;
                    while (nIT < ssg.InAdjSizes[n])
                    {
                        ni = ssg.InAdjList[n][nIT];
                        if (ni != n && domains_size[ni] > 1)
                        {

                            weights[ni][0] += i_query_e_weights[n][nIT];
                            weights[ni][1] -= i_query_e_weights[n][nIT];

                            if (node_flags[ni] == NS_UNV)
                            {
                                node_flags[ni] = NS_CNEIGH;
                                t_parent_node[ni] = n;
                                t_parent_type[ni] = PARENTTYPE_IN;
                                // add to queue
                                StateToQueryVertex[nqueueR] = ni;
                                QueryVertexToState[ni] = nqueueR;
                                nqueueR++;

                                nnIT = 0;
                                while (nnIT < ssg.InAdjSizes[ni])
                                {
                                    nni = ssg.InAdjList[ni][nnIT];
                                    weights[nni][1] += i_query_e_weights[n][nnIT];
                                    nnIT++;
                                }
                            }
                        }
                        nIT++;
                    }

                    si++;
                }
            }
        }
        nqueueL = nof_single_domains;

#ifdef MDEBUG
        std::cout << "others...\n";
#endif
        while (si < NumOfQueryVertex)
        {

            if (nqueueL == nqueueR)
            {
                // if queue is empty....
                maxi = -1;
                maxv = -1;
                nIT = 0;
                while (nIT < NumOfQueryVertex)
                {
                    if (node_flags[nIT] == NS_UNV && weights[nIT][2] > maxv)
                    {
                        maxv = weights[nIT][2];
                        maxi = nIT;
                    }
                    nIT++;
                }
                StateToQueryVertex[si] = maxi;
                QueryVertexToState[maxi] = si;
                t_parent_type[maxi] = PARENTTYPE_NULL;
                t_parent_node[maxi] = -1;

                nqueueR++;

                n = maxi;
                nIT = 0;
                while (nIT < ssg.OutAdjSizes[n])
                {
                    ni = ssg.OutAdjList[n][nIT];
                    if (ni != n)
                    {
                        weights[ni][1] += o_query_e_weights[n][nIT];
                    }
                    nIT++;
                }
                while (nIT < ssg.InAdjSizes[n])
                {
                    ni = ssg.InAdjList[n][nIT];
                    if (ni != n)
                    {
                        weights[ni][1] += i_query_e_weights[n][nIT];
                    }
                    nIT++;
                }
            }

            if (nqueueL != nqueueR - 1)
            {
                maxi = nqueueL;
                for (int mi = maxi + 1; mi < nqueueR; mi++)
                {
                    if (wcompare(StateToQueryVertex[mi], StateToQueryVertex[maxi], weights) < 0)
                    {
                        maxi = mi;
                    }
                }
                tmp = StateToQueryVertex[nqueueL];
                StateToQueryVertex[nqueueL] = StateToQueryVertex[maxi];
                StateToQueryVertex[maxi] = tmp;
            }

            n = StateToQueryVertex[si];
            QueryVertexToState[n] = si;

            // move queue left limit
            nqueueL++;
            // update nodes' flags & weights
            node_flags[n] = NS_CORE;
            nIT = 0;
            while (nIT < ssg.OutAdjSizes[n])
            {
                ni = ssg.OutAdjList[n][nIT];
                if (ni != n)
                {
                    weights[ni][0] += o_query_e_weights[n][nIT];
                    weights[ni][1] -= o_query_e_weights[n][nIT];

                    if (node_flags[ni] == NS_UNV)
                    {
                        node_flags[ni] = NS_CNEIGH;
                        t_parent_node[ni] = n;

                        t_parent_type[ni] = PARENTTYPE_OUT;

                        StateToQueryVertex[nqueueR] = ni;
                        QueryVertexToState[ni] = nqueueR;
                        nqueueR++;

                        nnIT = 0;
                        while (nnIT < ssg.OutAdjSizes[ni])
                        {
                            nni = ssg.OutAdjList[ni][nnIT];

                            weights[nni][1] += o_query_e_weights[n][nIT];
                            nnIT++;
                        }
                    }
                }
                nIT++;
            }

            nIT = 0;
            while (nIT < ssg.InAdjSizes[n])
            {
                ni = ssg.InAdjList[n][nIT];
                if (ni != n)
                {

                    weights[ni][0] += i_query_e_weights[n][nIT];
                    weights[ni][1] -= i_query_e_weights[n][nIT];

                    if (node_flags[ni] == NS_UNV)
                    {
                        node_flags[ni] = NS_CNEIGH;
                        t_parent_node[ni] = n;

                        t_parent_type[ni] = PARENTTYPE_IN;

                        // add to queue
                        StateToQueryVertex[nqueueR] = ni;
                        QueryVertexToState[ni] = nqueueR;
                        nqueueR++;

                        nnIT = 0;
                        while (nnIT < ssg.InAdjSizes[ni])
                        {
                            nni = ssg.InAdjList[ni][nnIT];
                            weights[nni][1] += i_query_e_weights[n][nIT];
                            nnIT++;
                        }
                    }
                }
                nIT++;
            }

            si++;
        }

        int e_count, o_e_count, i_e_count, nn;
        for (int si = 0; si < NumOfQueryVertex; si++)
        {

            n = StateToQueryVertex[si];

            if (t_parent_node[n] != -1)
                ParentState[si] = QueryVertexToState[t_parent_node[n]];
            else
                ParentState[si] = -1;
            ParentType[si] = t_parent_type[n];

#ifdef MDEBUG
            if (ParentType[si] != PARENTTYPE_NULL)
            {
                std::cout << "P:" << si << " " << n << " " << ParentState[si] << " " << StateToQueryVertex[ParentState[si]];
                if (ParentType[si] == PARENTTYPE_OUT)
                    std::cout << " out\n";
                else
                    std::cout << " in\n";
            }
            else
            {
                std::cout << "P:" << si << " " << n << " " << ParentState[si] << " NULL\n";;
            }
#endif
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

        delete[] node_flags;
        for (int i = 0; i < NumOfQueryVertex; i++)
            delete[] weights[i];
        delete[] weights;
        free(t_parent_node);
        delete[] t_parent_type;
    }

private:
    int wcompare(int i, int j, int** weights)
    {
        for (int w = 0; w < 3; w++)
        {
            if (weights[i][w] != weights[j][w])
            {
                return weights[j][w] - weights[i][w];
            }
        }
        return i - j;
    }
};
} // namespace rilib

#endif /* MAMACONSTRFIRST_H_ */
