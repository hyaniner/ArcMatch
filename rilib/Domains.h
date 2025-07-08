/*
 * Domains.h
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

#ifndef DOMAINS_H_
#define DOMAINS_H_

#include "AttributeComparator.h"
#include "Graph.h"
#include "sbitset.h"

#include <set>
#include <tuple>
#include <unordered_set>

namespace rilib
{
bool init_domains(FRiGraph& target, FRiGraph& pattern, FAmAttributeComparator& nodeComparator, FAmAttributeComparator& edgeComparator, FAmsbitset* domains, bool iso)
{

    if (iso)
    {
        for (int q = 0; q < pattern.NumOfVertex; q++)
        {
            for (int r = target.NumOfVertex - 1; r >= 0; r--)
            {
                if (target.OutAdjSizes[r] == pattern.OutAdjSizes[q] && target.InAdjSizes[r] == pattern.InAdjSizes[q] && nodeComparator.compare(target.VertexAttributes[r], pattern.VertexAttributes[q]))
                {
                    domains[q].set(r, true);
                }
            }
            if (domains[q].is_empty())
            {
                return false;
            }
        }
    }
    else
    {
        for (int q = 0; q < pattern.NumOfVertex; q++)
        {
            for (int r = target.NumOfVertex - 1; r >= 0; r--)
            {
                if (target.OutAdjSizes[r] >= pattern.OutAdjSizes[q] && target.InAdjSizes[r] >= pattern.InAdjSizes[q] && nodeComparator.compare(target.VertexAttributes[r], pattern.VertexAttributes[q]))
                {
                    domains[q].set(r, true);
                }
            }
            if (domains[q].is_empty())
            {
                return false;
            }
        }
    }

    int ra, qb, rb;
    bool notfound;

    // 1°-level neighborhood and edges labels
    for (int qa = 0; qa < pattern.NumOfVertex; qa++)
    {

        for (FAmsbitset::iterator qaIT = domains[qa].first_ones(); qaIT != domains[qa].end(); qaIT.next_ones())
        {
            ra = qaIT.first;
            // for each edge qa->qb  check if exists ra->rb
            for (int i_qb = 0; i_qb < pattern.OutAdjSizes[qa]; i_qb++)
            {
                qb = pattern.OutAdjList[qa][i_qb];
                notfound = true;

                for (int i_rb = 0; i_rb < target.OutAdjSizes[ra]; i_rb++)
                {
                    rb = target.OutAdjList[ra][i_rb];
                    if (domains[qb].get(rb) && edgeComparator.compare(pattern.OutAdjAttributes[qa][i_qb], target.OutAdjAttributes[ra][i_rb]))
                    {
                        notfound = false;
                        break;
                    }
                }

                if (notfound)
                {
                    domains[qa].set(ra, false);
                    break;
                }
            }
        }

        if (domains[qa].is_empty())
            return false;
    }

#ifdef NODE_D_CONV
    bool changes = true;
    while (changes)
    {
        changes = false;
        for (int qa = 0; qa < pattern.NumOfVertex; qa++)
        {
            for (FAmsbitset::iterator qaIT = domains[qa].first_ones(); qaIT != domains[qa].end(); qaIT.next_ones())
            {
                ra = qaIT.first;
                // fore each edge qa->qb  check if exists ra->rb
                for (int i_qb = 0; i_qb < pattern.OutAdjSizes[qa]; i_qb++)
                {
                    qb = pattern.OutAdjList[qa][i_qb];
                    notfound = true;
                    for (int i_rb = 0; i_rb < target.OutAdjSizes[ra]; i_rb++)
                    {
                        rb = target.OutAdjList[ra][i_rb];
                        if (domains[qb].get(rb) && edgeComparator.compare(pattern.OutAdjAttributes[qa][i_qb], target.OutAdjAttributes[ra][i_rb]))
                        {
                            notfound = false;
                            break;
                        }
                    }

                    if (notfound)
                    {
                        domains[qa].set(ra, false);
                        changes = true;
                    }
                }
            }
            if (domains[qa].is_empty())
                return false;
        }
    }
#endif

    return true;
};

struct pair_hash
{
    inline std::size_t operator()(const std::pair<int, int>& v) const
    {
        return v.first * 31 + v.second;
    }
};

typedef std::unordered_set<std::pair<int, int>, pair_hash> unordered_edge_set;

class FAmEdgeDomains
{
public:
    int** pattern_out_adj_eids;
    int* inv_pattern_edge_ids;
    int* inv_target_edge_ids;
    int nof_pattern_edges;
    int nof_target_edges;

    int** pattern_in_adj_eids;

    unordered_edge_set* domains;

    FAmEdgeDomains()
    {
    };
};

bool init_edomains(FRiGraph& target, FRiGraph& pattern, FAmsbitset* node_domains, FAmAttributeComparator& edgeComparator, FAmEdgeDomains& edomains)
{
    int nof_pedges = 0;
    for (int i = 0; i < pattern.NumOfVertex; i++)
    {
        nof_pedges += pattern.OutAdjSizes[i];
    }

    int nof_tedges = 0;
    for (int i = 0; i < target.NumOfVertex; i++)
    {
        nof_tedges += target.OutAdjSizes[i];
    }

    edomains.nof_pattern_edges = nof_pedges;
    edomains.nof_target_edges = nof_tedges;

    edomains.pattern_out_adj_eids = new int*[pattern.NumOfVertex];

    edomains.inv_pattern_edge_ids = new int[nof_pedges * 2];

    edomains.inv_target_edge_ids = new int[nof_tedges * 2];

    int n = 0;
    for (int i = 0; i < pattern.NumOfVertex; i++)
    {
        edomains.pattern_out_adj_eids[i] = new int[pattern.OutAdjSizes[i]];
        for (int j = 0; j < pattern.OutAdjSizes[i]; j++)
        {
            edomains.inv_pattern_edge_ids[n * 2] = i;
            edomains.inv_pattern_edge_ids[(n * 2) + 1] = j;
            edomains.pattern_out_adj_eids[i][j] = n;
            n++;
        }
    }

    edomains.pattern_in_adj_eids = new int*[pattern.NumOfVertex];
    for (int i = 0; i < pattern.NumOfVertex; i++)
    {
        edomains.pattern_in_adj_eids[i] = new int[pattern.InAdjSizes[i]];

        for (int j = 0; j < pattern.InAdjSizes[i]; j++)
        {
            int nn = pattern.InAdjList[i][j];

            for (int ni = 0; ni < pattern.OutAdjSizes[nn]; ni++)
            {
                if (pattern.OutAdjList[nn][ni] == i)
                {
                    edomains.pattern_in_adj_eids[i][j] = edomains.pattern_out_adj_eids[nn][ni];
                    break;
                }
            }
        }
    }

    n = 0;
    for (int i = 0; i < target.NumOfVertex; i++)
    {
        for (int j = 0; j < target.OutAdjSizes[i]; j++)
        {
            edomains.inv_target_edge_ids[n * 2] = i;
            edomains.inv_target_edge_ids[(n * 2) + 1] = j;
            n++;
        }
    }

    edomains.domains = new unordered_edge_set[nof_pedges];

    /*given a pattern edge, we search for compatible target edges.
    Given a pattern node p_s and a negihborn of it p_t, we look at the aleady computed domain of it.
    Given a target node t_s in the domain of p_s,
    we search for a neighborn t_t of t_s in the domain of p_t.
    The search is made by checking fo edge label compatibility.
    */

    for (int ps = 0; ps < pattern.NumOfVertex; ps++)
    {
        for (int ps_n = 0; ps_n < pattern.OutAdjSizes[ps]; ps_n++)
        {

            int pt = pattern.OutAdjList[ps][ps_n];

            for (FAmsbitset::iterator psIT = node_domains[ps].first_ones(); psIT != node_domains[ps].end(); psIT.next_ones())
            {

                int ts = psIT.first;

                for (FAmsbitset::iterator ptIT = node_domains[pt].first_ones(); ptIT != node_domains[pt].end(); ptIT.next_ones())
                {

                    int tt = ptIT.first;

                    for (int ts_n = 0; ts_n < target.OutAdjSizes[ts]; ts_n++)
                    {
                        if (tt == target.OutAdjList[ts][ts_n])
                        {
                            if (edgeComparator.compare(pattern.OutAdjAttributes[ps][ps_n], target.OutAdjAttributes[ts][ts_n]))
                            {

                                edomains.domains[edomains.pattern_out_adj_eids[ps][ps_n]].insert(std::pair<int, int>(ts, tt));
                            }
                        }
                    }
                }
            }
        }
    }

    return true;
};

class FAmDomainReduction
{
public:
    FRiGraph& query;
    FAmsbitset* node_domains;
    FAmEdgeDomains& edge_domains;
    int nof_target_nodes;

    FAmDomainReduction(FRiGraph& _query, FAmsbitset* ndomains, FAmEdgeDomains& edomains, int noftargetnodes)
        : query(_query)
        , node_domains(ndomains)
        , edge_domains(edomains)
        , nof_target_nodes(noftargetnodes)
    {
    };

    bool refine_domains(int altered_q_node)
    {
        bool erased = false;

        int n;
        for (int ni = 0; ni < query.OutAdjSizes[altered_q_node]; ni++)
        {
            n = query.OutAdjList[altered_q_node][ni];
            unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_out_adj_eids[altered_q_node][ni]]);

            for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end();)
            {
                if (node_domains[altered_q_node].get(eit->first) == false)
                {
                    eset->erase(eit++);
                    erased = true;
                }
                else
                {
                    ++eit;
                }
            }
        }

        for (int ni = 0; ni < query.InAdjSizes[altered_q_node]; ni++)
        {
            n = query.InAdjList[altered_q_node][ni];
            unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_in_adj_eids[altered_q_node][ni]]);

            for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end();)
            {
                if (node_domains[altered_q_node].get(eit->second) == false)
                {
                    eset->erase(eit++);
                    erased = true;
                }
                else
                {
                    ++eit;
                }
            }
        }

#ifdef EDGE_D_CONV
        if (erased)
        {
            while (erased)
            {
                erased = false;

                for (int s = 0; s < query.NumOfVertex; s++)
                {
                    for (int ni = 0; ni < query.OutAdjSizes[s]; ni++)
                    {
                        n = query.OutAdjList[s][ni];
                        unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_out_adj_eids[s][ni]]);
                        for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end();)
                        {
                            if ((node_domains[s].get(eit->first) == false) || (node_domains[n].get(eit->second) == false))
                            {
                                eset->erase(eit++);
                            }
                            else
                            {
                                ++eit;
                            }
                        }
                    }
                }

                for (int ni = 0; ni < query.InAdjSizes[altered_q_node]; ni++)
                {
                    n = query.InAdjList[altered_q_node][ni];
                    unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_in_adj_eids[altered_q_node][ni]]);

                    for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end();)
                    {
                        if (node_domains[altered_q_node].get(eit->second) == false)
                        {
                            eset->erase(eit++);
                            erased = true;
                        }
                        else
                        {
                            ++eit;
                        }
                    }
                }

                for (int s = 0; s < query.NumOfVertex; s++)
                {

                    for (int ni = 0; ni < query.OutAdjSizes[s]; ni++)
                    {
                        n = query.OutAdjList[s][ni];

                        std::set<int> source_set;
                        unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_out_adj_eids[s][ni]]);
                        for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++)
                        {
                            source_set.insert(eit->first);
                        }

                        for (FAmsbitset::iterator dit = node_domains[s].first_ones(); dit != node_domains[s].end(); dit.next_ones())
                        {

                            if (source_set.find(dit.first) == source_set.end())
                            {
                                node_domains[s].set(dit.first, false);
                                erased = true;
                            }
                        }
                    }

                    for (int ni = 0; ni < query.InAdjSizes[s]; ni++)
                    {
                        n = query.InAdjList[s][ni];

                        std::set<int> source_set;
                        unordered_edge_set* eset = &(edge_domains.domains[edge_domains.pattern_in_adj_eids[s][ni]]);
                        for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++)
                        {
                            source_set.insert(eit->second);
                        }

                        for (FAmsbitset::iterator dit = node_domains[s].first_ones(); dit != node_domains[s].end(); dit.next_ones())
                        {

                            if (source_set.find(dit.first) == source_set.end())
                            {
                                node_domains[s].set(dit.first, false);
                                erased = true;
                            }
                        }
                    }
                }
            }
        }
#endif
        return true;
    };

    bool verify_path_dfs(int* q_dfs, int* q_dfs_adji, int q_level, int* c_dfs, bool* c_visited, int c_level, bool circle)
    {
        if (c_level == q_level - 1)
        {
            if (circle)
            {
                unordered_edge_set eset = edge_domains.domains[edge_domains.pattern_out_adj_eids[q_dfs[c_level]][q_dfs_adji[c_level + 1]]];
                for (unordered_edge_set::iterator eit = eset.begin(); eit != eset.end(); eit++)
                {
                    if ((eit->first == c_dfs[c_level]) && (eit->second == c_dfs[0]))
                    {
                        return true;
                    }
                }
                return false;
            }
            else
            {
                unordered_edge_set eset = edge_domains.domains[edge_domains.pattern_out_adj_eids[q_dfs[c_level]][q_dfs_adji[c_level + 1]]];
                for (unordered_edge_set::iterator eit = eset.begin(); eit != eset.end(); eit++)
                {
                    if ((eit->first == c_dfs[c_level]) && (!c_visited[eit->second]))
                    {
                        return true;
                    }
                }
                return false;
            }
        }
        else
        {
            bool found = false;
            unordered_edge_set eset = edge_domains.domains[edge_domains.pattern_out_adj_eids[q_dfs[c_level]][q_dfs_adji[c_level + 1]]];
            for (unordered_edge_set::iterator eit = eset.begin(); eit != eset.end(); eit++)
            {
                if ((eit->first == c_dfs[c_level]) && (!c_visited[eit->second]))
                {
                    c_dfs[c_level + 1] = eit->second;
                    c_visited[eit->second] = true;
                    found |= verify_path_dfs(q_dfs, q_dfs_adji, q_level, c_dfs, c_visited, c_level + 1, circle);
                    c_visited[eit->second] = false;
                    if (found)
                    {
                        break;
                    }
                }
            }
            return found;
        }
        return false;
    };

    void verify_path(int* q_dfs, int* q_dfs_adji, int q_level, bool circle)
    {
        if (q_level > 1)
        {
            int* c_dfs = new int[query.NumOfVertex];
            bool* c_visited = new bool[nof_target_nodes];
            for (int i = 0; i < nof_target_nodes; i++)
            {
                c_visited[i] = false;
            }

            bool removed = false;
            int cnode;
            for (FAmsbitset::iterator dit = node_domains[q_dfs[0]].first_ones(); dit != node_domains[q_dfs[0]].end(); dit.next_ones())
            {
                cnode = dit.first;
                c_dfs[0] = cnode;
                c_visited[cnode] = true;

                if (!verify_path_dfs(q_dfs, q_dfs_adji, q_level, c_dfs, c_visited, 0, circle))
                {
                    node_domains[q_dfs[0]].set(cnode, false);
                    removed = true;
                }

                c_visited[cnode] = false;
            }
            if (removed)
            {
                refine_domains(q_dfs[0]);
            }

            delete[] c_visited;
            delete[] c_dfs;
        }
    };

    /*
    This is a temporary version which works for undirected graphs
    because only out edges are visited during the DFS visit.
    */
    void reduce_by_paths_dfs(int* dfs, int* dfs_adji, bool* visited, int level, int max_lp)
    {
        int n;
        int nof_p = 0;
        for (int ni = 0; ni < query.OutAdjSizes[dfs[level]]; ni++)
        {
            n = query.OutAdjList[dfs[level]][ni];

            if ((!visited[n]))
            {
                if (level == query.NumOfVertex - 2)
                {
                    nof_p++;

                    dfs[level + 1] = n;
                    dfs_adji[level + 1] = ni;
                    verify_path(dfs, dfs_adji, level + 1, false);
                }
                else
                {
                    nof_p++;
                    dfs[level + 1] = n;
                    dfs_adji[level + 1] = ni;

                    if (level == max_lp)
                    {
                        verify_path(dfs, dfs_adji, level + 1, false);
                    }
                    else
                    {
                        visited[n] = true;
                        reduce_by_paths_dfs(dfs, dfs_adji, visited, level + 1, max_lp);
                        visited[n] = false;
                    }
                }
            }
            else if ((level > 0) && (n != dfs[level - 1]) && (n == dfs[0]))
            {
                nof_p++;
                dfs[level + 1] = n;
                dfs_adji[level + 1] = ni;

                verify_path(dfs, dfs_adji, level + 1, true);
            }
        }
        if ((level > 0) && (nof_p == 0))
        {
            verify_path(dfs, dfs_adji, level - 1, false);
        }
    };

    void reduce_by_paths(int starting_node, int max_lp)
    {
        int* dfs = new int[query.NumOfVertex];
        int* dfs_adji = new int[query.NumOfVertex];
        bool* visited = new bool[query.NumOfVertex];
        for (int i = 0; i < query.NumOfVertex; i++)
        {
            visited[i] = false;
        }
        dfs[0] = starting_node;
        dfs_adji[0] = -1;
        visited[starting_node] = true;

        reduce_by_paths_dfs(dfs, dfs_adji, visited, 0, max_lp);

        delete[] visited;
        delete[] dfs;
        delete[] dfs_adji;
    };

    void reduce_by_paths(int max_lp)
    {
        for (int i = 0; i < query.NumOfVertex; i++)
        {
            reduce_by_paths(i, max_lp);
        }
    };

    void final_refinement()
    {

        for (int n = 0; n < query.NumOfVertex; n++)
        {

            for (FAmsbitset::iterator nnIT = node_domains[n].first_ones(); nnIT != node_domains[n].end(); nnIT.next_ones())
            {
                int dn = nnIT.first;

                bool found;

                for (int ni = 0; ni < query.OutAdjSizes[n]; ni++)
                {
                    found = false;

                    unordered_edge_set eset = edge_domains.domains[edge_domains.pattern_out_adj_eids[n][ni]];
                    for (unordered_edge_set::iterator eit = eset.begin(); eit != eset.end(); eit++)
                    {
                        if (dn == eit->first)
                        {
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        node_domains[n].set(dn, false);
                        break;
                    }
                }

                for (int ni = 0; ni < query.InAdjSizes[n]; ni++)
                {
                    found = false;

                    unordered_edge_set eset = edge_domains.domains[edge_domains.pattern_in_adj_eids[n][ni]];
                    for (unordered_edge_set::iterator eit = eset.begin(); eit != eset.end(); eit++)
                    {
                        if (dn == eit->second)
                        {
                            found = true;
                            break;
                        }
                    }

                    if (!found)
                    {
                        node_domains[n].set(dn, false);
                        break;
                    }
                }
            }
        }
    };
};

void print_domains(FRiGraph& query, FRiGraph& target, FAmsbitset* node_domains, FAmEdgeDomains& edge_domains)
{
    std::cout << "nof query nodes " << query.NumOfVertex << "\n";
    for (int i = 0; i < query.NumOfVertex; i++)
    {
        std::cout << "node domain " << i << ":" << node_domains[i].count_ones() << ": ";
        for (FAmsbitset::iterator it = node_domains[i].first_ones(); it != node_domains[i].end(); it.next_ones())
        {
        std:
            cout << it.first << " ";
        }
        std::cout << "\n";
    }
    for (int n = 0; n < query.NumOfVertex; n++)
    {
        for (int ni = 0; ni < query.OutAdjSizes[n]; ni++)
        {
            int eid = edge_domains.pattern_out_adj_eids[n][ni];
            std::cout << "edge domain: " << n << "-" << query.OutAdjList[n][ni] << ":eid " << eid << ":" << edge_domains.domains[eid].size() << "\n";
        }
    }
};

void print_domains_extended(FRiGraph& query, FRiGraph& target, FAmsbitset* node_domains, FAmEdgeDomains& edge_domains)
{
    std::cout << "nof query nodes " << query.NumOfVertex << "\n";
    for (int i = 0; i < query.NumOfVertex; i++)
    {
        std::cout << "node domain " << i << ":" << node_domains[i].count_ones() << ": ";
        for (FAmsbitset::iterator it = node_domains[i].first_ones(); it != node_domains[i].end(); it.next_ones())
        {
        std:
            cout << it.first << " ";
        }
        std::cout << "\n";
    }
    for (int n = 0; n < query.NumOfVertex; n++)
    {
        for (int ni = 0; ni < query.OutAdjSizes[n]; ni++)
        {
            int eid = edge_domains.pattern_out_adj_eids[n][ni];
            std::cout << "edge domain: " << n << "-" << query.OutAdjList[n][ni] << ":eid " << eid << ":" << edge_domains.domains[eid].size() << "\n";

            unordered_edge_set* eset = &(edge_domains.domains[eid]);

            for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++)
            {
                std::cout << "(" << (*eit).first << "," << (*eit).second << ")";
            }
            std::cout << "\n";
        }
    }
};
} // namespace rilib
#endif /* DOMAINS_H_ */
