/*
 * Solver.h
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

#ifndef SOLVER_H_
#define SOLVER_H_

// #define SOLVER_H_MDEBUG

#include "Domains.h"
#include "Graph.h"
#include "MatchingMachine.h"
#include "sbitset.h"

#include <unordered_map>

namespace rilib {

class Solver {
  public:
    MatchingMachine &mama;
    Graph &rgraph;
    Graph &qgraph;
    AttributeComparator &nodeComparator;
    AttributeComparator &edgeComparator;
    MatchListener &matchListener;
    sbitset *domains;
    int *domains_size;
    EdgeDomains &edomains;

    long steps;
    long triedcouples;
    long matchedcouples;

    long matchcount;

  public:
    Solver(MatchingMachine &_mama, Graph &_rgraph, Graph &_qgraph, AttributeComparator &_nodeComparator, AttributeComparator &_edgeComparator, MatchListener &_matchListener, sbitset *_domains, int *_domains_size, EdgeDomains &_edomains) : mama(_mama), rgraph(_rgraph), qgraph(_qgraph), nodeComparator(_nodeComparator), edgeComparator(_edgeComparator), matchListener(_matchListener), domains(_domains), domains_size(_domains_size), edomains(_edomains) {
        steps = 0;
        triedcouples = 0;
        matchedcouples = 0;

        matchcount = 0;
    }

    virtual ~Solver() {}

    void solve() {

        int ii;

        int nof_sn = mama.nof_sn;
        void **nodes_attrs = mama.nodes_attrs;           // indexed by state_id
        int *edges_sizes = mama.edges_sizes;             // indexed by state_id
        MaMaEdge **edges = mama.edges;                   // indexed by state_id
        int *map_node_to_state = mama.map_node_to_state; // indexed by node_id
        int *map_state_to_node = mama.map_state_to_node; // indexed by state_id
        int *parent_state = mama.parent_state;           // indexed by state_id
        MAMA_PARENTTYPE *parent_type = mama.parent_type; // indexed by state id

        int **candidates = new int *[nof_sn];  // indexed by state_id
        int *candidatesIT = new int[nof_sn];   // indexed by state_id
        int *candidatesSize = new int[nof_sn]; // indexed by state_id
        int *solution = new int[nof_sn];       // indexed by state_id
        for (ii = 0; ii < nof_sn; ii++)
            solution[ii] = -1;

        bool *matched = (bool *)calloc(rgraph.nof_nodes, sizeof(bool)); // indexed by node_id

        for (int i = 0; i < nof_sn; i++) {
            if (parent_type[i] == PARENTTYPE_NULL) {
                // std::cout<<"parent type null "<<i<<"\n";
                int n = map_state_to_node[i];
                candidates[i] = new int[domains_size[n]];

                int k = 0;
                for (sbitset::iterator IT = domains[n].first_ones(); IT != domains[n].end(); IT.next_ones()) {
                    candidates[i][k] = IT.first;
                    k++;
                }

                candidatesSize[i] = domains_size[n];
                candidatesIT[i] = -1;
            }
        }

        matchcount = 0;

        int psi = -1;
        int si = 0;
        int CandidateIndex = -1;
        int sip1;
        while (si != -1) {

            if (psi >= si) {
                matched[solution[si]] = false;
            }

            CandidateIndex = -1;
            candidatesIT[si]++;
            while (candidatesIT[si] < candidatesSize[si]) {

                CandidateIndex = candidates[si][candidatesIT[si]];
                solution[si] = CandidateIndex;

#ifdef MDEBUG
                if (!matched[CandidateIndex]) {
                    std::cout << "trying (" << map_state_to_node[si] << "," << CandidateIndex << ")\n";
                    if (!domains[map_state_to_node[si]].get(CandidateIndex))
                        std::cout << "\tfails on domains\n";
                    if (!edgesCheck(si, CandidateIndex, solution, matched))
                        std::cout << "\tfails on edges\n";
                }
#endif

                if ((!matched[CandidateIndex]) && domains[map_state_to_node[si]].get(CandidateIndex) && edgesCheck(si, CandidateIndex, solution, matched)) {
                    break;
                } else {
                    CandidateIndex = -1;
                }
                candidatesIT[si]++;
            }

            if (CandidateIndex == -1) {
                psi = si;
                si--;
            } else {
                matchedcouples++;

                if (si == nof_sn - 1) {
                    matchListener.match(nof_sn, map_state_to_node, solution);

                    matchcount++;

                    psi = si;
#ifdef FIRST_MATCH_ONLY
                    si = -1;
#endif
#ifdef FIRST_100k_MATCHES
                    if (matchcount >= 100000)
                        si = -1;
#endif
                } else {
                    matched[solution[si]] = true;
                    sip1 = si + 1;
                    if (parent_type[sip1] == PARENTTYPE_NULL) {
                    } else {
                        if (parent_type[sip1] == PARENTTYPE_IN) {
                            candidates[sip1] = rgraph.in_adj_list[solution[parent_state[sip1]]];
                            candidatesSize[sip1] = rgraph.in_adj_sizes[solution[parent_state[sip1]]];
                        } else { //(parent_type[sip1] == MAMA_PARENTTYPE::PARENTTYPE_OUT)
                            candidates[sip1] = rgraph.out_adj_list[solution[parent_state[sip1]]];
                            candidatesSize[sip1] = rgraph.out_adj_sizes[solution[parent_state[sip1]]];
                        }
                    }
                    candidatesIT[si + 1] = -1;

                    psi = si;
                    si++;
                }
            }
        }
    }

    void SolveEd() {

        int ii;

        int nof_sn = mama.nof_sn;
        int *map_node_to_state = mama.map_node_to_state; // indexed by node_id
        int *map_state_to_node = mama.map_state_to_node; // indexed by state_id
        int *parent_state = mama.parent_state;           // indexed by state_id
        MAMA_PARENTTYPE *parent_type = mama.parent_type; // indexed by state id

        int **candidates = new int *[nof_sn];         // indexed by state_id
        int **candidates_parents = new int *[nof_sn]; // indexed by state_id
        int *candidatesIT = new int[nof_sn];          // indexed by state_id
        int *candidatesSize = new int[nof_sn];        // indexed by state_id

        matchcount = 0;

        int *solution = new int[nof_sn]; // indexed by state_id
        for (ii = 0; ii < nof_sn; ii++)
            solution[ii] = -1;

        bool *matched = (bool *)calloc(rgraph.nof_nodes, sizeof(bool)); // indexed by node_id

        for (int i = 0; i < nof_sn; i++) {
#ifdef MDEBUG
            std::cout << "-----\n";
#endif
            if (parent_type[i] == PARENTTYPE_NULL) {
                int n = map_state_to_node[i];
                candidates[i] = new int[domains_size[n]];
                candidates_parents[i] = new int[domains_size[n]];

                int k = 0;
                for (sbitset::iterator IT = domains[n].first_ones(); IT != domains[n].end(); IT.next_ones()) {
                    candidates[i][k] = IT.first;
                    candidates_parents[i][k] = -1;
                    k++;
                }

                candidatesSize[i] = domains_size[n];
                candidatesIT[i] = -1;
            } else {
                int eid = -1;
                int s = 0, t = 0;
                if (parent_type[i] == PARENTTYPE_OUT) {
                    s = map_state_to_node[parent_state[i]];
                    t = map_state_to_node[i];
                } else {
                    t = map_state_to_node[parent_state[i]];
                    s = map_state_to_node[i];
                }

                for (int n = 0; n < qgraph.out_adj_sizes[s]; n++) {
                    if (qgraph.out_adj_list[s][n] == t) {
                        eid = edomains.pattern_out_adj_eids[s][n];
                    }
                }

#ifdef MDEBUG
                std::cout << s << " " << t << ":";
                if (parent_type[i] == PARENTTYPE_OUT)
                    std::cout << "out\n";
                else
                    std::cout << "in\n";
                std::cout << i << " " << eid << "<-------------------------\n";
#endif

                candidatesSize[i] = edomains.domains[eid].size();
                candidatesIT[i] = -1;
                candidates[i] = new int[candidatesSize[i]];
                candidates_parents[i] = new int[candidatesSize[i]];
                unordered_edge_set *eset = &(edomains.domains[eid]);
                int j = 0;
                for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++) {
#ifdef MDEBUG
                    std::cout << "E: " << eit->first << " -> " << eit->second << "\n";
#endif

                    if (parent_type[i] == PARENTTYPE_OUT) {
                        candidates_parents[i][j] = eit->first;
                        candidates[i][j] = eit->second;
                    } else {
                        candidates_parents[i][j] = eit->second;
                        candidates[i][j] = eit->first;
                    }
                    j++;
                }
            }
        }

#ifdef MDEBUG
        std::cout << "candidate sizes:\n";
        for (int i = 0; i < nof_sn; i++) {
            std::cout << i << ":" << map_state_to_node[i] << ":" << candidatesSize[i] << "\n";
        }
#endif

        int psi = -1;
        int si = 0;
        int CandidateIndex = -1;
        int sip1;
        while (si != -1) {

            if (psi >= si) {
                matched[solution[si]] = false;
            }

            CandidateIndex = -1;
            candidatesIT[si]++;
            while (candidatesIT[si] < candidatesSize[si]) {

                CandidateIndex = candidates[si][candidatesIT[si]];
                solution[si] = CandidateIndex;

#ifdef MDEBUG
                std::cout << "S: " << si << " " << CandidateIndex << "\n";
#endif

                if (((candidates_parents[si][candidatesIT[si]] == -1) || (candidates_parents[si][candidatesIT[si]] == solution[parent_state[si]])) && (!matched[candidates[si][candidatesIT[si]]] && edgesCheck(si, CandidateIndex, solution, matched))) {
                    break;
                } else {
                    CandidateIndex = -1;
                }

#ifdef MDEBUG
            std:
                cout << "CI: " << CandidateIndex << "\n";
#endif
                candidatesIT[si]++;
            }

            if (CandidateIndex == -1) {
                psi = si;
                si--;
            } else {
                matchedcouples++;

                if (si == nof_sn - 1) {

                    matchListener.match(nof_sn, map_state_to_node, solution);

                    matchcount++;

                    psi = si;
#ifdef FIRST_MATCH_ONLY
                    si = -1;
#endif
#ifdef FIRST_100k_MATCHES
                    if (matchcount >= 100000)
                        si = -1;
#endif
                } else {
                    matched[solution[si]] = true;
                    sip1 = si + 1;
                    candidatesIT[si + 1] = -1;
                    psi = si;
                    si++;
                }
            }
        }
    }

    // A hash function used to hash a pair of any kind
    struct hash_pair {
        template <class T1, class T2> size_t operator()(const pair<T1, T2> &p) const {
            auto hash1 = hash<T1>{}(p.first);
            auto hash2 = hash<T2>{}(p.second);
            return hash1 ^ hash2;
        }
    };

    void solve_rp() {

        int nof_sn = mama.nof_sn;
        int *map_node_to_state = mama.map_node_to_state; // indexed by node_id
        int *map_state_to_node = mama.map_state_to_node; // indexed by state_id

        typedef std::unordered_map<std::pair<int, int>, int, hash_pair> cand_ecount_t; // (tnodeid,eid) -> count

        cand_ecount_t ce_counter;
        cand_ecount_t ce_positions;

        matchcount = 0;

#ifdef MDEBUG
        std::cout << "ORDERED EDGE SETS\n";
#endif

        typedef std::set<std::pair<int, int>> ordered_edge_set;

        std::cout << "ORDERED EDGE SETS...\n";

        int **ordered_edge_domains = new int *[edomains.nof_pattern_edges];
        int *ordered_edge_domains_sizes = new int[edomains.nof_pattern_edges];

        for (int eid = 0; eid < edomains.nof_pattern_edges; eid++) {
            unordered_edge_set *eset = &(edomains.domains[eid]);
            ordered_edge_domains[eid] = new int[eset->size() * 2];
            ordered_edge_domains_sizes[eid] = eset->size();
        }
        for (int i = 0; i < nof_sn; i++) {
            for (int j = 0; j < mama.edges_sizes[i]; j++) {
                if (mama.edges[i][j].source == i) {
                    int eid = mama.edges[i][j].id;
                    ordered_edge_set tset;
                    unordered_edge_set *eset = &(edomains.domains[eid]);
                    for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++) {
                        tset.insert(std::pair<int, int>(eit->second, eit->first));

                        auto it = ce_counter.emplace(std::make_pair(eit->second, eid), 0);
                        it.first->second++;
                    }
                    int k = 0;
                    for (ordered_edge_set::iterator eit = tset.begin(); eit != tset.end(); eit++) {
                        ordered_edge_domains[eid][k] = eit->first;
                        ordered_edge_domains[eid][k + ordered_edge_domains_sizes[eid]] = eit->second;
#ifdef MDEBUG
                        std::cout << "OUT eid: " << eid << ";k " << k << ": " << eit->first << "-" << eit->second << "\n";
#endif

                        ce_positions.insert({std::make_pair(eit->first, eid), k});

                        k++;
                    }
                } else {
                    int eid = mama.edges[i][j].id;
                    ordered_edge_set tset;
                    unordered_edge_set *eset = &(edomains.domains[eid]);
                    for (unordered_edge_set::iterator eit = eset->begin(); eit != eset->end(); eit++) {
                        tset.insert(std::pair<int, int>(eit->first, eit->second));

                        auto it = ce_counter.emplace(std::make_pair(eit->first, eid), 0);
                        it.first->second++;
                    }
                    int k = 0;
                    for (ordered_edge_set::iterator eit = tset.begin(); eit != tset.end(); eit++) {
                        ordered_edge_domains[eid][k] = eit->first;
                        ordered_edge_domains[eid][k + ordered_edge_domains_sizes[eid]] = eit->second;
#ifdef MDEBUG
                        std::cout << "IN eid: " << eid << ";k " << k << ": " << eit->first << "-" << eit->second << "\n";
#endif
                        ce_positions.insert({std::make_pair(eit->first, eid), k});

                        k++;
                    }
                }
            }
        }

        std::cout << "ORDERED EDGE SETS: done\n";

        int **f_domains = new int *[nof_sn];
        for (int si = 0; si < nof_sn; si++) {
            if (mama.edges_sizes[si] == 0) {
                int n = map_state_to_node[si];
                f_domains[si] = new int[domains_size[n]];
                int k = 0;
                for (sbitset::iterator IT = domains[n].first_ones(); IT != domains[n].end(); IT.next_ones()) {
                    f_domains[si][k] = IT.first;
                    k++;
                }
            }
        }

        int *candidateIT = new int[nof_sn];
        int *candidateITeid = new int[nof_sn];
        int *candidateITpnode = new int[nof_sn];
        int *candidateITpstate = new int[nof_sn];
        int *candidateITsize = new int[nof_sn];

        for (int i = 0; i < nof_sn; i++) {
            candidateIT[i] = -1;
        }

        int *solution = new int[nof_sn]; // indexed by state_id
        for (int i = 0; i < nof_sn; i++)
            solution[i] = -1;

        bool *matched = (bool *)calloc(rgraph.nof_nodes, sizeof(bool)); // indexed by node_id

        int psi = -1;
        int si = 0;
        int CandidateIndex = -1;
        int sip1;
        int eid, pnode, maxp, maxe, maxpcount;
        while (si != -1) {

#ifdef MDEBUG
            std::cout << "----------------------------------------\n";
            std::cout << "SI " << si << " - PATTERN " << map_state_to_node[si] << "\n";
#endif

            if (psi >= si) {
                matched[solution[si]] = false;
            }

            CandidateIndex = -1;

            if (candidateIT[si] == -1) {
                if (mama.edges_sizes[si] == 0) {
                    candidateIT[si] = -1;
                    candidateITsize[si] = domains_size[map_state_to_node[si]];
                } else if (mama.edges_sizes[si] == 1) {
                    int pstate;
                    eid = mama.edges[si][0].id;

                    if (mama.edges[si][0].source == si) {
                        pnode = solution[mama.edges[si][0].target];
                        pstate = mama.edges[si][0].target;
                    } else {
                        pnode = solution[mama.edges[si][0].source];
                        pstate = mama.edges[si][0].source;
                    }

                    candidateITpnode[si] = pnode;
                    candidateITeid[si] = eid;
                    candidateIT[si] = ce_positions[std::make_pair(pnode, eid)] - 1;
                    candidateITsize[si] = ordered_edge_domains_sizes[eid];
#ifdef MDEBUG
                    std::cout << "single parent: eid " << candidateITeid[si] << "; pnode " << pnode << "; pstate " << pstate << ";ce position " << candidateIT[si] << "\n";
#endif
                } else {
                    maxp = -1;
                    int pstate;
                    for (int j = 0; j < mama.edges_sizes[si]; j++) {
                        eid = mama.edges[si][j].id;

                        if (mama.edges[si][j].source == si) {
                            pnode = solution[mama.edges[si][j].target];
                        } else {
                            pnode = solution[mama.edges[si][j].source];
                        }

                        if ((maxp == -1) || (ce_counter[std::make_pair(pnode, eid)] > maxpcount)) {
                            maxp = pnode;
                            maxe = eid;
                            maxpcount = ce_counter[std::make_pair(pnode, eid)];

                            if (mama.edges[si][j].source == si) {
                                pstate = mama.edges[si][j].target;
                            } else {
                                pstate = mama.edges[si][j].source;
                            }
                        }
                    }
                    candidateITpnode[si] = maxp;
                    candidateITpstate[si] = pstate;
                    candidateITeid[si] = maxe;
                    candidateIT[si] = ce_positions[std::make_pair(maxp, maxe)] - 1;
                    candidateITsize[si] = ordered_edge_domains_sizes[maxe];

#ifdef MDEBUG
                    std::cout << "multiple parents: eid " << candidateITeid[si] << "; pnode " << maxp << "; pstate " << pstate << ";ce position " << candidateIT[si] << "\n";
#endif
                }
            }

            if (mama.edges_sizes[si] == 0) {
                candidateIT[si]++;
                while (candidateIT[si] < candidateITsize[si]) {
                    if (!matched[f_domains[si][candidateIT[si]]]) {

                        CandidateIndex = f_domains[si][candidateIT[si]];

                        solution[si] = CandidateIndex;
                        if (edgesCheck(si, CandidateIndex, solution, matched)) {
                            break;
                        } else {
                            CandidateIndex = -1;
                        }
                    }
                    candidateIT[si]++;
                }
                if (candidateIT[si] >= candidateITsize[si]) {
                    CandidateIndex = -1;
                }
            } else {
                candidateIT[si]++;
                while (candidateIT[si] < candidateITsize[si]) {
#ifdef MDEBUG
                    std::cout << "pCI " << candidateIT[si] << "; size " << candidateITsize[si] << "; eid " << candidateITeid[si] << " " << ordered_edge_domains[candidateITeid[si]][candidateIT[si]] << "-" << ordered_edge_domains[candidateITeid[si]][candidateIT[si] + candidateITsize[si]] << ":" << ordered_edge_domains[candidateITeid[si]][candidateIT[si] + candidateITsize[si]] << "; pnode " << candidateITpnode[si] << "\n";
#endif
                    if (ordered_edge_domains[candidateITeid[si]][candidateIT[si]] == candidateITpnode[si]) {
                        if (!matched[ordered_edge_domains[candidateITeid[si]][candidateIT[si] + candidateITsize[si]]]) {
                            CandidateIndex = ordered_edge_domains[candidateITeid[si]][candidateIT[si] + candidateITsize[si]];
                            solution[si] = CandidateIndex;

                            if (mama.edges_sizes[si] > 1) {

                                bool checked = true;
                                for (int me = 0; me < mama.edges_sizes[si]; me++) {

                                    if (edomains.domains[mama.edges[si][me].id].count(std::pair<int, int>(solution[mama.edges[si][me].source], solution[mama.edges[si][me].target])) == 0) {
                                        checked = false;
                                        break;
                                    }
                                }
                                if (checked)
                                    checked &= edgesCheck(si, CandidateIndex, solution, matched);

                                if (checked) {
                                    break;
                                } else {
                                    CandidateIndex = -1;
                                }
                            } else {
                                break;
                            }

                        }
#ifdef MDEBUG
                        else {
                            std::cout << "already matched\n";
                        }
#endif
                    } else {
#ifdef MDEBUG
                        std::cout << "end of candidates\n";
#endif
                        CandidateIndex = -1;
                        break;
                    }
                    candidateIT[si]++;
                }
                if (candidateIT[si] >= candidateITsize[si]) {
#ifdef MDEBUG
                    std::cout << "end of candidate list\n";
#endif
                    CandidateIndex = -1;
                }
            }

#ifdef MDEBUG
            std::cout << "CI " << CandidateIndex << "\n";
#endif

            if (CandidateIndex == -1) {
                candidateIT[si] = -1;

                psi = si;
                si--;
            } else {
                matchedcouples++;

                if (si == nof_sn - 1) {
#ifdef PRINT_MATCHES
                    matchListener.match(nof_sn, map_state_to_node, solution);
#endif
                    matchcount++;

                    psi = si;
#ifdef FIRST_MATCH_ONLY
                    si = -1;
//					return IF U WANT JUST AN INSTANCE
#endif
#ifdef FIRST_100k_MATCHES
                    if (matchcount >= 100000)
                        si = -1;
#endif
                } else {
                    matched[solution[si]] = true;
                    sip1 = si + 1;
                    psi = si;
                    si++;
                }
            }
        }
    };

    void SolveLeafs()
    {

        int NumOfQueryVertex = mama.nof_sn;
        int* map_node_to_state = mama.map_node_to_state; // indexed by node_id
        int* map_state_to_node = mama.map_state_to_node; // indexed by state_id

        using cand_ecount_t = std::unordered_map<std::pair<int, int>, int, hash_pair>; // (tnodeid,eid) -> count

        cand_ecount_t ce_counter;
        cand_ecount_t ce_positions;

#ifdef SOLVER_H_MDEBUG
        std::cout << "ORDERED EDGE SETS\n";
#endif

        using ordered_edge_set = std::set<std::pair<int, int>>;

        std::cout << "ORDERED EDGE SETS...\n";

        auto ordered_edge_domains = new int*[edomains.nof_pattern_edges];
        auto ordered_edge_domains_sizes = new int[edomains.nof_pattern_edges];

        for (int eid = 0; eid < edomains.nof_pattern_edges; eid++)
        {
            unordered_edge_set* eset = &(edomains.domains[eid]);
            ordered_edge_domains[eid] = new int[eset->size() * 2];
            ordered_edge_domains_sizes[eid] = eset->size();
        }
        
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            for (int j = 0; j < mama.edges_sizes[i]; j++)
            {
                if (mama.edges[i][j].source == i)
                {
                    int eid = mama.edges[i][j].id;
                    ordered_edge_set tset;
                    unordered_edge_set* eset = &(edomains.domains[eid]);
                    for (auto eit = eset->begin(); eit != eset->end(); eit++)
                    {
                        tset.insert(std::pair<int, int>(eit->second, eit->first));

                        auto it = ce_counter.emplace(std::make_pair(eit->second, eid), 0);
                        it.first->second++;
                    }
                    int k = 0;
                    for (auto eit = tset.begin(); eit != tset.end(); eit++)
                    {
                        ordered_edge_domains[eid][k] = eit->first;
                        ordered_edge_domains[eid][k + ordered_edge_domains_sizes[eid]] = eit->second;
#ifdef SOLVER_H_MDEBUG
                        std::cout << "OUT eid: " << eid << ";k " << k << ": " << eit->first << "-" << eit->second << "\n";
#endif

                        ce_positions.insert({std::make_pair(eit->first, eid), k});

                        k++;
                    }
                }
                else
                {
                    int eid = mama.edges[i][j].id;
                    ordered_edge_set tset;
                    unordered_edge_set* eset = &(edomains.domains[eid]);
                    for (auto eit = eset->begin(); eit != eset->end(); eit++)
                    {
                        tset.insert(std::pair<int, int>(eit->first, eit->second));

                        auto it = ce_counter.emplace(std::make_pair(eit->first, eid), 0);
                        it.first->second++;
                    }
                    int k = 0;
                    for (auto eit = tset.begin(); eit != tset.end(); eit++)
                    {
                        ordered_edge_domains[eid][k] = eit->first;
                        ordered_edge_domains[eid][k + ordered_edge_domains_sizes[eid]] = eit->second;
#ifdef SOLVER_H_MDEBUG
                        std::cout << "IN eid: " << eid << ";k " << k << ": " << eit->first << "-" << eit->second << "\n";
#endif
                        ce_positions.insert({std::make_pair(eit->first, eid), k});

                        k++;
                    }
                }
            }
        }

        std::cout << "ORDERED EDGE SETS: done\n";

        auto f_domains = new int*[NumOfQueryVertex];
        for (int si = 0; si < NumOfQueryVertex; si++)
        {
            if (mama.edges_sizes[si] == 0)
            {
                int n = map_state_to_node[si];
                f_domains[si] = new int[domains_size[n]];
                int k = 0;
                for (sbitset::iterator IT = domains[n].first_ones(); IT != domains[n].end(); IT.next_ones())
                {
                    f_domains[si][k] = IT.first;
                    k++;
                }
            }
        }

        auto candidateIT = new int[NumOfQueryVertex];

        auto candidateITeid = new int[NumOfQueryVertex];
        auto candidateITpnode = new int[NumOfQueryVertex];
        auto candidateITpstate = new int[NumOfQueryVertex];
        auto candidateITsize = new int[NumOfQueryVertex];

        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            candidateIT[i] = -1;
        }

        auto solution = new int[NumOfQueryVertex]; // indexed by state_id
        for (int i = 0; i < NumOfQueryVertex; i++)
        {
            solution[i] = -1;
        }

        auto matched = (bool*)calloc(rgraph.nof_nodes, sizeof(bool)); // indexed by node_id

        int psi = -1;
        int StateIndex = 0;
        int CandidateIndex = -1;
        int StateIndexPlusOne;
        int eid, pnode, maxp, maxe, maxpcount;

#ifdef SOLVER_H_MDEBUG
        std::cout << "nof sn " << NumOfQueryVertex << "; nof leafs " << mama.nof_leafs << "; si " << StateIndex << "; to combinatorial escape " << (NumOfQueryVertex - mama.nof_leafs) << " \n";
#endif
        while (StateIndex != -1)
        {

#ifdef SOLVER_H_MDEBUG
            std::cout << "----------------------------------------\n";
            std::cout << "SI " << StateIndex << " - PATTERN " << map_state_to_node[StateIndex] << "\n";
#endif

            if (StateIndex == NumOfQueryVertex - mama.nof_leafs)
            {
                // we have leafs
                //Hyaniner: It doesn't seem necessary to check the calculation counting here.

                auto leaf_domains = new std::set<int>[mama.nof_leafs];

                long nof_leaf_solutions = 1;

                for (int l = 0; l < mama.nof_leafs; l++)
                {
                    int pstate;
                    eid = mama.edges[StateIndex + l][0].id;

                    if (mama.edges[StateIndex + l][0].source == StateIndex + l)
                    {
                        pnode = solution[mama.edges[StateIndex + l][0].target];
                        pstate = mama.edges[StateIndex + l][0].target;
                    }
                    else
                    {
                        pnode = solution[mama.edges[StateIndex + l][0].source];
                        pstate = mama.edges[StateIndex + l][0].source;
                    }

                    candidateITpnode[StateIndex + l] = pnode;
                    candidateITeid[StateIndex + l] = eid;
                    candidateIT[StateIndex + l] = ce_positions[std::make_pair(pnode, eid)] - 1;
                    candidateITsize[StateIndex + l] = ordered_edge_domains_sizes[eid];

                    candidateIT[StateIndex + l]++;
                    while (candidateIT[StateIndex + l] < candidateITsize[StateIndex + l])
                    {
                        if (ordered_edge_domains[candidateITeid[StateIndex + l]][candidateIT[StateIndex + l]] == candidateITpnode[StateIndex + l])
                        {
                            if (!matched[ordered_edge_domains[candidateITeid[StateIndex + l]][candidateIT[StateIndex + l] + candidateITsize[StateIndex + l]]])
                            {
                                CandidateIndex = ordered_edge_domains[candidateITeid[StateIndex + l]][candidateIT[StateIndex + l] + candidateITsize[StateIndex + l]];
                                leaf_domains[l].insert(CandidateIndex);
                            }
                        }
                        else
                        {
                            break;
                        }
                        candidateIT[StateIndex + l]++;
                    }

                    nof_leaf_solutions *= leaf_domains[l].size();
                }

                matchcount += nof_leaf_solutions;

#ifdef PRINT_MATCHES
                matchListener.match_multiple(NumOfQueryVertex, map_state_to_node, solution, StateIndex, leaf_domains);
#endif

                delete[] leaf_domains;

                psi = StateIndex;
                StateIndex--;

#ifdef FIRST_MATCH_ONLY
                StateIndex = -1;
                //					return IF U WANT JUST AN INSTANCE
#endif
#ifdef FIRST_100k_MATCHES
                if (matchcount >= 100000)
                    StateIndex = -1;
#endif
            }
            else
            {

                if (psi >= StateIndex)
                {
                    matched[solution[StateIndex]] = false;
                }

                CandidateIndex = -1;

                if (candidateIT[StateIndex] == -1)
                {
                    if (mama.edges_sizes[StateIndex] == 0)
                    {
                        candidateIT[StateIndex] = -1;
                        candidateITsize[StateIndex] = domains_size[map_state_to_node[StateIndex]];
                    }
                    else if (mama.edges_sizes[StateIndex] == 1)
                    {
                        int pstate;
                        eid = mama.edges[StateIndex][0].id;

                        if (mama.edges[StateIndex][0].source == StateIndex)
                        {
                            pnode = solution[mama.edges[StateIndex][0].target];
                            pstate = mama.edges[StateIndex][0].target;
                        }
                        else
                        {
                            pnode = solution[mama.edges[StateIndex][0].source];
                            pstate = mama.edges[StateIndex][0].source;
                        }

                        candidateITpnode[StateIndex] = pnode;
                        candidateITeid[StateIndex] = eid;
                        candidateIT[StateIndex] = ce_positions[std::make_pair(pnode, eid)] - 1;
                        candidateITsize[StateIndex] = ordered_edge_domains_sizes[eid];
#ifdef SOLVER_H_MDEBUG
                        std::cout << "single parent: eid " << candidateITeid[StateIndex] << "; pnode " << pnode << "; pstate " << pstate << ";ce position " << candidateIT[StateIndex] << "\n";
#endif
                    }
                    else
                    {
                        maxp = -1;
                        int pstate;
                        for (int j = 0; j < mama.edges_sizes[StateIndex]; j++)
                        {
                            eid = mama.edges[StateIndex][j].id;
                            if (mama.edges[StateIndex][j].source == StateIndex)
                            {
                                pnode = solution[mama.edges[StateIndex][j].target];
                            }
                            else
                            {
                                pnode = solution[mama.edges[StateIndex][j].source];
                            }

                            if ((maxp == -1) || (ce_counter[std::make_pair(pnode, eid)] > maxpcount))
                            {
                                maxp = pnode;
                                maxe = eid;
                                maxpcount = ce_counter[std::make_pair(pnode, eid)];

                                if (mama.edges[StateIndex][j].source == StateIndex)
                                {
                                    pstate = mama.edges[StateIndex][j].target;
                                }
                                else
                                {
                                    pstate = mama.edges[StateIndex][j].source;
                                }
                            }
                        }
                        candidateITpnode[StateIndex] = maxp;
                        candidateITpstate[StateIndex] = pstate;
                        candidateITeid[StateIndex] = maxe;
                        candidateIT[StateIndex] = ce_positions[std::make_pair(maxp, maxe)] - 1;
                        candidateITsize[StateIndex] = ordered_edge_domains_sizes[maxe];

#ifdef SOLVER_H_MDEBUG
                        std::cout << "multiple parents: eid " << candidateITeid[StateIndex] << "; pnode " << maxp << "; pstate " << pstate << ";ce position " << candidateIT[StateIndex] << "\n";
#endif
                    }
                }

                if (mama.edges_sizes[StateIndex] == 0)
                {
                    candidateIT[StateIndex]++;
                    while (candidateIT[StateIndex] < candidateITsize[StateIndex])
                    {
                        //Hyaniner: Calculation counting point A
                        if (!matched[f_domains[StateIndex][candidateIT[StateIndex]]])
                        {

                            CandidateIndex = f_domains[StateIndex][candidateIT[StateIndex]];

                            solution[StateIndex] = CandidateIndex;
                            if (edgesCheck(StateIndex, CandidateIndex, solution, matched))
                            {
                                break;
                            }
                            else
                            {
                                CandidateIndex = -1;
                            }
                        }
                        candidateIT[StateIndex]++;
                    }
                    if (candidateIT[StateIndex] >= candidateITsize[StateIndex])
                    {
                        CandidateIndex = -1;
                    }
                }
                else
                {
                    candidateIT[StateIndex]++;
                    while (candidateIT[StateIndex] < candidateITsize[StateIndex])
                    {
#ifdef SOLVER_H_MDEBUG
                        std::cout << "pCI " << candidateIT[StateIndex] << "; size " << candidateITsize[StateIndex] << "; eid " << candidateITeid[StateIndex] << " " << ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex]]
                            << "-" << ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex] + candidateITsize[StateIndex]] << ":" << ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex] +
                                candidateITsize[StateIndex]] << "; pnode " << candidateITpnode[StateIndex] << "\n";
#endif
                        if (ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex]] == candidateITpnode[StateIndex])
                        {
                            if (!matched[ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex] + candidateITsize[StateIndex]]])
                            {
                                CandidateIndex = ordered_edge_domains[candidateITeid[StateIndex]][candidateIT[StateIndex] + candidateITsize[StateIndex]];
                                solution[StateIndex] = CandidateIndex;

                                if (mama.edges_sizes[StateIndex] > 1)
                                {

                                    bool checked = true;
                                    for (int me = 0; me < mama.edges_sizes[StateIndex]; me++)
                                    {
                                        //Hyaniner: Calculation counting point B
                                        if (edomains.domains[mama.edges[StateIndex][me].id].count(std::pair<int, int>(solution[mama.edges[StateIndex][me].source], solution[mama.edges[StateIndex][me].target])) == 0)
                                        {
                                            checked = false;
                                            break;
                                        }
                                    }
                                    if (checked)
                                    {
                                        checked &= edgesCheck(StateIndex, CandidateIndex, solution, matched);
                                    }

                                    if (checked)
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        CandidateIndex = -1;
                                    }
                                }
                                else
                                {
                                    break;
                                }

                            }
#ifdef SOLVER_H_MDEBUG
                            else {
                                std::cout << "already matched\n";
                            }
#endif
                        }
                        else
                        {
#ifdef SOLVER_H_MDEBUG
                            std::cout << "end of candidates\n";
#endif
                            CandidateIndex = -1;
                            break;
                        }
                        candidateIT[StateIndex]++;
                    }
                    if (candidateIT[StateIndex] >= candidateITsize[StateIndex])
                    {
#ifdef SOLVER_H_MDEBUG
                        std::cout << "end of candidate list\n";
#endif
                        CandidateIndex = -1;
                    }
                }

#ifdef SOLVER_H_MDEBUG
                std::cout << "CI " << CandidateIndex << "\n";
#endif

                if (CandidateIndex == -1)
                {
                    //Hyaniner:Backtracking
                    candidateIT[StateIndex] = -1;

                    psi = StateIndex;
                    StateIndex--;
                }
                else
                {
                    //Hyaniner:matched pair found.
                    matchedcouples++;

                    if (StateIndex == NumOfQueryVertex - 1)
                    {
                        // there are no leafs
#ifdef PRINT_MATCHES
                        matchListener.match(NumOfQueryVertex, map_state_to_node, solution);
#endif
                        matchcount++;
                        psi = StateIndex;

#ifdef FIRST_MATCH_ONLY
                        StateIndex = -1;
#endif
#ifdef FIRST_100k_MATCHES
                        if (matchcount >= 100000)
                            StateIndex = -1;
#endif

                    }
                    else
                    {
                        matched[solution[StateIndex]] = true;
                        StateIndexPlusOne = StateIndex + 1;
                        psi = StateIndex;
                        StateIndex++;
                    }
                }
            }
        }
    };

    virtual bool edgesCheck(int si, int ci, int *solution, bool *matched) = 0;
};

} // namespace rilib

#endif /* SOLVER_H_ */
