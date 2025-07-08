/*
 * Graph.h
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

#ifndef GRAPH_H_
#define GRAPH_H_

namespace rilib
{
class FRiGraph
{
public:
    int NumOfVertex;

    void** VertexAttributes;

    int* OutAdjSizes;
    int* InAdjSizes;

    int** OutAdjList;
    int** InAdjList;
    void*** OutAdjAttributes;

    FRiGraph()
    {

        NumOfVertex = 0;
        VertexAttributes = NULL;
        OutAdjSizes = NULL;
        InAdjSizes = NULL;
        OutAdjList = NULL;
        InAdjList = NULL;
        OutAdjAttributes = NULL;
    }
};
} // namespace rilib

#endif /* REFERENCEGRAPH_H_ */
