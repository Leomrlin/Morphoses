

// Copyright (c) 2020 Mugilan Mariappan, Joanna Che and Keval Vora.
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifdef EDGEDATA
// NOTE: The edge data type header file should then be included as the first header
// file at the top of the user program.
#include "SSSP_edgeData.h"
#endif

#include "../core/common/utils.h"
#include "../core/morphoses/MorphosesEngine.h"
#include "../core/main.h"
#include <math.h>

#define MAX_DISTANCE 65535

// ======================================================================
// SSSPINFO
// ======================================================================
class SsspInfo {
public:
  uintV source_vertex;
  long weight_cap;

  SsspInfo() : source_vertex(0), weight_cap(3) {}

  SsspInfo(uintV _source_vertex, long _weight_cap)
      : source_vertex(_source_vertex), weight_cap(_weight_cap) {}


  uint16_t getWeight(uintV u, uintV v) {
    return (uint16_t)((u + v) % weight_cap + 1);
  }

  void copy(const SsspInfo &object) {
    source_vertex = object.source_vertex;
    weight_cap = object.weight_cap;
  }
  void init() {}

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {}

  void cleanup() {}
};

// ======================================================================
// VERTEXVALUE INITIALIZATION
// ======================================================================
template <class VertexValueType, class GlobalInfoType>
inline bool initializeVertexValue(const uintV &v,VertexValueType &v_vertex_value, GlobalInfoType &global_info){  //需要设置是否激活
  if (v != global_info.source_vertex) {
    v_vertex_value = MAX_DISTANCE;
  } else {
    v_vertex_value = 0;
  }
  return true;
}

template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const uintV &v,AggregationValueType &v_aggregation_value, GlobalInfoType &global_info){
  v_aggregation_value = MAX_DISTANCE;
}

template <class VertexValueType, class GlobalInfoType>
inline bool isVertexValueEqual(const VertexValueType &a, const VertexValueType &b, const GlobalInfoType &global_info){
  return (a == b);
}


// ======================================================================
// EDGE FUNCTION
// ======================================================================
template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline void edgeFunction(const uintV &u, const uintV &v, int iter,
                        const EdgeDataType &edge_weight,
                        const VertexValueType &u_value,
                        AggregationValueType &aggregate_value,
                        GlobalInfoType &global_info){
  if (u_value == MAX_DISTANCE) {
    return;
  } 
  else {
    AggregationValueType v_value = u_value + global_info.getWeight(u, v);
    aggregate_value = (v_value < aggregate_value) ? v_value : aggregate_value;
  }
}

template <class AggregationValueType, class VertexValueType, class GlobalInfoType>
inline bool computeFunction(const uintV &v,
                            int iter,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_prev,
                            VertexValueType &vertex_value_curr,
                            GlobalInfoType &global_info){//需要返回是否继续激活
  if(aggregation_value < vertex_value_prev){
    vertex_value_curr = aggregation_value;
    return true;
  }
  else{
    vertex_value_curr = vertex_value_prev;
    return false;
  }
}

template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline bool incrementalAggregationValue(const uintV &u, const uintV &v, int iter,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value_old,  bool u_active_old,
                         const VertexValueType &u_value_new,  bool u_active_new,
                         AggregationValueType &aggregate_value,
                         GlobalInfoType &global_info){
  if(u_active_old){//撤销旧值影响
    AggregationValueType v_value = u_value_old + global_info.getWeight(u, v);
    if(v_value <= aggregate_value) return false;
  }

  if(u_active_new){//添加新值影响
    AggregationValueType v_value = u_value_new + global_info.getWeight(u, v);
    aggregate_value = (v_value < aggregate_value) ? v_value : aggregate_value;
  }
  return true;
};

inline bool useIncrementalVertexValueRecompute(){
  return true;
}
inline bool theSameVertexVlaueTheSameMessage(){
  return true;
}
// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  long n = G.n;
  int source_vertex = config.getOptionLongValue("-source", 1);
  int weight_cap = config.getOptionLongValue("-weight_cap", 5);
  SsspInfo global_info(source_vertex, weight_cap);

  cout << "Initializing engine ....\n";
  //KickStarterEngine<vertex, uint16_t, SsspInfo> engine(G, global_info, config);
  MorphosesEngine<vertex, uint16_t, uint16_t, SsspInfo> engine(G, 200, global_info, config);
  engine.init();
  cout << "Finished initializing engine\n";
  engine.run();
}

