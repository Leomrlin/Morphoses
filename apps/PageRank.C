

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

#include "../core/common/utils.h"
#include "../core/morphoses/MorphosesEngine.h"
#include "../core/main.h"
#include <math.h>

// Morphoses
// template <class AggregationValueType, class GlobalInfoType>
// inline void initializeAggregationValue(const uintV &v,AggregationValueType &v_aggregation_value,const GlobalInfoType &global_info);
// template <class VertexValueType, class GlobalInfoType>
// inline bool initializeVertexValue(const uintV &v,VertexValueType &v_vertex_value,const GlobalInfoType &global_info);//需要设置是否激活
// template <class VertexValueType> 
// inline bool isVertexValueEqual(const VertexValueType &a, const VertexValueType &b);
// template <class AggregationValueType, class VertexValueType, class GlobalInfoType>
// inline bool computeFunction(const uintV &v,
//                             int iter,
//                             const AggregationValueType &aggregation_value,
//                             const VertexValueType &vertex_value_prev,
//                             VertexValueType &vertex_value_curr,
//                             GlobalInfoType &global_info);//需要返回是否继续激活
// template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
// inline void edgeFunction(const uintV &u, const uintV &v, int iter,
//                          const EdgeDataType &edge_weight,
//                          const VertexValueType &u_value,
//                          AggregationValueType &aggregate_value,
//                          GlobalInfoType &global_info);

// ======================================================================
// PAGERANKINFO
// ======================================================================
template <class vertex>
class PageRankInfo {
public:
  // Should I just use vectors for this?
  graph<vertex> *my_graph;
  uintV n;
  double epsilon;
  double damping;
  long *out_degrees;

  PageRankInfo() : my_graph(nullptr), n(0), epsilon(0), damping(0), out_degrees(nullptr) {}

  PageRankInfo(graph<vertex> *_my_graph, uintV _n, double _epsilon, double _damping)
      : my_graph(_my_graph), n(_n), epsilon(_epsilon), damping(_damping) {
    if (n > 0) {
      out_degrees = newA(long, n);
      parallel_for(uintV i = 0; i < n; i++) { out_degrees[i] = 0; }
    }
  }

  void init(){
    parallel_for(uintV i = 0; i < n; i++) {
      // out_degrees[i] = my_graph->V[i].getOutDegree();
      out_degrees[i] = my_graph->getOutDegree(i, 0);
    }    
  }

  void copy(const PageRankInfo &object) {
    if (object.n > n) {
      if (n == 0) {
        n = object.n;
        out_degrees = newA(long, n);
      } else {
        // realloc
        n = object.n;
        out_degrees = renewA(long, out_degrees, n);
      }
    }
    long min_n = std::min(object.n, n);
    parallel_for(uintV i = 0; i < min_n; i++) {
      out_degrees[i] = object.out_degrees[i];
    }
    epsilon = object.epsilon;
    damping = object.damping;
  }

  ~PageRankInfo() {
    if (n > 0)
      deleteA(out_degrees);
  }

  void processUpdates(edgeArray &edge_additions, edgeArray &edge_deletions) {
    // Increase out_degrees array size
    if (edge_additions.maxVertex >= n) {
      uintV n_old = n;
      n = edge_additions.maxVertex + 1;
      out_degrees = renewA(long, out_degrees, n);
      parallel_for(uintV i = n_old; i < n; i++) { out_degrees[i] = 0; }
    }

    parallel_for(long i = 0; i < edge_additions.size; i++) {
      uintV source = edge_additions.E[i].source;
      uintV destination = edge_additions.E[i].destination;
      writeAdd(&out_degrees[source], (long)1);
    }
    parallel_for(long i = 0; i < edge_deletions.size; i++) {
      uintV source = edge_deletions.E[i].source;
      uintV destination = edge_deletions.E[i].destination;
      writeAdd(&out_degrees[source], (long)-1);
    }
  }

  void cleanup() {}

  inline long getOutDegree(const uintV &v, const long version) {
    return my_graph->getOutDegree(v,version);
  }
};

struct MyPRVertexValue {
  double val;
  long outDegree;
};
// ======================================================================
// AGGREGATEVALUE AND VERTEXVALUE INITIALIZATION
// ======================================================================
double initial_aggregation_value = 0;
double initial_vertex_value = 1;

template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const uintV &v,AggregationValueType &v_aggregation_value,  long version, GlobalInfoType &global_info){
  v_aggregation_value = initial_aggregation_value;
}
template <class VertexValueType, class GlobalInfoType>
inline bool initializeVertexValue(const uintV &v,VertexValueType &v_vertex_value,  long version, GlobalInfoType &global_info){  //需要设置是否激活
  v_vertex_value.val = initial_vertex_value;
  v_vertex_value.outDegree = global_info.getOutDegree(v, version);
  return true;
}
template <class VertexValueType, class GlobalInfoType> 
inline bool isVertexValueEqual(const VertexValueType &a, const VertexValueType &b,  long version, const GlobalInfoType &global_info){
  return ((a.val-b.val)>-global_info.epsilon) && ((a.val-b.val)<global_info.epsilon) && a.outDegree == b.outDegree;
}

// ======================================================================
// VERTEX COMPUTE FUNCTION AND DETERMINE END OF COMPUTATION
// ======================================================================
template <class AggregationValueType, class VertexValueType, class GlobalInfoType>
inline bool computeFunction(const uintV &v,
                            long version,
                            int iter,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_prev,
                            VertexValueType &vertex_value_curr,
                            GlobalInfoType &global_info){//需要返回是否继续激活
  vertex_value_curr.val =
      (1 - global_info.damping) + (global_info.damping * aggregation_value);
  vertex_value_curr.outDegree = global_info.getOutDegree(v, version);
  return true;
}
// ======================================================================
// EDGE FUNCTIONS
// ======================================================================
template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline void edgeFunction(const uintV &u, const uintV &v, long version, int iter,
                        const EdgeDataType &edge_weight,
                        const VertexValueType &u_value,
                        AggregationValueType &aggregate_value,
                        GlobalInfoType &global_info){
  aggregate_value += (global_info.getOutDegree(u, version) != 0) ? ((u_value.val) / (global_info.getOutDegree(u, version))) : 0;
}


template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline bool incrementalAggregationValue(const uintV &u, const uintV &v, long version, int iter,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value_old,  bool u_active_old,
                         const VertexValueType &u_value_new,  bool u_active_new,
                         AggregationValueType &aggregate_value,
                         GlobalInfoType &global_info){

  // cout << "incre u: " <<u <<" v: " << v 
  // << " u_value_old: " <<u_value_old.val <<"\t"  <<u_value_old.outDegree
  // << " u_value_new: " <<u_value_new.val <<"\t"  <<u_value_new.outDegree
  // << " old_aggregate_value: " << aggregate_value << "\n";
  double oldVal = (u_value_old.outDegree != 0) ? ((u_value_old.val) / (u_value_old.outDegree)) : 0;
  double newVal = (u_value_new.outDegree != 0) ? ((u_value_new.val) / (u_value_new.outDegree)) : 0;
  aggregate_value += (newVal - oldVal);
  return true;
  // cout << "incre u: " <<u <<" v: " << v 
  // << " oldVal " <<oldVal <<" newVal: " << newVal
  // << " new_aggregate_value: " <<aggregate_value << "\n";
};

inline bool useIncrementalVertexValueRecompute(){
  return true;
}

inline bool theSameVertexVlaueTheSameMessage(){
  return false;
}
// ======================================================================
// HELPER FUNCTIONS
// ======================================================================
// template <class AggregationValueType, class VertexValueType,
//           class GlobalInfoType>
// void printHistory(const uintV &v, AggregationValueType **agg_values,
//                   VertexValueType **vertex_values, GlobalInfoType &info,
//                   int history_iterations) {
//   for (int iter = 0; iter < history_iterations; iter++) {
//     cout << iter << "," << agg_values[iter][v] << "," << vertex_values[iter][v]
//          << "\n";
//   }
// }

// template <class GlobalInfoType>
// void printAdditionalData(ofstream &output_file, const uintV &v,
//                          GlobalInfoType &info) {}

// ======================================================================
// COMPUTE FUNCTION
// ======================================================================
template <class vertex> void compute(graph<vertex> &G, commandLine config) {
  uintV n = G.n;
  int max_iters = config.getOptionLongValue("-maxIters", 10);
  double epsilon = config.getOptionDoubleValue("-epsilon", 0.01d);
  max_iters += 1;
  double damping = 0.85;

  PageRankInfo<vertex> global_info(&G, n, epsilon, damping);

  cout << "Initializing engine ....\n";
  MorphosesEngine<vertex, double, MyPRVertexValue, PageRankInfo<vertex>> engine(
      G, max_iters, global_info, config);
  engine.init();
  cout << "Finished initializing engine\n";

  engine.run();
}

