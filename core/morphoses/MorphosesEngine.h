
//MorphosesEngine.h
#ifndef MORPHOSES_ENGINE_H
#define MORPHOSES_ENGINE_H

#include "../graph/edgeDataType.h"
#include "../common/utils.h"
#include "../common/ThreadPool.h"
#include "../graph/versionGraph.h"
#include "../graph/versionDelta.h"
#include "../graph/vertexSubset.h"
#include <limits>
#include "ingestorMorphoses.h"
#include <vector>
#include <queue>
#include <thread>
#include <utility>

#define TEST_THREAD_NUM 1

#define CONVERGED_FLAG 0
#define ITERATION_FLAG 1

#define INVALID_FLAG 0
#define TODO_FLAG 1
#define RUNNING_FLAG 2
#define DONE_FLAG 3
#define RELEASED_FLAG 4

#define CONTINUE_FLAG 5
#define READY_FLAG 6
#define FREE_FLAG 7

#define AGGREGATION_MASK 1
#define AGGREGATION_INVALID 0
#define AGGREGATION_VALID 1
#define ACTIVE_MASK 2
#define IS_ACTIVE_FLAG 2
#define NOT_ACTIVE_FLAG 0

enum UpdateType { edge_addition_enum, edge_deletion_enum };

#ifdef EDGEDATA
#else
struct EmptyEdgeData {};
typedef EmptyEdgeData EdgeData;

EdgeData emptyEdgeData;
#endif

struct task{
  int flag;
  long version;
  long iteration;
  edgeArray addition;
  edgeArray deletion;
};

//主从线程交互
struct Pos{
  int row_pos; int col_pos;
  Pos(int _row_pos, int _col_pos):row_pos(_row_pos), col_pos(_col_pos){}
  Pos(const Pos & p):row_pos(p.row_pos), col_pos(p.col_pos){}
};
SafeQueue<Pos> execResultQueue;

// ======================================================================
// 应用实现的接口方法
// ======================================================================
//Morphoses
template <class AggregationValueType, class GlobalInfoType>
inline void initializeAggregationValue(const uintV &v,AggregationValueType &v_aggregation_value, long version, GlobalInfoType &global_info);
template <class VertexValueType, class GlobalInfoType>
inline bool initializeVertexValue(const uintV &v,VertexValueType &v_vertex_value, long version, GlobalInfoType &global_info);//需要设置是否激活
template <class VertexValueType, class GlobalInfoType>
inline bool isVertexValueEqual(const VertexValueType &a, const VertexValueType &b, long version, const GlobalInfoType &global_info);
template <class AggregationValueType, class VertexValueType, class GlobalInfoType>
inline bool computeFunction(const uintV &v,
                            long version,
                            int iter,
                            const AggregationValueType &aggregation_value,
                            const VertexValueType &vertex_value_prev,
                            VertexValueType &vertex_value_curr,
                            GlobalInfoType &global_info);//需要返回是否继续激活
template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline void edgeFunction(const uintV &u, const uintV &v, long version, int iter,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value,
                         AggregationValueType &aggregate_value,
                         GlobalInfoType &global_info);
template <class AggregationValueType, class VertexValueType, class EdgeDataType, class GlobalInfoType>
inline bool incrementalAggregationValue(const uintV &u, const uintV &v, long version, int iter,
                         const EdgeDataType &edge_weight,
                         const VertexValueType &u_value_old,  bool u_active_old,
                         const VertexValueType &u_value_new,  bool u_active_new,
                         AggregationValueType &aggregate_value,
                         GlobalInfoType &global_info);
inline bool theSameVertexVlaueTheSameMessage();
inline bool useIncrementalVertexValueRecompute();
// ======================================================================
// CLASS VSS：VertexStateType
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class VertexStateSet{

public:
  graph<vertex> & my_graph;
  uintEE batch_versionId;
  int iter;
  int state;
  uintV active_num;
  GlobalInfoType & global_info;
	VertexStateSet * baseVSS;					//作为基的VSS
	VertexStateSet * preIterationVSS;			//上轮迭代的VSS
	//基锚定点
	AggregationValueType * av;
	VertexValueType * vv;
  bool * active;
	uintEE * versionId;
	unsigned int * flags;
  //capacity
	uintV n;
  uintV c;
  uintV valid_n;
  uintV transactionAffectedVertexNum;

  bool theSameVertexVlaueTheSameMessageSwitch;
  bool useIncrementalVertexValueRecomputeSwitch;

  edgeArray & edge_additions;
  edgeArray & edge_deletions;
	vector<uintV> * affectedVertexSet;			//受到影响点集合锚定点
  vector<VertexValueType> * affectedVertexSetHistoryValue;			//受到影响点上一帧的状态
  vector<bool> * affectedVertexSetActive;			//受到影响点在两帧是否有活动

  vector<uintV> * nextItr_affectedVertexSet;			//受到影响点集合锚定点
  vector<VertexValueType> * nextItr_affectedVertexSetHistoryValue;			//受到影响点上一帧的状态
  vector<bool> * nextItr_affectedVertexSetActive;			//受到影响点在两帧是否有活动

	vector<uintV> * activeVertexList;			//本轮迭代激活点列表

  VertexStateSet(graph<vertex> &_my_graph, uintEE _batch_versionId, int _iter, uintV _n, GlobalInfoType &_global_info, 
      VertexStateSet * _baseVSS, VertexStateSet * _preIterationVSS, edgeArray & _edge_additions, edgeArray & _edge_deletions)
      : my_graph(_my_graph), batch_versionId(_batch_versionId), iter(_iter), state(ITERATION_FLAG), active_num(0), n(_n), c(0), valid_n(0), 
      transactionAffectedVertexNum(0), baseVSS(_baseVSS), global_info(_global_info),
      preIterationVSS(_preIterationVSS), edge_additions(_edge_additions), edge_deletions(_edge_deletions),
      affectedVertexSet(nullptr), affectedVertexSetHistoryValue(nullptr), affectedVertexSetActive(nullptr),
      nextItr_affectedVertexSet(nullptr), nextItr_affectedVertexSetHistoryValue(nullptr), nextItr_affectedVertexSetActive(nullptr),
      activeVertexList(nullptr), av(nullptr), vv(nullptr), active(nullptr), versionId(nullptr), flags(nullptr) {
        theSameVertexVlaueTheSameMessageSwitch = theSameVertexVlaueTheSameMessage();
        useIncrementalVertexValueRecomputeSwitch = useIncrementalVertexValueRecompute();
  }

	void run(){
    // cout <<"run():" << " batch_versionId: "<< batch_versionId << " iter: "<< iter << " \n ";
    //获取基
		if(baseVSS != nullptr && baseVSS->vv != nullptr){	//向上一帧拉取基
      // cout << "Pull base from frontier."  << "\n";
      uintV new_n = n;
			av = baseVSS->av;
      vv = baseVSS->vv;
      active = baseVSS->active;
      versionId = baseVSS->versionId;
      flags = baseVSS->flags;
      baseVSS->av = nullptr;
      baseVSS->vv = nullptr;
      baseVSS->active = nullptr;
      baseVSS->versionId = nullptr;
      baseVSS->flags = nullptr;
      n = baseVSS->n;
      c = baseVSS->c;
      valid_n = baseVSS->n;
      active_num = baseVSS->active_num;
      adjustBase(new_n, preIterationVSS);
		}
		else if(preIterationVSS != nullptr && preIterationVSS->vv != nullptr){	//拷贝上轮迭代的基
      // cout << "Copy base from pre-iteration."  << "\n";
			copyAndCreateBase(*preIterationVSS);
      valid_n = 0;
		}
		else{ //找不到基重新构造
      // cout << "Cannot find base and re-construct."  << "\n";
      createNewBase();
      valid_n = 0;
    }
    // cout <<" valid_n: "<< valid_n << " \n ";
    if(valid_n == 0){
      if(preIterationVSS == nullptr || (preIterationVSS != nullptr && preIterationVSS->valid_n > 0)){ //首次转入传统计算模式产生全激活点的集合
        if(affectedVertexSet != nullptr)  affectedVertexSet->clear();
        else {
          affectedVertexSet = new vector<uintV>();
          affectedVertexSet->reserve(n);
        }
        for(long v = 0; v < n; v++) {
          if(active[v]) affectedVertexSet->push_back(v);
        }
        if(affectedVertexSet->empty()){
          state = CONVERGED_FLAG;
          return;
        }
      }
      if(preIterationVSS != nullptr && preIterationVSS->affectedVertexSet != nullptr){
        affectedVertexSet = preIterationVSS->affectedVertexSet;
        preIterationVSS->affectedVertexSet = nullptr;
        // cout << "Print before traditionalIteration()."  << "\n";
        //  printVSS(false);
        traditionalIteration();//无参考历史信息，传统迭代计算          
      }  
    }
    else if(preIterationVSS == nullptr){
      affectedVertexSet = new vector<uintV>();
      affectedVertexSet->reserve((n >> 10)+1);
      affectedVertexSetHistoryValue = new vector<VertexValueType>();
      affectedVertexSetHistoryValue->reserve((n >> 10)+1);
      affectedVertexSetActive = new vector<bool>();
      affectedVertexSetActive->reserve((n >> 10)+1);
      VertexValueType emptyVertexVlue;
      for(long i = 0; i < edge_additions.size; i++) {
        uintV source = edge_additions.E[i].source;
        uintV destination = edge_additions.E[i].destination;
        // cout <<"ADD "<< " source: "<< source << " destination: "<< destination<< "\n";
        if(versionId[source] != batch_versionId){
          versionId[source] = batch_versionId;
          affectedVertexSet->push_back(source);
          transactionAffectedVertexNum++;
          uintV v = source;
          initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
          flags[v] = AGGREGATION_VALID;
          VertexValueType vv_old = vv[v];
          bool active_old = active[v];
          active[v] = initializeVertexValue<VertexValueType, GlobalInfoType>(v, vv[v], batch_versionId, global_info);
          if(v < valid_n){
            affectedVertexSetHistoryValue->push_back(vv_old);
            affectedVertexSetActive->push_back(active_old);
          }
          else{
            affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
            affectedVertexSetActive->push_back(true);            
          }
        }
        if(versionId[destination] != batch_versionId){
          versionId[destination] = batch_versionId;
          affectedVertexSet->push_back(destination);
          transactionAffectedVertexNum++;
          uintV v = destination;
          initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
          flags[v] = AGGREGATION_VALID;
          VertexValueType vv_old = vv[v];
          bool active_old = active[v];
          active[v] = initializeVertexValue<VertexValueType, GlobalInfoType>(v, vv[v], batch_versionId, global_info);
          if(v < valid_n){
            affectedVertexSetHistoryValue->push_back(vv_old);
            affectedVertexSetActive->push_back(active_old);
          }
          else{
            affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
            affectedVertexSetActive->push_back(true);            
          }
        }
      }
      for(long i = 0; i < edge_deletions.size; i++) {
        uintV source = edge_deletions.E[i].source;
        uintV destination = edge_deletions.E[i].destination;
        // cout <<"DEL "<< " source: "<< source << " destination: "<< destination<< "\n";
        if(versionId[source] != batch_versionId){
          versionId[source] = batch_versionId;
          affectedVertexSet->push_back(source);
          transactionAffectedVertexNum++;
          uintV v = source;
          initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
          flags[v] = AGGREGATION_VALID;
          VertexValueType vv_old = vv[v];
          bool active_old = active[v];
          active[v] = initializeVertexValue<VertexValueType, GlobalInfoType>(v, vv[v], batch_versionId, global_info);
          if(v < valid_n){
            affectedVertexSetHistoryValue->push_back(vv_old);
            affectedVertexSetActive->push_back(active_old);
          }
          else{
            affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
            affectedVertexSetActive->push_back(true);            
          }
        }
        if(versionId[destination] != batch_versionId){
          versionId[destination] = batch_versionId;
          affectedVertexSet->push_back(destination);
          transactionAffectedVertexNum++;
          uintV v = destination;
          initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
          flags[v] = AGGREGATION_VALID;
          VertexValueType vv_old = vv[v];
          bool active_old = active[v];
          active[v] = initializeVertexValue<VertexValueType, GlobalInfoType>(v, vv[v], batch_versionId, global_info);
          if(v < valid_n){
            affectedVertexSetHistoryValue->push_back(vv_old);
            affectedVertexSetActive->push_back(active_old);
          }
          else{
            affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
            affectedVertexSetActive->push_back(true);            
          }
        }
      }
      if(affectedVertexSet->empty()){
        state = CONVERGED_FLAG;
        return;
      }
    }
		else if(preIterationVSS != nullptr && preIterationVSS->affectedVertexSet != nullptr){//从上轮迭代拉取受影响点集合
			affectedVertexSet = preIterationVSS->affectedVertexSet;
      transactionAffectedVertexNum = preIterationVSS->transactionAffectedVertexNum;
      affectedVertexSetHistoryValue = preIterationVSS->affectedVertexSetHistoryValue;
      affectedVertexSetActive = preIterationVSS->affectedVertexSetActive;
			preIterationVSS->affectedVertexSet = nullptr;	
      preIterationVSS->affectedVertexSetHistoryValue = nullptr;		
      preIterationVSS->affectedVertexSetActive = nullptr;	
      // cout << "Print before computeVertexStateSet()."  << "\n";
      // printVSS(false);
      computeVertexStateSet();          	
		}
    else{//未果表示事务无影响，直接判定为收敛
      cout << "Inaffected transaction."  << "\n";
      if(affectedVertexSet != nullptr)  affectedVertexSet->clear();
      else { affectedVertexSet = new vector<uintV>(); }
      state = CONVERGED_FLAG;
      return;
    }

    // cout <<"run() over:" << " batch_versionId: "<< batch_versionId << " iter: "<< iter << " active_num: "<< active_num  << " \n ";

    // if(affectedVertexSet == nullptr)
    //   cout << " affectedVertexSet: null" << "\n";
    // else{
    //   cout << " affectedVertexSet: size() " <<affectedVertexSet->size() << "\n";
    // }

    // cout << "Print After run()."  << "\n";
    // printVSS(true);

    if(active_num == 0){
      state = CONVERGED_FLAG;
      return;
    }    
	}

  void initVertexState(long start_index, long end_index, VertexStateSet * preIterationVSS) {
    if(preIterationVSS == nullptr)
    for(long v = start_index; v < end_index; v++) {
      initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
      active[v] = initializeVertexValue<VertexValueType, GlobalInfoType>(v, vv[v], batch_versionId, global_info);
      if(active[v]) active_num++;
      versionId[v] = 0;
      flags[v] = AGGREGATION_VALID;
    }
    else
    for(long v = start_index; v < end_index; v++) {
      vv[v] = preIterationVSS->vv[v];
      av[v] = preIterationVSS->av[v];
      active[v] = false;//拷贝需要全部重算
      if(active[v]) active_num++;
      versionId[v] = 0;
      flags[v] = AGGREGATION_INVALID;
    }    
  }

	void createNewBase(){
    if(vv != nullptr){
      deleteA(av);
      deleteA(vv);
      deleteA(active);
      deleteA(versionId);
      deleteA(flags);
      av = nullptr;
      vv = nullptr;
      active = nullptr;
      versionId = nullptr;
      flags = nullptr;
    }
    //分配内存
    c = n + (n >> 1) + 1;
    av = newA(AggregationValueType, c);
    vv = newA(VertexValueType, c);
    active = newA(bool, c);
    versionId = newA(uintEE, c);
    flags = newA(unsigned int, c);
    //初始化VSS
    initVertexState(0,n,nullptr);
	}

	void adjustBase(uintV new_n, VertexStateSet * preIterationVSS){
    if(new_n <= n){
      n = new_n;
    }
		else if(new_n < c){
			//将点状态列表填充至最大编号
      initVertexState(n,new_n,preIterationVSS);
      n = new_n;
		}
    else{
      //扩容重新分配内存,增量初始化
      c = new_n + (new_n >> 1) + 1;
      av = renewA(AggregationValueType, av, c);
      vv = renewA(VertexValueType, vv, c);
      active = renewA(bool, active, c);
      versionId = renewA(uintEE, versionId, c);
      flags = renewA(unsigned int, flags, c);
      initVertexState(n,new_n,preIterationVSS);
      n = new_n;
    }
	}

	void copyAndCreateBase(VertexStateSet & vss){
    if(vv != nullptr){
      deleteA(av);
      deleteA(vv);
      deleteA(active);
      deleteA(versionId);
      deleteA(flags);
      av = nullptr;
      vv = nullptr;
      active = nullptr;
      versionId = nullptr;
      flags = nullptr;
    }
    if(vss.n > 0){
      n = vss.n;
      c = n + (n >> 1) + 1;
      av = newA(AggregationValueType, c);
      vv = newA(VertexValueType, c);
      active = newA(bool, c);
      versionId = newA(uintEE, c);
      flags = newA(unsigned int, c);
      for(long v = 0; v < n; v++) {
        vv[v] = vss.vv[v];
        av[v] = vss.av[v];
        active[v] = false;//拷贝需要全部重算
        versionId[v] = 0;
        flags[v] = AGGREGATION_INVALID;
      }    
    }
	}
  

  void activeOutNeighbors(uintV v, long affectedSetIdx){
// 激活邻居的方法：有增量计算方法，则用历史状态计算消息，再用新状态计算消息，将差异值增量推出，聚合消息仍设置为有效，放入激活点列表等待处理
// 没有增量计算方法，将聚合消息设置为无效，放入激活点列表等待处理
    // cout << "activeOutNeighbors ... v: "<<v <<" affectedSetIdx: " << affectedSetIdx  << "\n";    
    VersionGraphIterator<vertex> scanOutEdgeItr(my_graph, v, batch_versionId, true);
    versionEdge * tmpRef = nullptr;
    if(useIncrementalVertexValueRecomputeSwitch){
      //有增量计算方法可用
      while(scanOutEdgeItr.hasNext() && scanOutEdgeItr.next(&tmpRef)){
        uintV u = tmpRef->vertex;
        if((flags[u] & AGGREGATION_MASK)== AGGREGATION_VALID){
          bool incrementalRecomputeSuccess = incrementalAggregationValue<AggregationValueType, VertexValueType, EdgeData, GlobalInfoType>(
            v,u,batch_versionId,iter,tmpRef->weight,
            (*affectedVertexSetHistoryValue)[affectedSetIdx], (*affectedVertexSetActive)[affectedSetIdx],
            preIterationVSS->vv[v], preIterationVSS->active[v],
            av[u], global_info
          );
          if(!incrementalRecomputeSuccess) flags[u] = (flags[u] & (~AGGREGATION_MASK));
        }
        if(preIterationVSS->active[v])
          flags[u] = flags[u] | IS_ACTIVE_FLAG;        
        if(versionId[u] != batch_versionId){
          versionId[u] = batch_versionId;
          activeVertexList->push_back(u);
        }
      }   
    }
    else{
      while(scanOutEdgeItr.hasNext() && scanOutEdgeItr.next(&tmpRef)){
        uintV u = tmpRef->vertex;
        flags[u] = (flags[u] & (~AGGREGATION_MASK));
        if(preIterationVSS->active[v])
          flags[u] = flags[u] | IS_ACTIVE_FLAG;        
        if(versionId[u] != batch_versionId){
          versionId[u] = batch_versionId;
          activeVertexList->push_back(u);
        }
      } 
    }
  }

  bool recomputeVertex(uintV v, bool forcePushIntoAffectedVertexSet){
    //激活状态下的重算
    if ((flags[v] & ACTIVE_MASK) == IS_ACTIVE_FLAG){
      flags[v] = flags[v] & (~IS_ACTIVE_FLAG);
      //本轮聚合消息无效，重新拉取聚合消息
      if ((flags[v] & AGGREGATION_MASK) == AGGREGATION_INVALID){
        VersionGraphIterator<vertex> scanInEdgeItr(my_graph, v, batch_versionId, false);
        versionEdge * tmpRef = nullptr;
        initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);
        while(scanInEdgeItr.hasNext() && scanInEdgeItr.next(&tmpRef)){
          uintV u = tmpRef->vertex;
          if(preIterationVSS->active[u]){
            edgeFunction<AggregationValueType, VertexValueType, EdgeData, GlobalInfoType>
                (u, v, batch_versionId, iter, tmpRef->weight, preIterationVSS->vv[u], av[v], global_info);
          }
        }
        flags[v] = (flags[v] | AGGREGATION_VALID);
      }
      //从聚合消息重算
      VertexValueType vertex_value_tmp;
      bool active_this_iter = computeFunction<AggregationValueType, VertexValueType, GlobalInfoType>
                                  (v, batch_versionId, iter, av[v], preIterationVSS->vv[v], vertex_value_tmp, global_info);
      //将新产生的点状态与上一帧的结果比较，更新受影响的集合
      bool sameVertexValue = v < valid_n 
                              && isVertexValueEqual<VertexValueType, GlobalInfoType>(vv[v], vertex_value_tmp, batch_versionId, global_info) 
                              && active_this_iter == active[v];
      if(forcePushIntoAffectedVertexSet || !sameVertexValue){
        nextItr_affectedVertexSet->push_back(v);
        if(v < valid_n){
          nextItr_affectedVertexSetHistoryValue->push_back(vv[v]);
          nextItr_affectedVertexSetActive->push_back(active[v]);
        }
        else{
          VertexValueType emptyVertexVlue;
          nextItr_affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
          nextItr_affectedVertexSetActive->push_back(true);            
        }
      }
      vv[v] = vertex_value_tmp;
      if(active[v] && !active_this_iter)  active_num--;
      else if(!active[v] && active_this_iter) active_num++;
      active[v] = active_this_iter;
      return sameVertexValue;
    }
    //非激活状态下的重算等于拷贝上轮迭代状态
    else{
      flags[v] = AGGREGATION_INVALID;
      VertexValueType vertex_value_tmp = preIterationVSS->vv[v];
      bool active_this_iter = false;
      //将新产生的点状态与上一帧的结果比较，更新受影响的集合
      bool sameVertexValue = v < valid_n 
                              && isVertexValueEqual<VertexValueType, GlobalInfoType>(vv[v], vertex_value_tmp, batch_versionId, global_info) 
                              && active_this_iter == active[v];
      if(forcePushIntoAffectedVertexSet || !sameVertexValue){
        nextItr_affectedVertexSet->push_back(v);
        if(v < valid_n){
          nextItr_affectedVertexSetHistoryValue->push_back(vv[v]);
          nextItr_affectedVertexSetActive->push_back(active[v]);
        }
        else{
          VertexValueType emptyVertexVlue;
          nextItr_affectedVertexSetHistoryValue->push_back(emptyVertexVlue);
          nextItr_affectedVertexSetActive->push_back(true);            
        }
      }
      vv[v] = vertex_value_tmp;
      if(active[v] && !active_this_iter)  active_num--;
      else if(!active[v] && active_this_iter) active_num++;
      active[v] = active_this_iter;
      return sameVertexValue;
    }
  }

	void computeVertexStateSet(){
    activeVertexList = new vector<uintV>();
    activeVertexList->reserve((n >> 10)+1);
    nextItr_affectedVertexSet = new vector<uintV>();
    nextItr_affectedVertexSet->reserve((n >> 10)+1);
    nextItr_affectedVertexSetHistoryValue = new vector<VertexValueType>();
    nextItr_affectedVertexSetHistoryValue->reserve((n >> 10)+1);
    nextItr_affectedVertexSetActive = new vector<bool>();
    nextItr_affectedVertexSetActive->reserve((n >> 10)+1);
    // cout << "computeVertexStateSet ..."  << "\n";
    // cout << " affectedVertexSet->size() : " <<affectedVertexSet->size() << "\n";
    // cout << " affectedVertexSetHistoryValue->size() : " <<affectedVertexSetHistoryValue->size() << "\n";
    // cout << " affectedVertexSetActive->size() : " <<affectedVertexSetActive->size() << "\n";
    // for(long i = 0; i < affectedVertexSet->size();i++){
    //     cout << (*affectedVertexSet)[i]<<"\t"<< (*affectedVertexSetHistoryValue)[i].val
    //      <<"\t"<< (*affectedVertexSetHistoryValue)[i].outDegree <<"\t"<< (*affectedVertexSetActive)[i]<< "\n " ;
    // }
    // cout << "\n";
    // cout <<affectedVertexSet->size() << "\n";
    /*  消息传播阶段  */
    //结构变化点激活点强制拉取重算
    for(long i = 0; i < transactionAffectedVertexNum; i++){
      uintV v = (*affectedVertexSet)[i];
      if(versionId[v] != batch_versionId){
        if(preIterationVSS->active[v])
        flags[v] = AGGREGATION_INVALID | IS_ACTIVE_FLAG;
        else
        flags[v] = AGGREGATION_INVALID | NOT_ACTIVE_FLAG;
        versionId[v] = batch_versionId;
        activeVertexList->push_back(v);
      }
    }
    //结构变化点消息传播
    for(long i = 0; i < transactionAffectedVertexNum; i++){
      uintV v = (*affectedVertexSet)[i];
      if(v < valid_n && ((*affectedVertexSetActive)[i] != (preIterationVSS->active[v]))){
        activeOutNeighbors(v, i);
      }
      else if(v < valid_n && ((*affectedVertexSetActive)[i]) && (preIterationVSS->active[v])){
        //状态相同消息相同的免于激活
        if(isVertexValueEqual<VertexValueType, GlobalInfoType>(preIterationVSS->vv[v], (*affectedVertexSetHistoryValue)[i], batch_versionId, global_info) 
              && theSameVertexVlaueTheSameMessageSwitch)  ;
        else  activeOutNeighbors(v, i);
      }
    }
    //状态差异点消息传播
    for(long i = transactionAffectedVertexNum; i < affectedVertexSet->size(); i++){
      uintV v = (*affectedVertexSet)[i];
      if(((*affectedVertexSetActive)[i]) || (preIterationVSS->active[v])){
        activeOutNeighbors(v, i);
      }
    }
    /*  重算阶段  */
    //重算结构差异点 
    bool forcePushIntoAffectedVertexSet = true;
    for(long i = 0; i < transactionAffectedVertexNum; i++){
      uintV v = (*activeVertexList)[i];
      recomputeVertex(v, forcePushIntoAffectedVertexSet);
    }
    //重算状态差异点中未被消息激活的部分
    for(long i = transactionAffectedVertexNum; i < affectedVertexSet->size(); i++){
      uintV v = (*affectedVertexSet)[i];
      if(versionId[v] != batch_versionId){
        versionId[v] = batch_versionId;
        if(preIterationVSS->active[v]){
          flags[v] = flags[v] | IS_ACTIVE_FLAG;
          activeVertexList->push_back(v);
        }
        else{
          recomputeVertex(v, false);
        }
      }
    }
    //重算激活的节点
    for(long i = transactionAffectedVertexNum; i < activeVertexList->size(); i++){
      uintV v = (*activeVertexList)[i];
      recomputeVertex(v, false);
    } 
    // cout <<"activeVertexList->size(): "<<activeVertexList->size() << "\n";
    // if(activeVertexList == nullptr)
    //   cout << " activeVertexList: null" << "\n";
    // else{
    //   cout << " activeVertexList: size() " <<activeVertexList->size() << "\n";
    // }
    //析构本轮影响点集合，交换下一轮影响点集合
    if(affectedVertexSet != nullptr)  delete affectedVertexSet;
    affectedVertexSet = nextItr_affectedVertexSet;
    nextItr_affectedVertexSet = nullptr;
    if(affectedVertexSetHistoryValue != nullptr)  delete affectedVertexSetHistoryValue;
    affectedVertexSetHistoryValue = nextItr_affectedVertexSetHistoryValue;
    nextItr_affectedVertexSetHistoryValue = nullptr;
    if(affectedVertexSetActive != nullptr)  delete affectedVertexSetActive;
    affectedVertexSetActive = nextItr_affectedVertexSetActive;
    nextItr_affectedVertexSetActive = nullptr;
    if(activeVertexList != nullptr)  delete activeVertexList;
    activeVertexList = nullptr;
    if(affectedVertexSet->empty()){
      state = CONVERGED_FLAG;
    }
	}

  void traditionalIteration(){
    // cout << " traditional Iteration ..."  << "\n";
    //传播和计算叠加进行
    activeVertexList = new vector<uintV>();
    activeVertexList->reserve((n >> 10)+1);
    for(long i = 0; i < affectedVertexSet->size(); i++){
      uintV v = (*affectedVertexSet)[i];
      if(versionId[v] != batch_versionId){
        versionId[v] = batch_versionId;
        initializeAggregationValue<AggregationValueType, GlobalInfoType>(v, av[v], batch_versionId, global_info);             
        VersionGraphIterator<vertex> scanInEdgeItr(my_graph, v, batch_versionId, false);
        versionEdge * tmpRef = nullptr;
        while(scanInEdgeItr.hasNext() && scanInEdgeItr.next(&tmpRef)){
          uintV t = tmpRef->vertex;//遍历点的入边
          if(preIterationVSS->active[t]){
            edgeFunction<AggregationValueType, VertexValueType, EdgeData, GlobalInfoType>
                (t, v, batch_versionId, iter, tmpRef->weight, preIterationVSS->vv[t], av[v], global_info);
          }
        }
        flags[v] = AGGREGATION_VALID;
        VertexValueType vertex_value_tmp = preIterationVSS->vv[v];
        bool active_this_iter = computeFunction<AggregationValueType, VertexValueType, GlobalInfoType>
                                    (v, batch_versionId, iter, av[v], preIterationVSS->vv[v], vertex_value_tmp, global_info);
        vv[v] = vertex_value_tmp;
        active[v] = active_this_iter;
        if(active_this_iter){
          active_num++;
          activeVertexList->push_back(v);
        }
      }

      VersionGraphIterator<vertex> scanOutEdgeItr(my_graph, v, batch_versionId, true);
      versionEdge * vTmpRef = nullptr;     
      while(scanOutEdgeItr.hasNext() && scanOutEdgeItr.next(&vTmpRef)){
        uintV u = vTmpRef->vertex;//遍历激活点的邻居
        if(versionId[u] != batch_versionId){
          versionId[u] = batch_versionId;       
          initializeAggregationValue<AggregationValueType, GlobalInfoType>(u, av[u], batch_versionId, global_info);
          VersionGraphIterator<vertex> scanInEdgeItr(my_graph, u, batch_versionId, false);
          versionEdge * uTmpRef = nullptr; 
          while(scanInEdgeItr.hasNext() && scanInEdgeItr.next(&uTmpRef)){
            uintV t = uTmpRef->vertex;//遍历点的入边
            if(preIterationVSS->active[t]){
              edgeFunction<AggregationValueType, VertexValueType, EdgeData, GlobalInfoType>
                  (t, u, batch_versionId, iter, uTmpRef->weight, preIterationVSS->vv[t], av[u], global_info);
            }
          }
          flags[u] = AGGREGATION_VALID;
          VertexValueType vertex_value_tmp = preIterationVSS->vv[u];
          bool active_this_iter = computeFunction<AggregationValueType, VertexValueType, GlobalInfoType>
                                      (u, batch_versionId, iter, av[u], preIterationVSS->vv[u], vertex_value_tmp, global_info);
          vv[u] = vertex_value_tmp;
          active[u] = active_this_iter;
          if(active_this_iter){
            active_num++;
            activeVertexList->push_back(u);
          }
        }
      }
    }
    if(affectedVertexSet != nullptr)  delete affectedVertexSet;
    affectedVertexSet = activeVertexList;
    activeVertexList = nullptr;
    if(affectedVertexSet->empty()){
      state = CONVERGED_FLAG;
    }
    // cout << " affectedVertexSet.size(): "<< affectedVertexSet->size()<< " activeVertexList.size(): "<< activeVertexList->size()<< "\n";
	}

  void freeVSSData(){
    if(vv != nullptr){
      deleteA(av);
      deleteA(vv);
      deleteA(active);
      deleteA(versionId);
      deleteA(flags);
      av = nullptr;
      vv = nullptr;
      active = nullptr;
      versionId = nullptr;
      flags = nullptr;      
    }
    n = 0;
    c = 0;
    baseVSS = nullptr;
    preIterationVSS = nullptr;
    if(affectedVertexSet != nullptr){ delete affectedVertexSet;  affectedVertexSet = nullptr;}  
    if(activeVertexList != nullptr){  delete activeVertexList; activeVertexList = nullptr;}
    if(affectedVertexSetHistoryValue != nullptr){  delete affectedVertexSetHistoryValue; affectedVertexSetHistoryValue = nullptr;}
    if(affectedVertexSetActive != nullptr){  delete affectedVertexSetActive; affectedVertexSetActive = nullptr;}
  }

  void printVSS(bool detail){
    cout << "***********************************************" << "\n";
    cout << " batch_versionId: " <<batch_versionId <<" iter: "<<iter <<" state: "<<((state == CONVERGED_FLAG )? "Converged" : "Iterating") << "\n";
    if(baseVSS == nullptr)
      cout << " baseVss: null" << "\n";
    else
      cout << " baseVss: batch " <<baseVSS->batch_versionId <<" iter "<<baseVSS->iter <<" state: "
        <<((baseVSS->state) == CONVERGED_FLAG ? "Converged" : "Iterating" )<< "\n";

   if(preIterationVSS == nullptr)
      cout << " preIterationVSS: null" << "\n";
    else
      cout << " preIterationVSS: batch " <<preIterationVSS->batch_versionId <<" iter "<<preIterationVSS->iter 
        <<" state: "<<((preIterationVSS->state) == CONVERGED_FLAG ? "Converged" : "Iterating" )<< "\n";

    cout << " n: " <<n <<" c: "<<c <<" valid_n: "<<valid_n << "\n";

    if(affectedVertexSet == nullptr)
      cout << " affectedVertexSet: null" << "\n";
    else{
      cout << " affectedVertexSet: size() " <<affectedVertexSet->size() << "\n";
      if(detail)
      for(long i = 0; i < affectedVertexSet->size();i++){
        cout << (*affectedVertexSet)[i] << "; " ;
      }
      cout << "\n";
    }

    if(vv == nullptr)
      cout << " VSS content: null" << "\n";
    else{
      if(detail)
      for(long i = 0; i < n; i++){
        cout << "vv[" <<i <<"]\t"<<vv[i].val<< "\t" <<vv[i].outDegree  << "\tav[" <<i <<"]\t"<<av[i] 
          << "\t" <<active[i] << "\t" <<versionId[i] <<"\t"<<flags[i] << "\n";
      }
    }
    cout << "***********************************************" << "\n";
  }
};

template<class vertex, class AggregationValueType, class VertexValueType, class GlobalInfoType>
void exec(void * const VSS, const int x, const int y){
  ((VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType> * const)VSS)->run();
  Pos returnResult(x,y);
  execResultQueue.enqueue(returnResult);
}

// ======================================================================
// MORPHOSES ENGINE
// ======================================================================
template <class vertex, class AggregationValueType, class VertexValueType,
          class GlobalInfoType>
class MorphosesEngine {

public:
  //保存传入的参数
  graph<vertex> & my_graph;
  int max_iterations;
  int max_rows;
  int max_batchSize;
  GlobalInfoType &global_info;
  commandLine config;

  // Stream Ingestor
  Ingestor<vertex> ingestor;
  uintEE current_batch;

  //事务调度数据结构
  task * taskList;
  int taskPtr = 0;
  task * historyTask;
  int historyIteration = 0;
  VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType> ** historyVSS;
  VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType> ** VSSSchedularTable;
  int * flagSchedularTable;
  //线程数组
  int nWorkers;
  ThreadPool * pool;
  // ======================================================================
  // CONSTRUCTOR / INIT
  // ======================================================================
  MorphosesEngine(graph<vertex> &_my_graph, int _max_iter,
                  GlobalInfoType &_global_info, commandLine _config)
      : my_graph(_my_graph), max_iterations(_max_iter),
        global_info(_global_info), config(_config), 
        ingestor(_my_graph, _config), current_batch(0) {
    max_batchSize = 4;
    nWorkers = TEST_THREAD_NUM;
    max_rows = max_iterations + 3;
    pool = new ThreadPool(nWorkers);
  }

  void init() {}

  ~MorphosesEngine() {
    global_info.cleanup();
  }

  void initialCompute() {
    global_info.init();

    taskList = newA(task , max_batchSize);
    for(int i = 0; i < max_batchSize; i++){
      taskList[i].flag = INVALID_FLAG;
      taskList[i].version = -1;
      taskList[i].iteration = -1;
    }
    historyIteration = -1;
    historyVSS = (VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType> **)newA(void * ,max_rows);
    VSSSchedularTable = (VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType> **)newA(void * , max_batchSize * max_rows);
    flagSchedularTable = newA(int , max_batchSize * max_rows);
    for(int i = 0; i < max_rows; i++){
      historyVSS[i] = nullptr;
    }
    for(int i = 0; i < max_rows; i++){
      for(int j = 0; j < max_batchSize; j++){
        flagSchedularTable[i*max_batchSize + j] = INVALID_FLAG;
        VSSSchedularTable[i*max_batchSize + j] = nullptr;
      }
    }
    //启动线程池
    pool->init();
  }

  void run() {
    initialCompute();
    // ======================================================================
    // Incremental Compute - Get the next update batch from ingestor
    // ======================================================================
    ingestor.validateAndOpenFifo();

    while (ingestor.processNextBatch()) {
      current_batch++;
      edgeArray &edge_additions = ingestor.getEdgeAdditions();
      edgeArray &edge_deletions = ingestor.getEdgeDeletions();
      //将增加的边补写度数到delta
      for(long i = 0; i < edge_additions.size; i++) {
        uintV source = edge_additions.E[i].source;
        uintV destination = edge_additions.E[i].destination;
        long current_version = current_batch;

        my_graph.graph_delta -> insertVertexDelta(current_batch, source,
        my_graph.V[source].getOutDegree(), my_graph.V[source].getInDegree());

        my_graph.graph_delta -> insertVertexDelta(current_batch, destination,
        my_graph.V[source].getOutDegree(), my_graph.V[source].getInDegree());
      }
      //将删除的边补写到delta
      for(long i = 0; i < edge_deletions.size; i++) {
        uintV source = edge_deletions.E[i].source;
        uintV destination = edge_deletions.E[i].destination;
        long deletion_version = edge_deletions.edgeDataArray[i].edgeVersion;
        long current_version = current_batch;

        my_graph.graph_delta -> insertEdgeDelta(deletion_version, current_batch, 
        source, destination, true, edge_deletions.edgeDataArray[i],
        my_graph.V[source].getOutDegree(), my_graph.V[source].getInDegree());

        my_graph.graph_delta -> insertEdgeDelta(deletion_version, current_batch, 
        destination, source, false, edge_deletions.edgeDataArray[i],
        my_graph.V[destination].getOutDegree(), my_graph.V[destination].getInDegree());
      }
      //写入任务列表
      if ((edge_additions.size == 0) && (edge_deletions.size == 0)) {  
        continue;
      }
      if(taskPtr < max_batchSize){
        taskList[taskPtr].flag = TODO_FLAG;
        taskList[taskPtr].version = current_batch;
        taskList[taskPtr].iteration = 0;

        taskList[taskPtr].addition.copyEdgeArray(edge_additions);
        taskList[taskPtr].deletion.copyEdgeArray(edge_deletions);
        taskPtr++;
      }
      if(taskPtr >= max_batchSize){
        deltaCompute();
        cleanAndMerge();
      }

      // my_graph.graph_delta->print(true);
      // cout << "****************Graph Delta***************" << "\n";
      // for(long i = 0; i < edge_additions.size; i++) {
      //   uintV source = edge_additions.E[i].source;
      //   uintV destination = edge_additions.E[i].destination;
      //   long addition_version = edge_additions.edgeDataArray[i].edgeVersion;
      //   long current_version = current_batch;
      //   cout << "ADD " << source <<" TO " << destination <<" eV " << addition_version << " tV " << current_version << "\n";
      // }
      // for(long i = 0; i < edge_deletions.size; i++) {
      //   uintV source = edge_deletions.E[i].source;
      //   uintV destination = edge_deletions.E[i].destination;
      //   long deletion_version = edge_deletions.edgeDataArray[i].edgeVersion;
      //   long current_version = current_batch;
      //   cout << "DEL " << source <<" TO " << destination <<" eV " << deletion_version << " tV " << current_version << "\n";
      // }
      // printGraph(true);

    }

    cout << "Stop threads ..." << "\n";
    pool->shutdown();
  }

  void createVSSInSchedularTable(int row_pos, int col_pos){
    // cout << "createVSSInSchedularTable [" << row_pos <<"]["<<col_pos <<"] ..."<< "\n";
    VSSSchedularTable[row_pos*max_batchSize + col_pos] = new VertexStateSet<vertex, AggregationValueType, VertexValueType, GlobalInfoType>(
        my_graph,
        taskList[col_pos].version,
        taskList[col_pos].iteration,
        my_graph.n, 
        global_info, 
        col_pos > 0 ? VSSSchedularTable[row_pos*max_batchSize + col_pos - 1] : historyVSS[row_pos], 
        row_pos > 0 ? VSSSchedularTable[(row_pos - 1)*max_batchSize + col_pos] : nullptr, 
        taskList[col_pos].addition,
        taskList[col_pos].deletion
    );
    flagSchedularTable[row_pos*max_batchSize + col_pos] = READY_FLAG;
  }

  void freeVSSInSchedularTable(int row_pos, int col_pos){
    // cout << "freeVSSInSchedularTable [" << row_pos <<"]["<<col_pos <<"] ..."<< "\n";
    if(VSSSchedularTable[row_pos*max_batchSize + col_pos] != nullptr){
      VSSSchedularTable[row_pos*max_batchSize + col_pos] -> freeVSSData();
      delete VSSSchedularTable[row_pos*max_batchSize + col_pos];
      VSSSchedularTable[row_pos*max_batchSize + col_pos] = nullptr;      
    }
    flagSchedularTable[row_pos*max_batchSize + col_pos] = INVALID_FLAG;
  }  

  bool checkVSSCanSchedule(int row_pos, int col_pos){
    // cout << "checkVSSCanSchedule [" << row_pos <<"]["<<col_pos <<"] ..."<< "\n";
    if(row_pos < 0 || row_pos >= max_rows)  return false;
    if(col_pos < 0 || col_pos >= max_batchSize) return false;
    //事务存在,迭代数未超出且非收敛
    if(taskList[col_pos].flag == INVALID_FLAG || taskList[col_pos].flag == CONVERGED_FLAG || taskList[col_pos].iteration > max_iterations)  
      return false;
    //当前迭代处于无效状态
    if(flagSchedularTable[row_pos*max_batchSize + col_pos] != INVALID_FLAG) return false;
    //左侧VSS完成
    bool baseDone = col_pos == 0 || taskList[col_pos-1].flag == CONVERGED_FLAG || flagSchedularTable[row_pos*max_batchSize + col_pos - 1] == CONTINUE_FLAG;
    //上方VSS完成
    bool preItrDone = row_pos == 0 || flagSchedularTable[(row_pos-1)*max_batchSize + col_pos] == DONE_FLAG;
    return baseDone && preItrDone;
  }

  void deltaCompute(){
    long long parallelStastics[max_iterations] = {0};
    timer t;
    cout << "deltaCompute ... " << "\n";
    t.start();
    int runningVSSNum = 0;
    //调度左上角事务
    if(checkVSSCanSchedule(0,0)){
      // cout << "checkVSSCanSchedule [" << 0 <<"]["<<0 <<"] true"<< "\n";
      createVSSInSchedularTable(0,0);
      runningVSSNum++;
      pool->submit(
        exec<vertex, AggregationValueType, VertexValueType, GlobalInfoType>, 
        VSSSchedularTable[0*max_batchSize + 0], 0, 0
      );
      taskList[0].flag = RUNNING_FLAG;    
      taskList[0].iteration++;      
    }
    else{
      return;
    }
    Pos returnResult(0,0);
    while(true){
      // cout << "runningVSSNum: " <<runningVSSNum<< "\n";
      //检查完成的任务
      if(execResultQueue.dequeue(returnResult)){
        int row_pos = returnResult.row_pos;
        int col_pos = returnResult.col_pos;
        runningVSSNum--;
        // cout << "doneVSS [" << row_pos <<"]["<<col_pos <<"] ..."<< "\n";       
        flagSchedularTable[row_pos*max_batchSize + col_pos] = DONE_FLAG;
        if(row_pos > 0) flagSchedularTable[(row_pos-1)*max_batchSize + col_pos] = CONTINUE_FLAG;
        // VSSSchedularTable[row_pos*max_batchSize + col_pos]->printVSS(true);
        //检查迭代是否收敛
        if(row_pos >= max_iterations || VSSSchedularTable[row_pos*max_batchSize + col_pos]->state == CONVERGED_FLAG)
          taskList[col_pos].flag = CONVERGED_FLAG;
        //调度下一轮迭代
        if(checkVSSCanSchedule(row_pos + 1, col_pos)){
          createVSSInSchedularTable(row_pos + 1, col_pos);
          runningVSSNum++;
          pool->submit(
            exec<vertex, AggregationValueType, VertexValueType, GlobalInfoType>, 
            VSSSchedularTable[(row_pos + 1)*max_batchSize + col_pos], 
            row_pos + 1,
            col_pos
          );
          taskList[col_pos].flag = RUNNING_FLAG;    
          taskList[col_pos].iteration++;
        }
        //调度右上方迭代
        if(checkVSSCanSchedule(row_pos - 1, col_pos + 1)){
          createVSSInSchedularTable(row_pos - 1, col_pos + 1);
          runningVSSNum++;
          pool->submit(
            exec<vertex, AggregationValueType, VertexValueType, GlobalInfoType>, 
            VSSSchedularTable[(row_pos - 1)*max_batchSize + col_pos + 1], 
            row_pos - 1,
            col_pos + 1
          );
          taskList[col_pos + 1].flag = RUNNING_FLAG;    
          taskList[col_pos + 1].iteration++;
        }
      }
      else{
        int x; while(x < 1) x++;
      }
      if(runningVSSNum == 0)  break;

      parallelStastics[runningVSSNum]++;
    }
    cout << "deltaCompute use time:  " <<t.stop() << "\n";
    for(int i = 0; i < max_iterations; i++){
      cout << "para: " << i <<" num: " << parallelStastics[i] << endl;
    }
  }

  void cleanAndMerge(){
    cout << "cleanAndMerge ... " << "\n";
    cout << "historyIteration: "<<historyIteration << "\n";
    //清空历史
    for(int i = 0; i < historyIteration; i++){
      if(historyVSS[i] != nullptr){
        historyVSS[i] -> freeVSSData();
        delete historyVSS[i];
        historyVSS[i] = nullptr;      
      }
    }      
    historyIteration = 0;
    //将最新一轮迭代拷贝到历史
    historyIteration = taskList[max_batchSize - 1].iteration;
    for(int i = 0; i < historyIteration; i++){
      historyVSS[i] = VSSSchedularTable[i*max_batchSize + max_batchSize - 1];
      VSSSchedularTable[i*max_batchSize + max_batchSize - 1] = nullptr;
    } 
    cout << "historyIteration: "<<historyIteration << "\n";
    //释放所有块的内容
    for(int i = 0; i < max_rows; i++){
      for(int j = 0; j < max_batchSize; j++){
        freeVSSInSchedularTable(i, j);
      }
    }
    for(int i = 0; i < taskPtr; i++){
      // cout << "taskList["<<i <<"]: " << " version " << taskList[i].version << " iteration " << taskList[i].iteration  << "\n";
      taskList[i].flag = INVALID_FLAG;
      taskList[i].iteration = -1;
      taskList[i].version = -1;
    }
    taskPtr = 0;    
    my_graph.graph_delta->mergeVersion();

  } 

  void printGraph(long version, bool detail){
    cout << "****************Graph Start***************" << "\n";
    cout << "M Edges: " << my_graph.m << endl;
    cout << "N Vertex: " << my_graph.n << endl;
    cout << "version: " << version << endl;
    uintV numEdgesFromDegree = 0;
    for (uintV i = 0; i < my_graph.n; i++) {
      numEdgesFromDegree += my_graph.getOutDegree(i, version);
    }

    if (my_graph.m != numEdgesFromDegree) {
      cout << "~~~~~~~~~Edges ARE NOT EQUAL!!!!~~~~~~~~~" << endl;
      cout << "m: " << my_graph.m << " NumEdges: " << numEdgesFromDegree << endl;
    }
    if(detail){
    for (uintV i = 0; i < my_graph.n; ++i) {
        cout << i << ": outDegree: " <<my_graph.getOutDegree(i, version) << " inDegree: " <<my_graph.getInDegree(i, version) << "\n";
    }
    cout << "****************Graph Out***************" << "\n";
    for (uintV i = 0; i < my_graph.n; ++i) {
      VersionGraphIterator<vertex> scanOutEdgeItr(my_graph, i, version, true);
      versionEdge * tmpRef = nullptr;
      while(scanOutEdgeItr.hasNext() && scanOutEdgeItr.next(&tmpRef)){
        cout << i << "->" <<tmpRef->vertex << "\n";
      }
    }
    cout << "*****************Graph In*******************" << "\n";
    for (uintV i = 0; i < my_graph.n; ++i) {
      VersionGraphIterator<vertex> scanInEdgeItr(my_graph, i, version, false);
      versionEdge * tmpRef = nullptr;
      while(scanInEdgeItr.hasNext() && scanInEdgeItr.next(&tmpRef)){
        cout << i << "<-" <<tmpRef->vertex << "\n";
      }
    }
    }
    cout << "*******************Graph End****************" << "\n";
  }
};


#endif
