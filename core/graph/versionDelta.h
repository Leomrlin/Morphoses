
#ifndef VERSION_DELTA_H
#define VERSION_DELTA_H
#include "../graph/edgeDataType.h"
#include "../common/utils.h"
#include "../graph/versionGraph.h"
#include <vector>
struct versionEdge{
    long startVersion;
    long endVersion;
    uintV vertex;
    EdgeData weight;
};

struct versionVertex{
    long version;
    uintEE outDegree;
    uintEE inDegree;
};

class GraphVersionDelta{
public:
  vector<vector<versionVertex> > V;
  vector<vector<versionEdge> > outEdges;
  vector<vector<versionEdge> > inEdges;
  vector<long> deltaVersionTable;
  vector<uintV> deltaVertexSet;
  uintEE baseVersion;
  uintEE currentVersion;
  uintV n;

  GraphVersionDelta(uintEE _baseVersion, uintV _n) :
  baseVersion(_baseVersion), currentVersion(_baseVersion), n(_n){
    V.resize(n);
    outEdges.resize(n);
    inEdges.resize(n);
    deltaVersionTable.resize(n);
    for(uint i = 0; i < n ; i++)
        deltaVersionTable[i] = -1;
    deltaVertexSet.resize(0);
    
  }

  void insertEdgeDelta(long startVersion, long endVersion, uintV v,
              uintV u, bool outEdge, EdgeData weight, 
              uintEE outDegree, uintEE inDegree){
    //执行扩容操作
    if(v > n){
        V.resize(v);
        outEdges.resize(v);
        inEdges.resize(v);
        deltaVersionTable.resize(v);
        for(uint i = n; i < v ; i++)
            deltaVersionTable[i] = -1;
        n = v;
    }
    //插入点记录,并记录受影响的点
    insertVertexDelta(endVersion, v, outDegree, inDegree);
    //插入边记录
    versionEdge tmp;
    tmp.startVersion = startVersion;
    tmp.endVersion = endVersion;
    tmp.vertex = u;
    tmp.weight.setEdgeDataFromPtr(&weight);
    if(outEdge)
        outEdges[v].push_back(tmp);
    else
        inEdges[v].push_back(tmp);
    if(endVersion > currentVersion) currentVersion = endVersion;
  }

  void insertVertexDelta(long version, uintV v, uintEE outDegree, uintEE inDegree){
    // cout <<"insertVertexDelta v: " << v<<" version: " <<version<<" outD: " <<outDegree <<" inD: "<<inDegree << endl;
    if(version > currentVersion) currentVersion = version;
    //如果记录比点原始记录更新，则更新点度数记录
    if(deltaVersionTable[v] < version){
        versionVertex tmpv;
        tmpv.version = version;
        tmpv.outDegree = outDegree;
        tmpv.inDegree = inDegree;
        V[v].push_back(tmpv);
        deltaVersionTable[v] = version;
        deltaVertexSet.push_back(v);
    }
  }

  void mergeVersion(){
    int vertexNum = deltaVertexSet.size();
    for(int i = 0; i < vertexNum; i++){
        uintV v = deltaVertexSet[i];
        V[v].clear();
        outEdges[v].clear();
        inEdges[v].clear();
        deltaVersionTable[v] = -1;
    }
    deltaVertexSet.clear();
    baseVersion = currentVersion;
  }

  void print(bool detail){
    cout << "****************Graph Delta Start***************" << "\n";
    cout << "baseVersion: " << baseVersion << endl;
    cout << "currentVersion: " << currentVersion << endl;
    cout << "n: " << n << endl;

    if(detail){
        int vertexNum = deltaVertexSet.size();
        for(int i = 0; i < vertexNum; i++){
            uintV v = deltaVertexSet[i];
            cout << "v: " << v << " deltaVersion[v]: " << deltaVersionTable[v] << endl;
            for(versionVertex x : V[v])
                cout << x.version <<"\t"<< x.outDegree <<"\t"<< x.inDegree <<endl;
            for(versionEdge x : outEdges[v])
               cout << x.startVersion <<"\t"<< x.endVersion <<"\t"<< x.vertex<<"\t"<< x.weight.edgeVersion <<endl;
            for(versionEdge x : inEdges[v])
                cout << x.startVersion <<"\t"<< x.endVersion <<"\t"<< x.vertex<<"\t"<< x.weight.edgeVersion <<endl;
        }
    }

    cout << "*******************Graph Delta End****************" << "\n";
  }
  
  ~GraphVersionDelta(){
      for(uintV i = 0; i < n; i++){
          V[i].clear();
          outEdges[i].clear();
          inEdges[i].clear();
          deltaVersionTable[i] = currentVersion;
      }
      deltaVertexSet.clear();
  }

};
#endif
