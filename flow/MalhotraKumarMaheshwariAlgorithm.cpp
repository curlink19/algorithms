/// to be: Tvertex -> size_t
#include <bits/stdc++.h>

template<typename Tvertex, typename Tflow>
class Edge {
private:
  Tvertex start;
  Tvertex finish;
  Tflow capacity;
  Tflow flow;

public:
  Edge(Tvertex u, Tvertex v, Tflow c, Tflow f):
    start(u),
    finish(v),
    capacity(c),
    flow(f) {}

  Tvertex getStart() const {
    return start;
  }

  Tvertex getFinish() const {
    return finish;
  }

  Tflow getCapacity() const {
    return capacity;
  }

  Tflow getFlow() const {
    return flow;
  }

  Tflow getResidualCapacity() const {
    return capacity - flow;
  }

  void pushFlow(Tflow value) {
    flow += value;
  }
};

template<typename Tvertex, typename Tflow>
class Network;

template<typename Tvertex, typename Tflow>
class outgoingEdgesIterator {
  using Tnetwork = Network<Tvertex, Tflow>;
  using Tedge = Edge<Tvertex, Tflow>;

private:
  Tnetwork& network;
  int position;
  
public:
  outgoingEdgesIterator(Tnetwork& net, Tvertex v):
    network(net),
    position(network.lastOutgoingEdge[v]) {}

  outgoingEdgesIterator(Tnetwork& net):
    network(net),
    position(-1) {}

  outgoingEdgesIterator(const outgoingEdgesIterator& other):
    network(other.network),
    position(other.position) {}

  outgoingEdgesIterator& operator=(const outgoingEdgesIterator& other) {
    position = other.position;
    return *this;
  }

  bool valid() const {
    return position != -1;
  }

  void next() {
    position = network.previousEdge[position];
  }

  Tedge getEdge() const {
    return network.edges[position];
  }

  Tvertex getStart() const {
    return getEdge().getStart();
  }

  Tvertex getFinish() const {
    return getEdge().getFinish();
  }

  Tflow getCapacity() const {
    return getEdge().getCapacity();
  }

  Tflow getFlow() const {
    return getEdge().getFlow();
  }

  Tflow getResidualCapacity() const {
    return getEdge().getResidualCapacity();
  }

  void pushFlow(Tflow value) {
    network.edges[position].pushFlow(value);
    network.edges[position ^ 1].pushFlow(-1 * value);
  }

  friend Network<Tvertex, Tflow>;
};


template<typename Tvertex, typename Tflow>
class Network {
public:
  using Tedge = Edge<Tvertex, Tflow>;

private:
  int vertexCount;
  std::vector<Tedge> edges;
  std::vector<int> previousEdge;
  std::vector<int> lastOutgoingEdge;

public:
  Network(int n):
    vertexCount(n),
    edges(),
    previousEdge(),
    lastOutgoingEdge(std::vector<int>(n, -1)) {}

public:
  using Titerator = outgoingEdgesIterator<Tvertex, Tflow>;

  Titerator getEdgesIterator(Tvertex v) {
    return Titerator(*this, v);
  }

  int getVertexCount() const {
    return vertexCount;
  }

  int getEdgesNumber() const {
    return static_cast<int>(edges.size());
  }

private:
  void insertEdgeLocal(Tvertex newStart, Tvertex newFinish, Tflow newCapacity, Tflow newFlow) {
    Tedge newEdge(newStart, newFinish, newCapacity, newFlow);
    int newIndex = getEdgesNumber();
    edges.push_back(newEdge);
    previousEdge.push_back(lastOutgoingEdge[newStart]);
    lastOutgoingEdge[newStart] = newIndex;
  }

public:
  void insertEdge(Tvertex newStart, Tvertex newFinish, Tflow newCapacity) {
    insertEdgeLocal(newStart, newFinish, newCapacity, 0);
    insertEdgeLocal(newFinish, newStart, 0, 0);
  }

  friend outgoingEdgesIterator<Tvertex, Tflow>;
};

template<typename Tvertex, typename Tflow>
class blockPathFlowAlgorithm {
private:
  using Tnetwork = Network<Tvertex, Tflow>;
  using Tedge = Edge<Tvertex, Tflow>;
  using Titerator = outgoingEdgesIterator<Tvertex, Tflow>;

protected:
  Tnetwork& network;
  Tvertex source;
  Tvertex sink;
  std::vector<int> dist;

public:
  blockPathFlowAlgorithm(Tnetwork& net, Tvertex newSource, Tvertex newSink):
    network(net),
    source(newSource),
    sink(newSink),
    dist(std::vector<int>(net.getVertexCount(), -1)) {}
    
    Tvertex getSource() const {
      return source;
    }

    Tvertex getSink() const {
      return sink;
    }

    int getVertexCount() const {
      return network.getVertexCount();
    }

private:
  /// is any way to sink 
  bool calculateDist() {
    dist = std::vector<int>(getVertexCount(), -1);
    dist[source] = 0;
    std::queue<Tvertex> que;
    que.push(source);
    while (!que.empty()) {
      Tvertex startVertex = que.front();
      que.pop();
      for (auto edge = network.getEdgesIterator(startVertex); edge.valid(); edge.next()) {
        if (edge.getResidualCapacity() == 0) {
          continue;
        }
        int currentDist = dist[edge.getFinish()];
        if ((currentDist == -1) || (dist[startVertex] + 1 < currentDist)) {
          dist[edge.getFinish()] = dist[startVertex] + 1;
          que.push(edge.getFinish());
        }
      }
    }
    if (dist[sink] != -1) {
      for (Tvertex startVertex = 0; startVertex < getVertexCount(); ++startVertex) {
        if (startVertex == sink) {
          continue;
        }
        if (dist[startVertex] >= dist[sink]) {
          dist[startVertex] = -1;
        }
      } /// not full
    }
    return dist[sink] != -1;
  }

protected:
  virtual void findBlockFlow() = 0;

public:
  Tflow getMaximalFlow() {
    while (calculateDist()) {
      findBlockFlow();
    }
    Tflow count = 0;
    for (auto edge = network.getEdgesIterator(sink); edge.valid(); edge.next()) {
      count -= edge.getFlow();
    }
    return count;
  }
};

template<typename Tvertex, typename Tflow>
class MalhotraKumarMaheshwariAlgorithm: public blockPathFlowAlgorithm<Tvertex, Tflow> {
  using blockPathFlowAlgorithm<Tvertex, Tflow>::blockPathFlowAlgorithm;

private:
  using Tnetwork = Network<Tvertex, Tflow>;
  using Tedge = Edge<Tvertex, Tflow>;
  using Titerator = outgoingEdgesIterator<Tvertex, Tflow>;

private:
  std::vector<Titerator> currentEdge;
  std::vector<bool> isDeleted;
  std::vector<Tflow> inPotential;
  std::vector<Tflow> outPotential;

  std::vector<std::vector<Titerator>> backwardEdges;
  std::vector<int> currentBackwardEdge;

  bool isInDistNetwork(const Titerator& edge) const {
    return ((this->dist[edge.getStart()] + 1) == this->dist[edge.getFinish()]) &&
            (this->dist[edge.getStart()] != -1) && (this->dist[edge.getFinish()] != -1);
  }

  Tflow getPotential(Tvertex v) const {
    if (v == this->getSource()) {
      return outPotential[v];
    }
    if (v == this->getSink()) {
      return inPotential[v];
    }
    return std::min(inPotential[v], outPotential[v]);
  }

  void initialize() {
    /// to be currentEdge!!!
    currentEdge.clear();
    isDeleted = std::vector<bool>(this->getVertexCount(), false);
    inPotential = std::vector<Tflow>(this->getVertexCount(), 0);
    outPotential = std::vector<Tflow>(this->getVertexCount(), 0);
    for (Tvertex startVertex = 0; startVertex < this->getVertexCount(); ++startVertex) {
      currentEdge.push_back(this->network.getEdgesIterator(startVertex));
      for (auto edge = this->network.getEdgesIterator(startVertex); edge.valid(); edge.next()) { 
        if (!isInDistNetwork(edge)) {
          continue;
        }
        inPotential[edge.getFinish()] += edge.getResidualCapacity();
        outPotential[edge.getStart()] += edge.getResidualCapacity();
      }
    }
    for (Tvertex startVertex = 0; startVertex < this->getVertexCount(); ++startVertex) {
      if (this->dist[startVertex] == -1) {
        isDeleted[startVertex] = true;
      }
    }
    /// Now initialize backward edges (in dist network!!!) 
    backwardEdges = std::vector<std::vector<Titerator>>(this->getVertexCount());
    currentBackwardEdge = std::vector<int>(this->getVertexCount(), 0);  /// bad if >= size
    for (Tvertex startVertex = 0; startVertex < this->getVertexCount(); ++startVertex) {
      for (auto edge = this->network.getEdgesIterator(startVertex); edge.valid(); edge.next()) {
        if (!isInDistNetwork(edge)) {
          continue;
        }
        backwardEdges[edge.getFinish()].push_back(edge);
      }
    }
  }

  /// returns vertexCount if ALL vertex are DELETED
  Tvertex getMinimalPotentialVertex() const {
    Tvertex result = this->getVertexCount();
    for (Tvertex startVertex = 0; startVertex < this->getVertexCount(); ++startVertex) {
      if (isDeleted[startVertex]) {
        continue;
      }
      if ((result == this->getVertexCount()) || (getPotential(startVertex) < getPotential(result))) {
        result = startVertex;
      }
    }
    return result;
  }
 
  Tflow propogate(Tvertex startVertex, Tflow flowValue) {
    if ((flowValue == 0) || (startVertex == this->getSink())) {
      return flowValue;
    }
    if (isDeleted[startVertex]) {
      return 0;
    }
    Tflow result = 0;
    while ((currentEdge[startVertex].valid()) && (flowValue > 0)) {
      if (!isInDistNetwork(currentEdge[startVertex])) {
        currentEdge[startVertex].next();
        continue;
      }
      Tflow limit = std::min(flowValue, currentEdge[startVertex].getResidualCapacity());
      Tflow propogateValue = propogate(currentEdge[startVertex].getFinish(), limit);
      currentEdge[startVertex].pushFlow(propogateValue);
      result += propogateValue;
      flowValue -= propogateValue;
      inPotential[currentEdge[startVertex].getFinish()] -= propogateValue;
      outPotential[currentEdge[startVertex].getStart()] -= propogateValue;
      if ((flowValue > 0) || (currentEdge[startVertex].getResidualCapacity() == 0)) {
        currentEdge[startVertex].next();
        continue;
      }
    }
    return result;
  }

  int getBackwardEdgesCount(Tvertex startVertex) const {
    return static_cast<int>(backwardEdges[startVertex].size());
  }

  Tflow backwardPropogate(Tvertex startVertex, Tflow flowValue) {
    if ((flowValue == 0) || (startVertex == this->getSource())) {
      return flowValue;
    }
    if (isDeleted[startVertex]) {
      return 0;
    }
    Tflow result = 0;
    while ((currentBackwardEdge[startVertex] < getBackwardEdgesCount(startVertex)) && (flowValue > 0)) {
      Tflow limit = std::min(flowValue, backwardEdges[startVertex][currentBackwardEdge[startVertex]].getResidualCapacity());
      Tflow propogateValue = backwardPropogate(backwardEdges[startVertex][currentBackwardEdge[startVertex]].getStart(), limit);
      backwardEdges[startVertex][currentBackwardEdge[startVertex]].pushFlow(propogateValue);
      result += propogateValue;
      flowValue -= propogateValue;
      inPotential[backwardEdges[startVertex][currentBackwardEdge[startVertex]].getFinish()] -= propogateValue;
      outPotential[backwardEdges[startVertex][currentBackwardEdge[startVertex]].getStart()] -= propogateValue;
      if ((flowValue > 0) || (backwardEdges[startVertex][currentBackwardEdge[startVertex]].getResidualCapacity() == 0)) {
        ++currentBackwardEdge[startVertex];
        continue;
      }
    }
    return result;
  }

  void safeDelete(Tvertex startVertex) {
    isDeleted[startVertex] = true;
    for (auto edge = this->network.getEdgesIterator(startVertex); edge.valid(); edge.next()) {
      if (!isInDistNetwork(edge)) {
        continue;
      }
      inPotential[edge.getFinish()] -= edge.getResidualCapacity();
    }
    for (int i = 0; i < getBackwardEdgesCount(startVertex); ++i) {
      auto edge = backwardEdges[startVertex][i];
      outPotential[edge.getStart()] -= edge.getResidualCapacity();
    }
  }

/// find all to be deleted!!!
/// and calculate potential
protected:
  void findBlockFlow() override {
    initialize();
    Tvertex startVertex = getMinimalPotentialVertex(); 
    while (startVertex != this->getVertexCount()) {
      Tflow flowValue = getPotential(startVertex);
      propogate(startVertex, flowValue);
      backwardPropogate(startVertex, flowValue);
      safeDelete(startVertex);
      startVertex = getMinimalPotentialVertex();
    }
    return;
  }
};

int main() {
  return 0;
}
