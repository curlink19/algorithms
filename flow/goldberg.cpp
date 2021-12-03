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
class preflowMaximalFlowAlgorithm {
private:
  using Tnetwork = Network<Tvertex, Tflow>;
  using Tedge = Edge<Tvertex, Tflow>;
  using Titerator = outgoingEdgesIterator<Tvertex, Tflow>;

private:
  Tnetwork& network;
  Tvertex source;
  Tvertex sink;
  std::vector<Tflow> excess;
  std::vector<int> height;
  std::vector<Titerator> dischargePointer;

public:
  preflowMaximalFlowAlgorithm(Tnetwork& net, Tvertex newSource, Tvertex newSink):
    network(net),
    source(newSource),
    sink(newSink),
    excess(std::vector<Tflow>(net.getVertexCount(), 0)), 
    height(std::vector<int>(net.getVertexCount(), 0)),
    dischargePointer() 
    {
      height[source] = network.getVertexCount();
      for (auto it = network.getEdgesIterator(source); it.valid(); it.next()) {
        excess[it.getFinish()] += it.getResidualCapacity();
        it.pushFlow(it.getResidualCapacity());
      }
      for (Tvertex v = 0; v < network.getVertexCount(); ++v) {
        dischargePointer.push_back(network.getEdgesIterator(v));
      }
    }

    Tvertex getSource() const {
      return source;
    }

    Tvertex getSink() const {
      return sink;
    }

    Tflow getExcess(Tvertex v) const {
      return excess[v];
    }

    int getVertexCount() const {
      return network.getVertexCount();
    }

protected:
   /// 0 - no push
   /// 1 - part of push
   /// 2 - full push
   int push(Titerator edge) {
    if (height[edge.getStart()] != (height[edge.getFinish()] + 1)) {
      return 0;
    }
    int result;
    Tflow residualFlow;
    if (excess[edge.getStart()] < edge.getResidualCapacity()) {
      residualFlow = excess[edge.getStart()];
      result = 1;
    } else {
      residualFlow = edge.getResidualCapacity();
      result = 2;
    }
    if (residualFlow <= 0) {
      return 0;
    }
    edge.pushFlow(residualFlow);
    excess[edge.getStart()] -= residualFlow;
    excess[edge.getFinish()] += residualFlow;
    return result;
  }

  /// was relabel or not
  bool relabel(Tvertex startVertex) {
    if (excess[startVertex] <= 0) {
      return false;
    }
    int minLevel = -1;
    for (auto it = network.getEdgesIterator(startVertex); it.valid(); it.next()) {
      if (it.getResidualCapacity() <= 0) {
        continue;
      } 
      int currentLevel = height[it.getFinish()] + 1;
      if (currentLevel <= height[startVertex]) {
        continue;
      }
      if ((minLevel == -1) || (currentLevel < minLevel)) {
        minLevel = currentLevel;
      }
    }
    if (minLevel == -1) {
      return false;
    }
    height[startVertex] = minLevel;
    dischargePointer[startVertex] = network.getEdgesIterator(startVertex);
    return true;
  }

  /// is it to be repeated
  bool discharge(Tvertex startVertex) {
    while (dischargePointer[startVertex].valid()) {
      if (excess[startVertex] <= 0) {
        break;
      }
      if (push(dischargePointer[startVertex]) != 1) {
        dischargePointer[startVertex].next();
      }
    }
    return relabel(startVertex);
  }

  // was relable
  bool maximalDischarge(Tvertex startVertex) {
    int counter = 0;
    while (discharge(startVertex)) {
      ++counter;
    } 
    return counter > 0;
  }
};

template<typename Tvertex, typename Tflow>
class GoldbergAlgorithm: public preflowMaximalFlowAlgorithm<Tvertex, Tflow> {
  using preflowMaximalFlowAlgorithm<Tvertex, Tflow>::preflowMaximalFlowAlgorithm;

public:
  Tflow getMaximalFlow() {
    std::list<Tvertex> que;
    for (Tvertex v = 0; v < this->getVertexCount(); ++v) {
      if (v == this->getSource()) {
        continue;
      }
      if (v == this->getSink()) {
        continue;
      }
      que.push_back(v);
    }
    auto current = que.begin();
    while (current != que.end()) {
      Tvertex v = *current;
      if (this->maximalDischarge(v)) {
        que.erase(current);
        que.push_front(v);
        current = que.begin();
      } else {
        ++current;
      }
    }
    return this->getExcess(this->getSink());
  }
};

int main() {
  return 0;
}
