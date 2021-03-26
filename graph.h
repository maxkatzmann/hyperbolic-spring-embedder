#pragma once

#include <string>
#include <unordered_map>
#include <vector>

using std::ostream;
using std::string;
using std::unordered_map;
using std::vector;

class Graph {
 public:
  /**
   * The number of vertices in the graph.
   */
  int n;

  /**
   * The edges in the graph.
   */
  vector<vector<int>> edges;

  /**
   * The labels of the nodes.
   */
  vector<std::string> labels;

  /**
   * Yields a vector that has at position i the number of nodes with
   * degree i.
   */
  vector<int> degHisto() const;

  /**
   * Constructor for Graph.
   *
   * The graph will have n vertices and no edges.
   *
   * @param n The number of nodes of the graph.
   */
  Graph(int n) {
    this->n = n;
    edges.resize(n);
    for (int i = 0; i < n; ++i) {
      labels.push_back(std::to_string(i));
    }
  }

  Graph(const Graph& G);
  virtual ~Graph() {}

  /**
   * Reads a graph from an edge list stored in a file.
   *
   * @param  filename      The name of the file where the graph is to
   *                       be read from.
   * @param  label_to_node Maps a node index to each labeled node.
   * @return               The resulting graph.
   */
  static Graph fromFile(const string& filename,
                        unordered_map<std::string, int>* label_to_node);

  /**
   * Changes the number of nodes.
   *
   * Adds zero-degree vertices, if n is larger than the number of vertices
   * in the graph.
   *
   * Deletes the last vertices in the graph if n is smaller than the
   * number of vertices in the graph.
   *
   * Warning: Does not delete the edges to the deleted nodes!
   *
   * @param n The number of nodes that graph has after the resizing.
   */
  void resize(int n);

  /**
   * Deletes the edge from node1 to node2.
   *
   * Runtime linear in the degree of node1 + node2.
   *
   * @param node1 One node of the edge.
   * @param node2 The other node of the edge.
   */
  void deleteEdge(int node1, int node2);

  /**
   * Get the subgraph of a graph containing only the vertices that are
   * marked as 'true' in the passed vertices vector.
   *
   * @param  vertices A vector of bool values indicating for each vertex
   *                  whether it is included in the subgraph or not.
   * @return          The resulting subgraph.
   */
  Graph subgraph(vector<bool> vertices);

  /**
   * Finds the subgraph of the largest connected component in the graph.
   * @return The corresponding subgraph.
   */
  Graph giantSubgraph();

  /**
   * A simple subgraph, i.e. no self-loops or double edges.
   * @return [description]
   */
  Graph simpleSubgraph();

  /**
   * Calculates the average degree of the graph.
   * @return The average degree.
   */
  double averageDegree() const;

  /**
   * The number of edges in the graph.
   * @return The number of edges.
   */
  int numEdges() const;

  /**
   * The minimum degree of the graph.
   * @return The minimum degree.
   */
  int minDegree();

  /**
   * The maximum degree of the graph.
   * @return The maximum degree.
   */
  int maxDegree();

  /**
   * Estimates the power-law exponent of the graph. Uses a cumulative
   * degree sequence (i.e., a rank sequence)
   * @param  xmin [description]
   * @return      The estimated power-law exponent of the graph.
   */
  double powerLawExponent(int* xmin = nullptr) const;

  /**
   * The sizes of the components in the graph.
   * @return	A vector containin the sizes of the components of the graph.
   */
  vector<int> components() const;

  /**
   * Determines the giant component of a graph.
   * @param components	A vector containing the sizes of the components.
   * @return 				The bool vector indicating which
   * vertices belong to the giant component.
   */
  vector<bool> giant(vector<int> components) const;

  /**
   * The component histogram of all components except the giant one.
   * @return       A vector where the value at index i denotes the number
   *               of components of size i.
   */
  vector<int> smallComponentHisto(vector<int> components, vector<bool> giant);

  /**
   * The number of vertices in the giant component.
   * @param  giant The bool vector indicating which vertices belong to
   *               the giant component.
   * @return       The number of vertices in the giant component.
   */
  int sizeOfGiant(vector<bool> giant) const;

  /**
   * An approximation of the clustering coefficient.
   * @param  eps   [description]
   * @param  delta [description]
   * @return       An approximated clustering coefficient.
   */
  double approximateClusterCoeff(double eps, double delta);

  /**
   * Determines the clustering coefficient of the graph.
   * @return The clustering coefficient of the graph.
   */
  double clusterCoeff();

  /**
   * Determines the distance between u and v in the graph.
   * @param  u             A vertex in the graph.
   * @param  v             Another vertex in the graph.
   * @param  pointsVisited An array containing the vertices visited
   *                       when determining the distance between u
   *                       and v.
   * @return               The distance between u and v.
   */
  int distance(int u, int v, int* pointsVisited = NULL) const;

  /**
   * A vector containing the distances from v to all other vertices.
   * @param v			The vertex whose distance to all other vertices
   * is to be determined.
   * @param parent	An array containg the parent of each vertex obtained
   *               	when determining the distances between v and the
   *               	remaining vertices in the graph.
   * @return          A vector where the value at index i is the distance
   *                  between vertex v and the node i.
   */
  vector<int> distances(int v, vector<int>* parent = NULL);

  /**
   * Determines the diameter of the graph.
   * @param  nodes_on_diameter A bool vector indicating which vertices belong to
   * the path representing the diameter of the graph.
   * @return The diameter of the graph.
   */
  int diameter(vector<bool>* nodes_on_diameter);

  /**
   * A vector where the value at index i indicates the number of nodes
   * that have distance i to vertex v.
   * @param v		The vertex whose distance histogram is to be determined.
   * @return The distance histogram of vertex v.
   */
  vector<int> distanceHisto(int v);

  /**
   * Samples the distance histogram of the whole graph.
   * @param numsamples The number of samples used during the sampling.
   * @return A vector where the value at index i contains an approximation
   *         of how often the distance i occured between two vertices in
   *         the graph.
   */
  vector<double> sampleDistanceHisto(int numsamples);

  /**
   * An approximation of the average distance in the graph.
   * @param  eps                [description]
   * @param  delta              [description]
   * @param  diameterUpperBound [description]
   * @return                    [description]
   */
  double approximateAverageDistance(double eps, double delta,
                                    int diameterUpperBound);

  /**
   * Samples the average distance of the graph.
   * @param  numSamples The number of samples to use during the sampling.
   * @return            An approximation of the average distance in the
   *                    graph.
   */
  double sampleAverageDistance(int numSamples);

  /**
   * Samples an interval containing a lower and upper bound on the diameter.
   * @param numSamples The number of samples to use during the sampling.
   * @param upperBound The resulting upper bound on the diameter.
   * @param lowerBound The resulting lower bound on the diameter.
   */
  void diameterIntervalSampler(int numSamples, int& upperBound,
                               int& lowerBound);

  int push(int src);

  void pushpull(int src, double alpha1, double alpha2, int& time1, int& time2);

  double asynch_pushpull(int src, double alpha = 0.0);

  void asynch_pushpull(int src, double alpha1, double alpha2, double& time1,
                       double& time2);

  int pushpull(int src, double alpha, vector<int>* evol = NULL);

  /**
   * Determines whether the vertices u and v are adjacent.
   * @param  u A vertex in the graph.
   * @param  v Another vertex in the graph.
   * @return   true if u and v are adjacent. false if not.
   */
  bool adjacent(int u, int v) const;

  /**
   * Permutes the nodes on the graph to obtain an isomorphic graph.
   * Discards each node independently with probability 1-q.
   * @param q The inverse probability of discarding a node.
   */
  vector<int> permuteAndSubsampleGraph(double q);

  /**
   * Permutes nodes according to the given permutation.
   * If permutation[i] == -1 node i is discarded.
   * @param permutation [description]
   */
  void permuteGraph(const vector<int>& permutation);

  /**
   * Subsamples the graph's edges, keeping each edge with probability p.
   * @param  p The probability that an edge is kept.
   * @return   The subsampled graph.
   */
  Graph subsample(double p);

  /**
   * Reorders the adjacency lists s.t. they are ordered according to the
   * given permutation of nodes. Does not change the graph structure.
   *
   * E.g. if the permutation is [4,3,5,1,2] then a node with neighbors
   * 1,2,3 stores the neighbors as [3,1,2].
   * @param permutation The permutation that the node orders should be
   *                    adapted to.
   */
  void ReorderEdges(const vector<int>& permutation);

  /**
   * Determines whether the graph is directed or not.
   * @return true if and only if there exists an edge u-> but not v->u.
   */
  bool isDirected();

  /**
   * Reorders the nodes by decreasing degree. edges[0] is the node with
   * the largest degree. new_permutation will contain the reordering of
   * the nodes.
   * permutation[i] = j means that old node i is now node j.
   * @param new_permutation [description]
   */
  virtual void sortByDegrees(vector<int>* new_permutation);

  /**
   * Saves this graph to a file with the given name. Saves each edge
   * only once, i.e. consideres the graph to be undirected. Ignores
   * all edges (u, v) where either ignore_nodes[u] == true or
   * ignore_nodes[v] == true.
   * @param filename     The file that the graph should be written to.
   * @param ignore_nodes The nodes to not consider in the written graph.
   * @param use_labels   Denotes whether the indices of the nodes should
   *                     be used or their original labels.
   */
  virtual void printToFile(const char* filename,
                           const vector<bool>* ignore_nodes = NULL,
                           bool use_labels = true);

  /**
   * Determines a random edge. The two endpoints of the edge are passed
   * through.
   * @param u One endpoint of the random edge.
   * @param v The other endpoint of the random edge.
   */
  void getRandomEdge(int* u, int* v);

  /**
   * Swaps (w.h.p.) all edges at random, effectifely creating a new graph
   * with the same degree distribution but independent edge probabilities.
   * @return The graph with the swapped edges.
   */
  Graph randomSwap();

  /**
   * Counts the number of common neighbors of nodes u and v.
   * @param  u One vertex in the graph.
   * @param  v Another vertex in the graph.
   * @return   The number of common neighbors of u and v.
   */
  int commonNeighbors(int u, int v);

  /**
   * ignore_nodes indicates which nodes are indistinguishable from their
   * neighborhoods. ignore_nodes[i] is true if there is a node j in the
   * graph that has the exact same neighborhood as i, and ignore_nodes[j] =
   * false.
   * @param ignore_nodes The bool vector indicating whether node i should be
   *                     ignored or not.
   */
  void indistinguishableNodes(vector<bool>* ignore_nodes);

  /**
   * Performes a lexicographic breadth-first search on the graph.
   * @return The lexicographic ordering of the nodes.
   */
  vector<int> LexBFS();

  /**
   * Tests whether the graph is chordal, i.e. all cycles of length > 3
   * have a shortcut.
   * @return true if the graph is chordal, false if not.
   */
  bool isChordal();
};

ostream& operator<<(ostream& out, Graph& g);
