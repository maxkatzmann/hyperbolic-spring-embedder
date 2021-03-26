#include "graph.h"

#include <glog/logging.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <unordered_set>

#include "powerlawCommon.h"
#include "random.h"

#define FOR(i, n) for (int(i) = 0; (i) < (n); (i)++)
#define FORB(i, a, n) for (int(i) = (a); (i) < (n); (i)++)
#define FOREACH(it, c) \
  for (__typeof((c).begin()) it = (c).begin(); it != (c).end(); ++it)
#define PB push_back
#define MP make_pair

using std::cout;
using std::endl;
using std::map;
using std::max;
using std::min;
using std::queue;
using std::set;
using std::stack;
using std::swap;

Graph::Graph(const Graph& G) {
  n = G.n;
  edges.resize(n);
  FOR(i, n)
  FOR(j, G.edges[i].size())
  edges[i].PB(G.edges[i][j]);
  labels = G.labels;
}

Graph Graph::fromFile(const string& filename,
                      unordered_map<string, int>* label_to_node) {
  Graph g(0);
  std::ifstream file;
  file.open(filename.c_str());
  CHECK(file.good()) << "Maybe the file " << filename << " doesn't exist?";
  label_to_node->clear();
  g.labels.clear();

  // read all nodes
  int current_node = 0;
  string current_label;
  while (file >> current_label) {
    if ((file >> std::ws).peek() == '#') {
      // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }

    // ignore directory structure
    current_label = current_label.substr(0, current_label.find("/"));

    if (label_to_node->count(current_label) == 0) {
      // new node
      (*label_to_node)[current_label] = current_node;
      ++current_node;
      g.labels.PB(current_label);
    }
  }

  file.close();

  // now construct graph
  g.resize(current_node);

  file.open(filename.c_str());
  string u, v;

  while (file >> u && file >> v) {
    if ((file >> std::ws).peek() == '#') {  // ignore rest of the line
      std::getline(file, current_label);
      continue;
    }

    u = u.substr(0, u.find("/"));
    v = v.substr(0, v.find("/"));
    if ((*label_to_node).find(u) != label_to_node->end() &&
        (*label_to_node).find(v) != label_to_node->end()) {
      g.edges[(*label_to_node)[u]].PB((*label_to_node)[v]);
      g.edges[(*label_to_node)[v]].PB((*label_to_node)[u]);
      //      LOG(INFO) << "u=" << u << "               v=" << v;
      //      LOG(INFO) << (*label_to_node)[u] << "-" << (*label_to_node)[v];
    }
  }

  file.close();

  return g;
}

void Graph::resize(int n) {
  int old_n = this->n;
  this->n = n;
  edges.resize(n);
  labels.resize(n, "");

  FORB(i, old_n, n) if (labels[i].empty()) labels[i] = std::to_string(i);
}

int Graph::push(int src) {
  vector<int> infoset;
  vector<bool> informed(n, false);
  infoset.push_back(src);
  informed[src] = true;
  int numinformed = 1;

  int T = 0;
  vector<int> newly;
  while (numinformed < n) {
    newly.clear();

    FOREACH(it, infoset) {
      int i = *it;
      int j = edges[i][randint((int)edges[i].size())];
      if (!informed[j]) {
        informed[j] = true;
        numinformed++;
        newly.push_back(j);
      }
    }

    FOREACH(it, newly) infoset.push_back(*it);

    ++T;
  }

  return T;
}

void Graph::pushpull(int src, double alpha1, double alpha2, int& time1,
                     int& time2) {
  vector<bool> informed(n, false);
  vector<bool> newinf(n, false);
  informed[src] = true;
  newinf[src] = true;
  int numinformed = 1;

  int T = 0;
  bool reached1 = false, reached2 = false;
  while (!reached1 || !reached2) {
    if (numinformed >= (1. - alpha1) * n && !reached1) {
      reached1 = true;
      time1 = T;
    }
    if (numinformed >= (1. - alpha2) * n && !reached2) {
      reached2 = true;
      time2 = T;
    }

    FOR(i, n) {
      if (informed[i]) {
        // push
        int j = edges[i][randint((int)edges[i].size())];
        if (!newinf[j]) {
          newinf[j] = true;
          numinformed++;
        }
      } else if (!newinf[i]) {
        // pull
        int j = edges[i][randint((int)edges[i].size())];
        if (informed[j]) {
          newinf[i] = true;
          numinformed++;
        }
      }
    }
    FOR(i, n) informed[i] = newinf[i];

    ++T;
  }
}

double Graph::asynch_pushpull(int src, double alpha) {
  vector<bool> informed(n, false);
  informed[src] = true;
  int numinformed = 1;

  int T = 0;
  while (numinformed < (1. - alpha) * n) {
    int i = randint(n);
    int j = edges[i][randint((int)edges[i].size())];
    if (informed[i]) {
      if (!informed[j]) {
        informed[j] = true;
        numinformed++;
      }
    } else {
      if (informed[j]) {
        informed[i] = true;
        numinformed++;
      }
    }
    T++;
  }

  return ((double)T) / ((double)n);
}

void Graph::asynch_pushpull(int src, double alpha1, double alpha2,
                            double& time1, double& time2) {
  vector<bool> informed(n, false);
  informed[src] = true;
  int numinformed = 1;

  int T = 0;
  bool reached1 = false, reached2 = false;
  while (!reached1 || !reached2) {
    if (numinformed >= (1. - alpha1) * n && !reached1) {
      reached1 = true;
      time1 = ((double)T) / ((double)n);
    }
    if (numinformed >= (1. - alpha2) * n && !reached2) {
      reached2 = true;
      time2 = ((double)T) / ((double)n);
    }

    int i = randint(n);
    int j = edges[i][randint((int)edges[i].size())];
    if (informed[i]) {
      if (!informed[j]) {
        informed[j] = true;
        numinformed++;
      }
    } else {
      if (informed[j]) {
        informed[i] = true;
        numinformed++;
      }
    }
    T++;
  }
}

int Graph::pushpull(int src, double alpha, vector<int>* evol) {
  vector<bool> informed(n, false);
  vector<bool> newinf(n, false);
  informed[src] = true;
  newinf[src] = true;
  int numinformed = 1;

  if (evol != NULL) evol->clear();

  int T = 0;
  while (numinformed < (1. - alpha) * n) {
    int thisround = 0;
    FOR(i, n) {
      if (informed[i]) {
        // push
        int j = edges[i][randint((int)edges[i].size())];
        if (!newinf[j]) {
          newinf[j] = true;
          numinformed++;
          thisround++;
        }
      } else {
        // pull
        int j = edges[i][randint((int)edges[i].size())];
        if (informed[j]) {
          newinf[i] = true;
          numinformed++;
          thisround++;
        }
      }
    }
    FOR(i, n) informed[i] = newinf[i];

    if (evol != NULL) evol->push_back(thisround);
    ++T;
    std::cout << "Round " << T << " informed " << thisround << " nodes."
              << std::endl;
  }

  return T;
}

int Graph::minDegree() {
  int res = (int)edges[0].size();
  FOR(i, n) res = min(res, (int)edges[i].size());
  return res;
}

int Graph::maxDegree() {
  int res = (int)edges[0].size();
  FOR(i, n) res = max(res, (int)edges[i].size());
  return res;
}

Graph Graph::subgraph(vector<bool> vertices) {
  vector<int> index(n + 1);
  index[0] = 0;
  FORB(i, 1, n + 1) index[i] = index[i - 1] + (vertices[i - 1] ? 1 : 0);
  int m = index[n];
  Graph G(m);
  FOR(i, n) {
    if (vertices[i]) {
      G.labels[index[i]] = labels[i];
      FOR(j, (int)edges[i].size()) {
        int jj = edges[i][j];
        if (vertices[jj]) G.edges[index[i]].PB(index[jj]);
      }
    }
  }
  return G;
}

int Graph::distance(int u, int v, int* pointsVisited) const {
  if (u == v) return 0;
  int ind[2];
  ind[0] = u;
  ind[1] = v;
  set<int> reached[2];
  set<int> level[2];
  FOR(t, 2) {
    reached[t].insert(ind[t]);
    level[t].insert(ind[t]);
  }
  int dist = 0;
  while (level[0].size() != 0 && level[1].size() != 0) {
    int t = (level[0].size() < level[1].size()) ? 0 : 1;
    dist++;
    set<int> newlvl;
    FOREACH(it, level[t])
    FOREACH(nei, edges[*it])
    if (reached[t].find(*nei) == reached[t].end()) {
      reached[t].insert(*nei);
      newlvl.insert(*nei);
      if (reached[1 - t].find(*nei) != reached[1 - t].end()) {
        if (pointsVisited != NULL)
          *pointsVisited = (int)(reached[0].size() + reached[1].size());
        return dist;
      }
    }
    level[t] = newlvl;
  }

  return n + 1;
}

vector<int> Graph::distances(int v, vector<int>* parent) {
  vector<int> dist(n, -1);
  dist[v] = 0;
  if (parent != NULL) (*parent)[v] = -1;

  queue<int> q;
  q.push(v);
  int max_dist = 0;

  while (!q.empty()) {
    int curr = q.front();
    q.pop();

    vector<int>::iterator I;
    for (I = edges[curr].begin(); I != edges[curr].end(); ++I)
      if (dist[*I] == -1) {
        q.push(*I);
        dist[*I] = 1 + dist[curr];
        if (parent != NULL) (*parent)[*I] = curr;

        max_dist = max(max_dist, dist[*I]);
      }
  }

  return dist;
}

int Graph::diameter(vector<bool>* nodes_on_diameter) {
  int totalmaxdist = 0;
  nodes_on_diameter->clear();
  nodes_on_diameter->resize(n);
  FOR(i, n) {
    vector<int> parent(n);
    vector<int> dists = distances(i, &parent);

    FOR(j, n) {
      if (totalmaxdist < dists[j]) {
        totalmaxdist = dists[j];

        std::fill(nodes_on_diameter->begin(), nodes_on_diameter->end(), false);
        int curr = j;
        while (curr != i) {
          (*nodes_on_diameter)[curr] = true;
          curr = parent[curr];
        }
        (*nodes_on_diameter)[i] = true;
      }
    }
  }

  // LOG(INFO) << "Diam:" << totalmaxdist;
  return totalmaxdist;
}

void Graph::diameterIntervalSampler(int numSamples, int& upperBound,
                                    int& lowerBound) {
  upperBound = n;
  lowerBound = 0;
  FOR(t, numSamples) {
    int start = randint(n);
    if (t == 0) FOR(i, n) if (edges[i].size() > edges[start].size()) start = i;

    vector<int> parent(n);
    vector<int> dist = distances(start, &parent);
    int worst = 0;
    FOR(i, n) if (dist[i] > dist[worst]) worst = i;
    upperBound = min(upperBound, 2 * dist[worst]);
    lowerBound = max(lowerBound, dist[worst]);

    dist = distances(worst);
    int maxdist = 0;
    FOR(i, n) maxdist = max(maxdist, dist[i]);
    upperBound = min(upperBound, 2 * maxdist);
    lowerBound = max(lowerBound, maxdist);

    // todo: suche besseren Mittelpunkt in dem Baum!
    vector<int> degree(n, 0);
    vector<int> reports(n, 0);
    vector<int> largest(n, 0);
    vector<int> largestindex(n, 0);
    vector<int> seclargest(n, 0);

    FOR(i, n) if (i != start) degree[parent[i]]++;
    FOR(ii, n) {
      int i = ii;
      while (degree[i] == reports[i] && i != start) {
        // report at parent:
        int p = parent[i];
        if (reports[p] == 0) {
          largest[p] = largest[i] + 1;
          largestindex[p] = i;
        } else if (reports[p] == 1) {
          seclargest[p] = largest[i] + 1;
          if (seclargest[p] > largest[p]) {
            swap(largest[p], seclargest[p]);
            largestindex[p] = i;
          }
        } else {
          int len = largest[i] + 1;
          if (len >= largest[p]) {
            seclargest[p] = largest[p];
            largest[p] = len;
            largestindex[p] = i;
          } else if (len >= seclargest[p])
            seclargest[p] = len;
        }
        reports[p]++;

        i = p;
      }
    }

    int v = start;
    int over = 0;
    while (max(largest[v] + 1, over) <
           max(max(seclargest[v] + 2, over + 1), largest[v])) {
      over = max(seclargest[v] + 2, over + 1);
      v = largestindex[v];
    }
    if (v != start) {
      vector<int> dist2 = distances(v);
      maxdist = 0;
      FOR(i, n) maxdist = max(maxdist, dist2[i]);
      upperBound = min(upperBound, 2 * maxdist);
      lowerBound = max(lowerBound, maxdist);
    }
  }
}

vector<int> Graph::distanceHisto(int v) {
  vector<int> dist = distances(v);

  int max_dist = 0;
  FOR(i, n) max_dist = max(max_dist, dist[i]);

  vector<int> histo(max_dist + 1, 0);
  for (int i = 0; i < n; ++i)
    if (dist[i] != -1) ++histo[dist[i]];

  return histo;
}

vector<double> Graph::sampleDistanceHisto(int numsamples) {
  vector<double> histo;
  FOR(T, numsamples) {
    int v = randint(n);
    vector<int> histoV = distanceHisto(v);
    if (histoV.size() > histo.size()) {
      int h = (int)histo.size();
      histo.resize(histoV.size());
      for (; h < (int)histo.size(); h++) histo[h] = 0.;
    }
    FOR(h, (int)histoV.size())
    histo[h] += histoV[h];
  }
  FOR(h, (int)histo.size()) histo[h] /= numsamples;
  return histo;
}

double Graph::approximateAverageDistance(double eps, double delta,
                                         int diameterUpperBound) {
  int numsamples = max(1, (int)ceil(0.5 * log(2. / delta) * diameterUpperBound *
                                    diameterUpperBound / (eps * eps)));
  double sum = 0.;
  int visited = 0;
  FOR(T, numsamples) {
    int numVisited;
    sum += distance(randint(n), randint(n), &numVisited);
    visited += numVisited;
  }
  return sum / numsamples;
}

double Graph::sampleAverageDistance(int numSamples) {
  double sum = 0.;
  int visited = 0;
  FOR(T, numSamples) {
    int numVisited;
    sum += distance(randint(n), randint(n), &numVisited);
    visited += numVisited;
  }
  return sum / numSamples;
}

ostream& operator<<(ostream& out, Graph& g) {
  for (int i = 0; i < g.n; ++i) {
    out << i << ": ";
    for (vector<int>::iterator I = g.edges[i].begin(); I != g.edges[i].end();
         ++I)
      out << *I << " ";
    out << endl;
  }

  return out;
}

double Graph::averageDegree() const { return 2.0 * numEdges() / ((double)n); }

int Graph::numEdges() const {
  int m = 0;
  FOR(i, n) m += (int)edges[i].size();
  return m / 2;
}

/*void Graph::mark(int i, int T, vector<int> & comp) {
 if (comp[i] == -1) {
 comp[i] = T;
 FOREACH(it,edges[i])
 mark(*it,T,comp);
 }
 }*/

vector<int> Graph::components() const {
  vector<int> comp(n, -1);
  int T = 0;
  FOR(i, n)
  if (comp[i] == -1) {
    stack<int> S;
    S.push(i);
    while (!S.empty()) {
      int v = S.top();
      S.pop();
      if (comp[v] == -1) {
        comp[v] = T;
        FOREACH(it, edges[v])
        S.push(*it);
      }
    }

    T++;
  }
  return comp;
}

vector<bool> Graph::giant(vector<int> components) const {
  vector<int> size(n, 0);
  FOR(i, n) size[components[i]]++;
  int maxcomp = 0;
  FOR(i, n)
  if (size[i] > size[maxcomp]) maxcomp = i;
  vector<bool> giant(n);
  FOR(i, n)
  giant[i] = (components[i] == maxcomp);
  return giant;
}

vector<int> Graph::smallComponentHisto(vector<int> components,
                                       vector<bool> giant) {
  vector<int> size(n, 0);
  int maxsize = 0;
  FOR(i, n)
  if (!giant[i]) {
    size[components[i]]++;
    maxsize = max(maxsize, size[components[i]]);
  }
  vector<int> histo(maxsize + 1, 0);
  FOR(i, n) histo[size[i]]++;
  histo[0] = 0;
  return histo;
}

int Graph::sizeOfGiant(vector<bool> giant) const {
  int size = 0;
  FOR(i, n) size += (giant[i]) ? 1 : 0;
  return size;
}

double Graph::clusterCoeff() {
  int numclosed = 0;
  int numopen = 0;
  FOR(i, n)
  FOR(j, (int)edges[i].size())
  FOR(k, j) {
    int jj = edges[i][j], kk = edges[i][k];
    if (edges[jj].size() > edges[kk].size()) swap(jj, kk);
    bool triangle = false;
    FOR(l, (int)edges[jj].size())
    if (edges[jj][l] == kk) {
      triangle = true;
      break;
    }
    if (triangle)
      numclosed++;
    else
      numopen++;
  }
  numclosed /= 3;
  return ((double)numclosed) / ((double)numclosed + numopen);
}

bool Graph::adjacent(int u, int v) const {
  if (edges[u].size() > edges[v].size()) return adjacent(v, u);
  FOREACH(w, edges[u])
  if (*w == v) return true;
  return false;
}

double Graph::approximateClusterCoeff(double eps, double delta) {
  vector<int> twonodes;
  FOR(i, n) if (edges[i].size() >= 2) twonodes.PB(i);
  int k = (int)ceil(log(2. / delta) / (2. * eps * eps));
  int l = 0;
  if (twonodes.size() == 0) return 0.;
  FOR(i, k) {
    int r = randint((int)twonodes.size());
    int j = twonodes[r];
    int u = edges[j][randint((int)edges[j].size())];
    int w;
    do {
      w = edges[j][randint((int)edges[j].size())];
    } while (u == w);
    if (adjacent(u, w)) l++;
  }

  return ((double)l) / ((double)k);
}

vector<int> Graph::degHisto() const {
  int maxDeg = 0;
  FOR(i, n) maxDeg = max(maxDeg, (int)edges[i].size());
  vector<int> histo(maxDeg + 1);
  FOR(i, n) histo[edges[i].size()]++;
  return histo;
}

// double Graph::powerLawExponent(int dmin) const {
//  vector<int> deg = degHisto();
//  int maxDeg = (int)(deg.size()-1);
//
////  if (dmin < 0) {
//  for(int i = maxDeg-1; i>=0; i--) {
//    deg[i] += deg[i+1];
//  }
//  double meanx = 0.;
//  int cnt = 0;
//  FORB(i,dmin,maxDeg) if (deg[i] != 0) cnt++;
//  FORB(i,dmin,maxDeg) if (deg[i] != 0) meanx += log(i);
//  meanx /= cnt;
//  double meany = 0.;
//  FORB(i,dmin,maxDeg) if (deg[i] != 0) meany += log(deg[i]);
//  meany /= cnt;
//  double num = 0.;
//  FORB(i,dmin,maxDeg) if (deg[i] != 0) num += (log(i) - meanx)*(log(deg[i]) -
//  meany); double den = 0.; FORB(i,dmin,maxDeg) if (deg[i] != 0) den += (log(i)
//  - meanx)*(log(i) - meanx); return -num/den + 1.; // LEAST SQUARES FIT:
//  http://mathworld.wolfram.com/LeastSquaresFitting.html
//}

double Graph::powerLawExponent(int* xmin) const {
  plfit::VectorType degs;
  FOR(i, n)
  degs.PB(edges[i].size());
  plfit::VectorType results;
  plfit::Powerlaw::SingleFit(degs, results, false, false, 1.5, 0.01, 4);

  if (xmin) (*xmin) = results[1];

  return results[0];
}

Graph Graph::giantSubgraph() {
  vector<int> comp = components();
  vector<bool> gcc = giant(comp);
  return subgraph(gcc);
}

Graph Graph::simpleSubgraph() {
  vector<bool> seen(n, false);
  Graph G(n);
  FOR(i, n) {
    G.labels[i] = labels[i];
    for (int node : edges[i])
      if (!seen[node] && (node != i)) {
        seen[node] = true;
        G.edges[i].PB(node);
      }
    std::fill(seen.begin(), seen.end(), false);
  }
  return G;
}

vector<int> Graph::permuteAndSubsampleGraph(double q) {
  vector<int> permutation(n);
  for (int i = 0, j = 0; i < n; i++) {
    if (randdblpos() < q) {
      permutation[i] = j;
      ++j;
    } else {
      permutation[i] = -1;
    }
  }

  std::random_device rd;
  std::mt19937 g(rd());

  std::shuffle(permutation.begin(), permutation.end(), g);
  permuteGraph(permutation);
  return permutation;
}

void Graph::permuteGraph(const vector<int>& permutation) {
  int num_subsampled_nodes = 0;
  FOR(i, permutation.size())
  if (permutation[i] >= 0) ++num_subsampled_nodes;

  vector<vector<int>> new_edges(num_subsampled_nodes);
  vector<std::string> new_labels(num_subsampled_nodes);
  FOR(i, n) {
    if (permutation[i] >= 0) {
      new_edges[permutation[i]].swap(edges[i]);
      new_labels[permutation[i]] = labels[i];
    }
  }

  FOR(i, num_subsampled_nodes) {
    vector<int> new_edges_i;
    FOR(j, new_edges[i].size())
    if (permutation[new_edges[i][j]] >= 0)
      new_edges_i.PB(permutation[new_edges[i][j]]);
    new_edges[i].swap(new_edges_i);
  }

  edges.swap(new_edges);
  labels.swap(new_labels);
  n = (int)edges.size();
}

Graph Graph::subsample(double p) {
  Graph H(n);

  FOR(i, n) {
    FOR(j, edges[i].size()) {
      if (edges[i][j] <= i) continue;  // Saw this edge already
      if (randdblpos() < p) {
        H.edges[i].PB(edges[i][j]);
        H.edges[edges[i][j]].PB(i);
      }
    }
  }

  return H;
}

bool Graph::isDirected() {
  FOR(i, n)
  FOR(j, edges[i].size())
  if (std::find(edges[edges[i][j]].begin(), edges[edges[i][j]].end(), i) ==
      edges[edges[i][j]].end()) {
    cout << "Didn't find edge back: " << i << "<->" << edges[i][j];
    return true;
  }

  return false;
}

void Graph::sortByDegrees(vector<int>* new_permutation) {
  // calculate the histogram of key frequencies:
  vector<int> count(n, 0);
  FOR(i, n)
  count[edges[i].size()]++;

  // calculate the starting index for each key:
  int total = 0;
  for (int i = n - 1; i >= 0; --i) {
    int oldCount = count[i];
    count[i] = total;
    total += oldCount;
  }

  new_permutation->clear();
  new_permutation->resize(n);

  // copy to output array, preserving order of inputs with equal keys:
  FOR(i, n) {
    int deg = (int)edges[i].size();
    (*new_permutation)[i] = count[deg];
    count[deg] += 1;
  }

  permuteGraph(*new_permutation);
}

void Graph::printToFile(const char* filename, const vector<bool>* ignore_nodes,
                        bool use_labels) {
  std::ofstream file;
  file.open(filename);

  FOR(i, n) {
    for (int j : edges[i]) {
      if (ignore_nodes == NULL ||
          (!ignore_nodes->at(i) && !ignore_nodes->at(j)))
        if (j >= i) {
          if (use_labels)
            file << labels[i] << "\t" << labels[j] << std::endl;
          else
            file << i << "\t" << j << std::endl;
        }
    }
  }

  file.close();
}

void Graph::getRandomEdge(int* u, int* v) {
  *u = -1;
  *v = -1;
  int m = 0;
  FOR(i, n) m += (int)edges[i].size();
  int e = randint(m);

  FOR(i, n) {
    if (e >= edges[i].size()) {
      e -= edges[i].size();
    } else {
      *u = i;
      *v = edges[i][e];
      break;
    }
  }

  CHECK(*u >= 0);
  CHECK(*v >= 0);
}

Graph Graph::randomSwap() {
  Graph g(n);
  g.edges = edges;

  FOR(j, averageDegree() * n * log(averageDegree() * n)) {
    int u1, u2, v1, v2;
    g.getRandomEdge(&u1, &u2);
    g.getRandomEdge(&v1, &v2);

    if (u1 == v1 || u1 == v2 || u2 == v1 || u2 == v2 || g.adjacent(u1, v2) ||
        g.adjacent(v1, u2))
      continue;

    FOR(i, g.edges[u1].size())
    if (g.edges[u1][i] == u2) g.edges[u1][i] = v2;

    FOR(i, g.edges[u2].size())
    if (g.edges[u2][i] == u1) g.edges[u2][i] = v1;

    FOR(i, g.edges[v1].size())
    if (g.edges[v1][i] == v2) g.edges[v1][i] = u2;

    FOR(i, g.edges[v2].size())
    if (g.edges[v2][i] == v1) g.edges[v2][i] = u1;

    CHECK(g.adjacent(u1, v2)) << u1 << " " << u2 << " " << v1 << " " << v2;
    CHECK(g.adjacent(u2, v1)) << u1 << " " << u2 << " " << v1 << " " << v2;
  }

  LOG(INFO) << "Done swapping.";

  return g;
}

int Graph::commonNeighbors(int u, int v) {
  int common_neighbors = 0;
  FOR(i, edges[u].size())
  FOR(j, edges[v].size())
  if (edges[u][i] == edges[v][j]) ++common_neighbors;

  return common_neighbors;
}

void Graph::indistinguishableNodes(vector<bool>* ignore_nodes) {
  ignore_nodes->clear();
  ignore_nodes->resize(n, false);
  // sort all edge vectors
  FOR(i, n) { std::sort(edges[i].begin(), edges[i].end()); }

  FOR(i, n) {
    if (!edges[i].empty() && !(*ignore_nodes)[i]) {
      FORB(j, i + 1, n) {
        if (edges[i].size() == edges[j].size()) {
          bool same = true;
          for (int w : edges[i]) {
            if (w != j && !adjacent(j, w)) {
              same = false;
              break;
            }
          }

          if (same) {
            (*ignore_nodes)[j] = true;
          }
        }
      }
    }
  }
}

void Graph::ReorderEdges(const vector<int>& permutation) {
  vector<vector<int>> new_edges(n);

  FOR(i, n) {
    FOR(j, edges[permutation[i]].size()) {
      new_edges[edges[permutation[i]][j]].PB(permutation[i]);
    }
  }

  edges.swap(new_edges);
}

namespace {
struct lexcomp {
  bool operator()(const vector<int>& lhs, const vector<int>& rhs) const {
    int i = 0;
    while (i < rhs.size()) {
      if (i >= lhs.size() || rhs[i] < lhs[i])
        return false;
      else if (lhs[i] < rhs[i])
        return true;
      ++i;
    }
    return (i != lhs.size());
  }
};
}  // namespace

vector<int> Graph::LexBFS() {
  // Create a ordered mapping that maps labels of vertices (represented as
  // integer lists) to the nodes having these labels.
  std::map<vector<int>, std::unordered_set<int>, lexcomp> queue;

  std::unordered_set<int> all_nodes;
  FOR(i, n)
  all_nodes.insert(i);
  queue[vector<int>()] = all_nodes;

  // We also need to map nodes back to their labels.
  vector<vector<int>> labels(n, vector<int>());

  vector<int> lexbfs;

  std::cout << "LexBFS ordering: ";
  int current_label = 0;

  // Now proceed using a BFS in the order given by the queue
  while (1) {
    // if there are no nodes with that label, remove the key
    while (queue.begin() != queue.end() &&
           queue.begin()->second.begin() == queue.begin()->second.end())
      queue.erase(queue.begin());
    if (queue.empty()) break;

    int node = *((*queue.begin()).second.begin());
    lexbfs.push_back(node);

    //     print label
    std::cout << "Node " << node << " [";
    for (int j : queue.begin()->first) std::cout << j << ", ";
    std::cout << "]" << std::endl;

    // remove the node from the first set
    queue.begin()->second.erase(node);

    // update all neighbors of node
    for (int neighbor : edges[node]) {
      if (queue.find(labels[neighbor]) != queue.end()) {
        int num_erased = (int)queue[labels[neighbor]].erase(neighbor);
        if (num_erased > 0) {
          labels[neighbor].PB(current_label);
          queue[labels[neighbor]].insert(neighbor);
        }
      }
    }

    ++current_label;
  }
  return lexbfs;
}

bool Graph::isChordal() {
  vector<int> lexbfs = LexBFS();
  std::reverse(lexbfs.begin(), lexbfs.end());
  // construct inverse permutation, i.e. node i gets mapped to position
  // lexbfs_inv[i].
  vector<int> lexbfs_inv(n, -1);
  FOR(i, n)
  lexbfs_inv[lexbfs[i]] = i;
  ReorderEdges(lexbfs);

  FOR(i, n) {
    int node = lexbfs[i];
    std::cout << "Node " << node << std::endl;
    // Find out the neighbor of node that appears first after node in the lexbfs
    // ordering
    int nextNeighborIndex = 0;
    while (nextNeighborIndex < edges[node].size() &&
           lexbfs_inv[edges[node][nextNeighborIndex]] < lexbfs_inv[node])
      ++nextNeighborIndex;

    if (nextNeighborIndex == edges[node].size()) continue;

    int nextNeighbor = edges[node][nextNeighborIndex];
    std::cout << "Next neighbor: " << nextNeighbor << std::endl;

    // Test whether all subsequent edges of node are also edges of nextNeighbor
    int neighborIndex = nextNeighborIndex + 1;
    int otherNeighborIndex = 0;
    while (neighborIndex < edges[node].size()) {
      int neighbor = edges[node][neighborIndex];
      std::cout << "Checking neighbor " << neighbor << std::endl;

      if (edges[nextNeighbor].empty()) return false;
      while (edges[nextNeighbor][otherNeighborIndex] != neighbor) {
        ++otherNeighborIndex;
        if (otherNeighborIndex >= edges[nextNeighbor].size()) return false;
      }

      ++neighborIndex;
    }
  }

  return true;
}

void Graph::deleteEdge(int node1, int node2) {
  edges[node1].erase(
      std::remove(edges[node1].begin(), edges[node1].end(), node2),
      edges[node1].end());
  edges[node2].erase(
      std::remove(edges[node2].begin(), edges[node2].end(), node1),
      edges[node2].end());
}
