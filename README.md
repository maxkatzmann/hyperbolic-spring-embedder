# Hyperbolic Spring Embedder

Given the adjacency information of a network, the hyperbolic spring embedder
tries to assign each vertex a hyperbolic coordinate such that adjacent vertices
are placed close to each other, while non-adjacent vertices are positioned
farther apart.  This is done by simulating physical forces between the vertices.

For more details see our publication: https://doi.org/10.4230/LIPIcs.SEA.2021.22
> **Force-Directed Embedding of Scale-Free Networks in the Hyperbolic Plane**
> by Thomas Bläsius, Tobias Friedrich, Maximilian Katzmann

The paper will be presented at the [19th Symposium on Experimental Algorithms
(SEA 2021)](https://sea2021.i3s.unice.fr).

## Usage

The project is built using [Bazel](https://bazel.build).  Consequently, running
the embedder should be as simple as navigating to the root directory of the
project and executing the following command
```
bazel run -c opt main -- --graph /path/to/edgelist.txt --embed /path/to/output/embedding.txt --draw /path/to/output/drawing.svg
```

## Flags

All supported flags can be listed by running
```
bazel run -c opt main -- --help
```

The most important flags are:
- `graph`: The path to the edge list of the graph that should be embedded.  In
  it each line is expected to contain two strings that denote vertex labels.
- `embed`: The path to the file in which the resulting should be stored.  In it
  each line will contain the label of a vertex followed by its hyperbolic
  coordinates (radius and angle in radians).
- `draw`: The path to the SVG file that should contain the visualized embedding.
