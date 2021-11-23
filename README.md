# pagerank-julia

## Term project for 21-241 Matrices and Linear Transformations

## Siyuan "Daniel" Chen, Michael Zhou

### TODO

- Add more, simpler graphs
  - Write in raw format
  - Translate using `generate.jl`
- Need larger test sets
- Visualize graphs for presentation
- ~~Generate matrix base on nodes~~
- ~~Implement damping~~
- ~~Display PageRank results~~
- What's next? - HITS
- Do some reading: once we have the matrix, how/what do we do with eigenstuff?
  - Apparently simply doing eigen would not work
  - Need to use damping factor since has cycles and sinkholes
- How should we present PageRank and HITS
  - Is adding HITS enough for complexity score?
- Search for datasets
- Optimization for sparse networks?
- Optimization for network of networks?

### Data used

[name](https://snap.stanford.edu/data/soc-Epinions1.html)
[name](https://snap.stanford.edu/data/soc-Slashdot0811.html)
[name](https://snap.stanford.edu/data/wiki-Vote.html)
[name](https://snap.stanford.edu/data/email-EuAll.html)
[name](https://snap.stanford.edu/data/wiki-Talk.html)
[name](https://snap.stanford.edu/data/web-Google.html)
[name](https://snap.stanford.edu/data/p2p-Gnutella04.html)
[name](https://snap.stanford.edu/data/p2p-Gnutella06.html)
[name](https://snap.stanford.edu/data/p2p-Gnutella09.html)
