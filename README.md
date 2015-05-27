<snippet>
  <content><![CDATA[
# ${1:Markovian}

A simple Markov Network implementation that supports brute force partitioning, variable elimination, and Loopy Belief Propagation.

## Usage

Network input format specified [here.](http://www.hlt.utdallas.edu/~vgogate/uai14-competition/modelformat.html)

To perform a brute force calculation of a given network's partition function do:

```
java -jar brute-force-partition.jar [path/to/network/file]
```

Note that the above is impractical for even medium sized networks. To more efficiently calculate a network's partition function using variable elimination do:
 
```
java -jar variable-elimination-partition.jar [path/to/network/file]
```

To compute variable marginals using the Loopy Belief Propagation algorithm do:

```
java -jar loopy-bp.jar [path/to/network/file]
```


