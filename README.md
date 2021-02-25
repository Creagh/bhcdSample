# Bayesian Hierarchical Community Discovery

## Summary

This repository contains a sample of work, for an ongoing project, which includes an implementation of a fully Bayesian Hierarchical Community Discovery model for networks. This includes multiple model variaitions, aimed at learning clusters within a network (or graph). These models are based on the work by Clauset et al. [1]. 

The code framework is implemented in Java, and interfaces with the probabilistic programming language, [Blang](https://www.stat.ubc.ca/~bouchard/blang/).

## Examples

The main codebase can be found in `src\main\java\bhcd`.

`BasicDendrogram.java` contains the class for a dendrogram structure (rooted binary tree with labelled leaves), which is used to represent the hieararchical clustering structure. 

The method `interchange(BasicTreeNode node, BasicTreeNode child)` performs a Nearest-Neighbour Interchange move, permuting the structure of the dendrogram, used within the Markov Chain Monte Carlo sampler that generates random samples from the posterior distribution.

The method `sampleUniform(Random rand)` samples a dendrogram topology, uniformly at random (in place). This is a novel algorithm, developed by leveraging Graph Theory, which shows that there exists a bijective function mapping unrooted binary trees with n+1 leaves to the space of rooted binary trees with n leaves. This method uses data structures like Linked Lists, Queues and Hash tables, as well as Recursion based algorithms.

# References

[1] A. Clauset, C. Moore, and M. E. Newman. [Hierarchical structure and the prediction of missing links in networks](https://arxiv.org/pdf/0811.0484.pdf). Nature, 453(7191):98, 2008.



