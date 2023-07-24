## A "Naive" Implementation of the Wavefront Algorithm for Sequence Alignment with Gap-Affine Scoring

This repository contains some simple code that I wrote to understand the wavefront sequence alignment algorithm (Fast gap-affine pairwise alignment using the wavefront algorithm, Bioinformatics, 37(4), 2021, 456–463).

The code is not heavily optimized, and I relied heavily on the FxHashMap for bookkeeping in the implementation. This may not be the most efficient in terms of computational speed. However, using fxHashMap, the code becomes close to the pseudo-code in the paper, although I used a slightly different convention for the sequence coordinates. There are some opportunities to parallelize the computation for further speed improvements too.

Jason Chin, July 23, 2023
