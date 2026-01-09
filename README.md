# thalContCrtx_seq
code for 'Thalamic Control of Cortical Dynamics\\in a Model of Flexible Motor Sequencing'

Logiaco, Abbott, Escola

We provide:

1. an extensively commented matlab code to understand and reproduce our main results related to using a fully optimized thalamocortical loop to robustly produce a single motif when starting from the right initial conditions Note that we used this code with optimizers that leverage evolutionary algorithms to facilitate the discovery of nonlocal minima.

2. an optimized python code and associated notebooks to show how eigenvalue control for motif production depends on various parameters, and to quantify how well the optimization of the preparatory network works. Here the optimization of our cost functions uses gradient-based methods to quickly discover local minima. Our results show that gradient-based optimization can reach similar results as evolutionary optimizers as long as the network size is large enough. We notably illustrate these techniques in code that optimize a "preparatory" thalamic module that allows transitions between motifs by setting the right initial state for each motif.

Do not hesitate to contact me (laureline.logiaco@gmail.com) if you need clarifications.
