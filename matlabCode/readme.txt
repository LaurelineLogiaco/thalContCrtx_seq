% This code was originally created by Laureline Logiaco

% Here, the optimization of the cost functions uses CMA-ES (an evolutionary algorithm, http://cma.gforge.inria.fr/)
% because in general our cost functions admit several solutions. This is generally less fast (because the optimization is a gradient-free method)
% but this method searches for lower minima and gives robust results.

% This code is still being commented to facilitate its use and is expected to be fully ready by the end of April.
