# GL-AR MATLAB Codebase
This codebase's organization is as follows:

- `main` holds 2 scripts utilizing the various methods in `GL_classes` and is where most readers interested in GL-AR should start
    - `main_synth.m` is an example script of creating a synthetic time-series graph signal using the Graph Tensor Method and then recovering its groundtruth graph structure with GL-AR
    - `main_real.m` is an example script of importing real-world data and using GL-AR to solve for its groundtruth graph structure
- `GL_classes` is where the classes used for computing GL-AR and its subsidiaries are defined
    - `GL.m` defines the graph_learning class where the GL-AR algorithm resides
    - `graphs.m` defines numerous functions related to creating and evaluating graph structures
    - `signals.m` defines signal creation functions, specifically synthetic time-series graph signal creation using the Graph Tensor Method
- `Formal_XXX` folders were used for the official evaluation of GL-AR used for research
- `deprecated_again` is just here for legacy purposes and should provide nothing of interest for the average viewer
