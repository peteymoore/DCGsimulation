# DCGsimulation
Numerical models that simulate coupled ablation and debris transport on debris-covered glaciers (DCG). The model treats ablation as a 1D problem in the vertical dimension, and debris transport as a 1D problem in a single horizontal dimension. Initial conditions for domain dimensions, ice surface profile, initial debris thickness, and englacial debris concentration must be supplied. Parameters for condictive heat transfer and graviational debris movement must also be supplied. The model then runs the following steps in an explicit, finite-difference framework:
* compute (from initial condition or previous iteration) debris thickness
* compute ablation beneath debris
* adjust ice and debris surface
* compute slopes and debris diffusivities
* solve for increment of debris movement across surface
* update debris thickness

The model continues to run either until the specified model time is complete or until all ice is removed from the domain.

# Running the model
A description of the model, the types of simulations available, and example results are provided in Moore (2021). The code is written in Matlab as script files, and each runs "as is" with default parameters. By default, a real-(model)-time visualization of DCG surface profile and debris thickness profile is displayed as the simulation runs. Code for various post-hoc visualizations is also provided in commented sections below the model code. Changing model parameters and solution methods can yield interesting behavior and patterns, as discussed in Moore (2021).

## Reference
Moore, P.L., 2021. Numerical simulation of supraglacial debris mobility: implications for ablation and landform genesis. Frontiers in Earth Science, doi: 10.3389/feart.2021.710131. [https://www.frontiersin.org/articles/10.3389/feart.2021.710131/abstract]
