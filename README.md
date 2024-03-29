# ISM
Numerical simulation of the Inertial Spin Model of active particles

Written by Davide Venturelli & Emanuele Loffredo

Last update: 25 April 2022

Original version: 14 March 2019 (Davide Venturelli, Rome),
used in my Master Thesis,

https://zenodo.org/record/7615248

and more recently in this article called "Collective response to local perturbations":

https://doi.org/10.1088/1478-3975/acc5cc

Physical motivation and description of the model can be found, e.g., here:

https://doi.org/10.1007/s10955-014-1119-3

https://doi.org/10.1016/j.physrep.2017.11.003

The program simulates a 2-dimensional system of polar active particles ("flock") with spin alignment interactions, and subject to noise. It can be used in particular to study collective turns induced by an external perturbing field.

It can be adjusted to account for either metric or topological interactions (via a combination of Verlet list + cell list methods). Different boundary conditions (periodic or open) can be chosen.

It runs on a single core, and uses Ciccotti/Vanden-Eijnden method for the integration of the coupled stochastic differential equations, as in

https://doi.org/10.1016/j.cplett.2006.07.086

For further info, feel free to contact me!
