# Supplementary data for a paper on the contact angle of water on free-standing graphene

This repository contains supplementary data supporting the findings of the paper:

> Revealing Strain Effects on the Graphene-Water Contact Angle Using a Machine Learning Potential ([arXiv:2601.20134](https://arxiv.org/abs/2601.20134))

by Darren Wayne Lim, Xavier R. Advincula, William C. Witt, Fabian L. Thiemann, and Christoph Schran.


### Contents

- `00-common`: Common files for various software, as well as the MACE model.

- `01-unstrained`: Simulations for water droplets (of various sizes) on unstrained graphene, either free-standing or frozen, using either the MACE model or SPC/E force field.

- `02-strained`: Simulations for water droplets (of a fixed size) on strained, free-standing graphene, using the MACE model.

- `03-others`: Other simulations for various benchmarking purposes or property characterization.

- `04-plots`: Graphical plots for the paper.

Note that the simulations are designed to run on [ARCHER2](https://www.archer2.ac.uk/) using [LAMMPS](https://www.lammps.org/) patched by [symmetrix](https://github.com/wcwitt/symmetrix). Simulation trajectories have been omitted from the repository due to large filesizes; only processed and analyzed data is present within the repository. All data presented in the paper is reproducible from the saved processed data.

