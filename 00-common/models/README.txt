The file "ch2o-dens-inv_swa.model" contains the MACE machine learning potential as a PyTorch model. This machine learning potential was trained on revPBE-D3 data, for various configurations of water interacting with graphene and carbon nanotubes amongst other training data.

As a PyTorch model, the "ch2o-dens-inv_swa.model" file can be loaded directly into the MACE Python module, e.g. if using within ASE. It can also be loaded into LAMMPS through the MLIAP implementation (although this is not the intended usage).

However, the file CANNOT be used by symmetrix as-is. It needs to be converted into a JSON format first, in order to load it into LAMMPS+symmetrix. The JSON file was excluded from the repository due to large filesize.

To perform this conversion, install both the MACE and symmetrix Python modules, and then run:

>    symmetrix_extract_mace ch2o-dens-inv_swa.model --atomic-numbers 1 6 8

For more details, see the symmetrix documentation: https://github.com/wcwitt/symmetrix/blob/main/symmetrix/README.md
