# MALAT
matlab codes for malat speckles analysis

These scripts are used to quantify the number and sizes of MALAT and NEAT speckles labeled by FISH and imaged using a DeltaVsion deconvolution microscope. The image dataset should be a z-stack of 3 channels (DAPI, MALAT FISH and NEAT FISH).

malat2.m computes the different features of MALAT and NEAT speckles in the cytoplasmic and nuclear compartment of each cell

MalatCellCycle.m identifies mitotic vs interphasic nuclei

WatershedDAPI.m is used to segment nuclei

The MATLAB bioformat package (bfmatlab) should also be installed. It can be downloaded at:
https://downloads.openmicroscopy.org/bio-formats/5.3.4/

