Code related to the publication: 

*The CALIPR framework for highly accelerated myelin water imaging with improved precision and sensitivity*. 

# Supporting Data
Supporting data is available at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8049179.svg)](https://doi.org/10.5281/zenodo.8049179)

and includes brain and spinal cord image data for:
- conventional structural images
- CALIPR multi-echo image volumes
- CALIPR quantitative metric maps (MWF, IET2)
as well as the data required to reproduce the CALIPR simulations results (*FS_dat.mat*).


# Outline of the CALIPR Framework

### Acquisition
![calipr_fig8_acq](https://github.com/avdvorak/CALIPR/assets/24612184/7713e00d-a742-4b56-bdf0-9c7c042df5d7)


### Reconstruction
![calipr_fig8_recon](https://github.com/avdvorak/CALIPR/assets/24612184/a8c59b0c-b1f0-48ab-a24f-41d3e36e6c9d)


# Dependencies for Image Reconstruction


### Gyrotools ReconFrame

For Philips scanner data, the Gyrotools ReconFrame environment (https://www.gyrotools.com/gt/index.php/products/reconframe) is used to read the raw k-space data and perform hardware corrections, specifically using the Matlab Reconstruction Library (MRecon, version 5.1.0).


### Berkeley Advanced Reconstruction Toolbox (BART)

Image reconstruction is performed using the Berkeley Advanced Reconstruction Toolbox (BART, version v0.7.00)(https://github.com/mrirecon/bart/releases/tag/v0.7.00).

Information on installing and using BART can be found here:
https://mrirecon.github.io/bart/


# Dependencies for Image Processing and Analysis

The following dependencies are used for additional processing, but not required to reconstruct CALIPR data.


### Advanced Normalization Tools (ANTs)

The pipeline used the Advanced Normalization Tools software (https://github.com/ANTsX/ANTs). 
Instructions for compiling and setting up ANTs can be found here:
https://github.com/ANTsX/ANTs/wiki/Compiling-ANTs-on-Linux-and-Mac-OS


### FSL
The script also uses some generic functions from FSL. Instructions for installing FSL can be found here:
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation


### Spinal Cord Toolbox
The spinal cord processing script also uses some functions from Spinal Cord Toolbox (SCT). Instructions for installing SCT can be found here:
https://spinalcordtoolbox.com/user_section/installation.html


### dcm2niix

The dcm2niix tools is used to convert data from the DICOM or PAR/REC format to the NIfTI format:
https://github.com/rordenlab/dcm2niix


### OASIS Template

The OASIS brain template is used with ANTs to perform brain extraction, and is available at:
https://figshare.com/articles/dataset/ANTs_ANTsR_Brain_Templates/915436


### Matlab Kernel for Jupyter Notebooks

Simulations were performed using Jupyter notebooks running a Matlab (2019b) kernel, with installation instructions available here:
https://am111.readthedocs.io/en/latest/jmatlab_install.html

With Matlab 2020b or later, the newly available Jupyter-Matlab proxy from Mathworks can be used:
https://github.com/mathworks/jupyter-matlab-proxy
