# STModule R package installation

STModule is developed based on torch.

We recommend to create a conda environment for STModule:

```
conda create -n STModule python=3.9
conda activate STModule
conda install conda-forge::r-base=4.4.1
```

&nbsp;

## Installation of torch for R

**1. installing CUDA and cuDNN (torch now requires CUDA 11.7):**

```
conda install cuda -c nvidia/label/cuda-11.7.0
python3 -m pip install nvidia-cudnn-cu11
```

Export the path if necessary:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:your_path_to_conda/envs/STModule/lib
```

**2. installing torch in R:**

```r
install.packages("torch")
```

Additional installation of torch in R:

```r
Sys.setenv(CUDA="11.7")
library(torch)
```

Please input "yes" to install additional requirements for torch.

&nbsp;


## Installation of Seurat and devtools

We recommend to install Seurat and devtools using conda:

```
conda install conda-forge::r-seurat
conda install conda-forge::r-devtools
```

## Installation of STModule

```r
library(devtools)
devtools::install_github('rwang-z/STModule_package')
```
