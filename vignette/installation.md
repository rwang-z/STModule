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
[Instructions for torch installation](https://torch.mlverse.org/docs/articles/installation.html)

You can install torch using the following method or [install from pre-built binaries](https://torch.mlverse.org/docs/articles/installation.html#pre-built)

**1. Install CUDA and cuDNN (torch now requires CUDA 11.7):**

```
conda install cuda -c nvidia/label/cuda-11.7.0
python3 -m pip install nvidia-cudnn-cu11
```

Export the path if necessary:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:your_path_to_conda/envs/STModule/lib
```

**2. Install torch in R:**

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
devtools::install_github('rwang-z/STModule')
```

STModule uses GPUmatrix 1.0.2. Please **do not** update GPUmatrix during installation.

Please input your GitHub token using `gitcreds::gitcreds_set()` if a GitHub API rate limit error occurs.
