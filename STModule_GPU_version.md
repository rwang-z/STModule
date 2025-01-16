# STModule GPU version

STModule GPU version is developed based on torch and GPUmatrix.

We recommend to create a conda environment for STModule:

```
conda create -n STModule python=3.9
conda activate STModule
conda install conda-forge::r-base=4.4.1
```

&nbsp;

## Installation of r-torch

Frist, to install CUDA and cuDNN (torch now requires CUDA 11.7):

```
conda install cuda -c nvidia/label/cuda-11.7.0
python3 -m pip install nvidia-cudnn-cu11
```

Sometimes, we need to export path:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:path_to_your_conda/envs/STModule/lib
```

Then, installing torch using `install.packages()` in R:

```r
install.packages("torch")
```

&nbsp;

## Installation of GPUmatrix

Download GPUmatrix source package from <a href="https://cran.r-project.org/web/packages/GPUmatrix/index.html">GPUmatrix</a>

Then run the following code to install GPUmatrix in R:

```r
Sys.setenv(CUDA="11.7")
library(torch)
install.packages("your_path_to_the_file/GPUmatrix_1.0.2.tar.gz", repos = NULL, type = "source")
```

&nbsp;


## Installation of other dependent packages

```
conda install conda-forge::r-seurat
```

```r
install.packages('ggplot2')
install.packages('ggrepel')
install.packages('viridis')
```

&nbsp;

## Usage of STModule GPU version



```r

```











