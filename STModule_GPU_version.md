# STModule GPU version

STModule GPU version is developed based on torch and GPUmatrix.

We recommend to create a conda environment for STModule:

```
conda create -n STModule python=3.9
conda activate STModule
conda install conda-forge::r-base=4.4.1
```

&nbsp;

## Installation of torch for R

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

Installation of GPUmatrix using source package in R:

```r
Sys.setenv(CUDA="11.7")
library(torch)
install.packages("your_path_to_the_file/GPUmatrix_1.0.2.tar.gz", repos = NULL, type = "source")
```

&nbsp;


## Installation of other required packages

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

For high-resolution data, e.g., Slide-seqV2 and Stereo-seq, we recommend using the GPU version of STModule. 

To use the GPU version, we only need to set the parameter `version` to `'gpu'`. For example, for the Slide-seqV2 MOB data:

```r
source('STModule.r')

# data pre-processing 
count_file = 'data/slide_seq_v2_mob_count_matrix.txt'
loc_file = 'data/slide_seq_v2_mob_locations.txt'
data <- data_preprocessing(count_file, loc_file, high_resolution = TRUE, gene_filtering = 50)

# run STModule
res = run_STModule(data, 10, high_resolution = TRUE, max_iter = 100, version = 'gpu')

# visualize spatial maps
plots <- spatial_map_visualization(res, 'log', point_size = 0.3)
pdf('plots/slide_seqv2_mob_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

To apply the tissue modules identified from ST MOB data to Slide-seqV2 data:

```r
source('STModule.r')

# load results of tissue modules identified from ST MOB data
load('results/STModule_res_st_mob.RData')

# estimate spatial maps of the tissue modules for the Slide-seqV2 data
count_file = 'data/slide_seq_v2_mob_count_matrix.txt'
loc_file = 'data/slide_seq_v2_mob_locations.txt'
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution = TRUE, version = 'gpu')

plots <- spatial_map_visualization(spatial_maps)
pdf('plots/spatial_maps_new_mob_data.pdf')
for(p in plots){
    print(p)
}
dev.off()
```









