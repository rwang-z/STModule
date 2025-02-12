# STModule tutorial

All data used in this tutorial is available <a href="https://drive.google.com/drive/folders/15jKtTqfeDMtPaXJgDYeg55aiD-HcvtQs?usp=sharing">here</a>.

In this tutorial, we first demonstrate the step-by-step usage of STModule with a breast cancer dataset profiled by ST *(<a href="https://www.science.org/doi/10.1126/science.aaf2403">Ståhl et al., 2016</a>)*, using the following files:

- st_bc2_count_matrix.txt

- st_bc2_locations.txt

- st_bc1_count_matrix.txt

- st_bc1_locations.txt

Then we provide three additional examples of analysis:

- human dorsolateral prefrontal cortex data profiled by 10x Visium *(<a href="https://www.nature.com/articles/s41593-020-00787-0">Maynard et al., 2021</a>)*;

- mouse olfactory bulb data profiled by Slide-seqV2 *(<a href="https://www.nature.com/articles/s41587-020-0739-1">Stickels et al., 2020</a>)*;

- applying the tissue modules identified from ST MOB data *(<a href="https://www.science.org/doi/10.1126/science.aaf2403">Ståhl et al., 2016</a>)* to other data *(e.g., Slide-seqV2 MOB)*.


&nbsp;

## Outline
[Load STModule](#load-STModule)

[Data requirement](#data-requirement)

[Data pre-processing](#data-pre-processing)

[Run STModule](#run-STModule)

[Visualize spatial maps of the tissue modules](#visualize-spatial-maps)

[Associated genes of the tissue modules](#associated_genes)

[Visualize spatial expression of interested genes](#visualize-spatial_expression)

[Applying the tissue modules to another tissue section](#apply-to-others)

[Example analysis: human dorsolateral prefrontal cortex profiled by 10x Visium](#example-dlpfc)

[Example analysis: mouse olfactory bulb profiled by Slide-seqV2](#example-mob)

[Example analysis: applying the tissue modules identified from ST MOB data to your own data](#example-mob-application)

&nbsp;

## Load STModule<a id='load-STModule'></a>

```r
> library(STModule)
```

&nbsp;

## Data requirements<a id='data-requirement'></a>

STModule requires the following input data of spatially resolved transcriptomics:

**1. Count matrix: profiled gene expression at different sptos/cells**

- rownames: id of spots/cells

- colnames: gene names

```r
          GAPDH USP4 MAPKAPK2 CPEB1 LANCL2
loc_1      1    1        1     0      0
loc_2      7    0        1     0      0
loc_3      5    0        0     0      0
loc_4      1    0        0     0      0
loc_5      2    0        1     0      0
loc_6      9    0        0     0      0
loc_7      4    0        0     0      0
loc_8      4    0        1     0      0
loc_9      3    0        0     0      0
loc_10     4    1        0     0      0
```


**2. Spatial information: spatial coordinates of the spots/cells**

- rownames: id of spots/cells (same ids with the count matrix)

- colnames: coordinates of the spots/cells, including two columns `x` and `y`

```r
          x     y
loc_1  17.907 4.967
loc_2  18.965 5.003
loc_3  18.954 5.995
loc_4  17.846 5.993
loc_5  20.016 6.019
loc_6  20.889 6.956
loc_7  20.062 6.974
loc_8  16.980 6.989
loc_9  17.918 6.991
loc_10 18.877 6.984
```

&nbsp;
&nbsp;

## Data pre-processing<a id='data-pre-processing'></a>

Pre-processing and generating input data for STModule using the function `data_preprocessing()`:

```r
data_preprocessing(count_file, loc_file, high_resolution = FALSE, gene_mode = 'HVG', top_hvg = 2000, gene_list = c(), 
                   file_sep = '\t', gene_filtering = 0.1, cell_thresh = 100)
```


**Parameters**
- `count_file`: path and file name of the count matrix

- `loc_file`: path and file name of the spatial information

- `high_resolution`: a boolean value to indicate spatial resolution of the data.
  - `high_resolution = TRUE`: used for data profiled by Slide-seq, Slide-seqV2, Stereo-seq, etc.
  - `high_resolution = FALSE`: used for data profiled by ST, 10x Visium, etc. (default)
  
- `gene_mode`: the way to select genes for the analysis. 
  - `gene_mode = 'HVG'`: to use top highly variable genes selected by Seurat (default)
  - `gene_mode = 'selected'`: to use user-selected genes included in the filtered data, provided by the `gene_list` parameter
  - `gene_mode = 'combined'`: to use the combination of top HVGs and user-selected genes

- `top_hvg`: the number of top HVGs to use in the analysis, default 2000

- `gene_list`: a vector of genes provided by the user to use in the analysis, default: c()

- `file_sep`: the field separator character, default '\t'
  
- `gene_filtering`: parameter to filter the genes, removing genes expressed in less than gene_filtering locations/cells
  - when `high_resolution = FALSE`: gene_filtering indicates the threshold of percentage of locations (default 0.1)
  - when `high_resolution = TRUE`: gene_filtering indicates the threshold of number of cells
  
- `cell_thresh`: parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts (default 100)
  - used when `high_resolution = TRUE`


**Output**

- A list containing the processed expression matrix, distance matrix and location information


**Usage**

- For data profiled by **ST** and **10x Visium**, we recommend using the default setting:

```r
data <- data_preprocessing(count_file, loc_file)
```

- For **high-resolution data** profiled by Slide-seq, Slide-seqV2, Stereo-seq, etc., we recommend using the following setting to exclude genes and cells with too sparse expression:

```r
data <- data_preprocessing(count_file, loc_file, high_resolution = TRUE, gene_filtering = 50, cell_thresh = 100)
```


**Example**

Example for breast cancer layer 2 profiled by ST:

```r
count_file = 'st_bc2_count_matrix.txt'
loc_file = 'st_bc2_locations.txt'
data <- data_preprocessing(count_file, loc_file)
```

&nbsp;
&nbsp;

## Run STModule<a id='run-STModule'></a>

Run STModule on the generated data using the function `run_STModule()`:

```r
run_STModule(data, num_modules, high_resolution = FALSE, max_iter = 2000, version = 'cpu')
```

**Parameters**
- `data`: result of the function `data_preprocessing()`, generated data after data preprocessing
  
- `num_modules`: number of tissue modules to identify.
  
  To identify major expression components of a tissue section, we recommend setting `num_modules = 10`. If you want to identify more detailed components, please use a larger number for this parameter.
  
- `high_resolution`: a boolean value to indicate spatial resolution of the data
  - `high_resolution = TRUE`: used for data profiled by Slide-seq, Slide-seqV2, Stereo-seq, etc.
  - `high_resolution = FALSE`: used for data profiled by ST, 10x Visium, etc. (default)

- `max_iter`: number of maximum iterations, default 2000

  For this parameter, we recommend using the default setting for data profiled by ST and 10x Visium, and setting `max_iter = 100` for high-resolution data for efficiency.

- `version`: to use CPU or GPU version of STModule
  - `version = 'cpu'`: used for ST data (default)
  - `version = 'gpu'`: used for 10x Visium and high-resolution data

**Output**

- A list containing estimated results for variables in the model, as well as parameters, gene list, location list and location information.

**Usage**

- For data profiled by **ST**, we recommend using the default setting:

```r
res <- run_STModule(data, 10)
```

- For data profiled by **10x Visium**, we recommend using the GPU version:

```r
res = run_STModule(data, 10, version = 'gpu')
```

- For **high-resolution data** profiled Slide-seq, Slide-seqV2, Stereo-seq, etc., we recommend using the following setting:

```r
res = run_STModule(data, 10, high_resolution = TRUE, max_iter = 100, version = 'gpu')
```

**Example**

Example for the breast cancer data:

```r
res <- run_STModule(data, 10)
```

&nbsp;
&nbsp;

## Visualize spatial maps of the tissue modules<a id='visualize-spatial-maps'></a>

Plotting spatial maps of the tissue modules using the function `spatial_map_visualization()`:

```r
spatial_map_visualization(res, normalization = 'log', point_size = 6)
```

**Parameters**
- `res`: the result of function `run_STModule()`
  
- `normalization`: normalization method for spatial maps
  - `normalization = 'log'`: log-transformation while keeping the direction of activities (default)
  - `normalization = 'scale'`: scale the spatial map of each module to a vector with zero-mean and unit-variance
    
- `point_size`: size of points in the plots. Recommended point size for data profiled by different SRT technologies:
  - ST data: `point_size = 6`
  - 10x Visium data: `point_size = 2`
  - Slide-seqV2 and Stereo-seq data: `point_size = 0.3`

**Output**

- A list of plots each for a tissue module

**Usage**

- For data profiled by **ST**, we recommend using the following setting:

```r
plots <- spatial_map_visualization(res, 'log', point_size = 6)
```

- For data profiled by **10x Visium**, we recommend using the following setting:

```r
plots <- spatial_map_visualization(res, 'log', point_size = 2)
```

- For **high-resolution data** profiled Slide-seq, Slide-seqV2, Stereo-seq, etc., we recommend using the following setting:

```r
plots <- spatial_map_visualization(res, 'log', point_size = 0.3)
```

**Example**

Example for the breast cancer data:

```r
plots <- spatial_map_visualization(res, 'log', point_size = 6)
pdf('st_bc2_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Spatial map of the tissue module representing ductal carcinoma in situ (DCIS):

<img width="300" alt="example_spatial_map" src="https://github.com/rwang-z/STModule/assets/57746198/55977a89-7029-44bd-a27e-1076adf05540">


&nbsp;
&nbsp;


## Associated genes of the tissue modules<a id='associated_genes'></a>

### Get associated genes and their activities of the tissue modules

Use the function `get_assocaited_genes()` to get the associated genes of the tissue modules:

```r
get_assocaited_genes(res)
```

**Parameter**
- `res`: the result of function `run_STModule()`

**Output**

- A data frame with three columns:

  - `gene`: gene name

  - `module`: the index of the tissue module

  - `activity`: loading of the gene in the tissue module

**Example**

Example for the breast cancer data:

```r
module_genes <- get_assocaited_genes(res)
write.table(module_genes, file = 'st_bc2_module_associated_genes.txt', sep = '\t', row.names = FALSE)
```

&nbsp;

### Visualize the activities of associated genes:

Plot the activities of the associated genes of the tissue modules using the function `associated_gene_visualization()`:

```r
associated_gene_visualization(res, quantile_thresh = 0.75)
```

**Parameters**

- `res`: the result of `run_STModule` running on a tissue section
  
- `quantile_thresh`: names will be demonstrated for genes with activity higher than the quantile threshold


**Output**

- A list of plots each demonstrating the associated genes of a tissue module

**Example**

Example for the breast cancer data:

```r
plots <- associated_gene_visualization(res)
pdf('st_bc2_associated_gene_activity.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Gene activities of the tissue module representing DCIS:

<img width="500" alt="example_gene_activity" src="https://github.com/rwang-z/STModule/assets/57746198/fad94f9e-8b74-4f03-b43d-d14efb0c2dd8">


&nbsp;
&nbsp;

## Visualize spatial expression of interested genes<a id='visualize-spatial_expression'></a>

We also provide a function `spatial_expression_visualization()` to plot the spatial expression of a list of interested genes:

```r
spatial_expression_visualization(count_file, loc_file, gene_list, file_sep = '\t', point_size = 6)
```

**Parameters**

- `count_file`: path and file name of the count matrix

- `loc_file`: path and file name of the spatial information

- `gene_list`: genes to visualize

- `file_sep`: the field separator character, default '\t'

- `point_size`: size of points in the plots. Recommended point size for some SRT technologies:
  - ST data: `point_size = 6`
  - 10x Visium data: `point_size = 2`
  - Slide-seqV2 and Stereo-seq data: `point_size = 0.3`

**Output**

- PDF files, each demonstrating the spatial expression of a query gene

**Usage**

- For data profiled by **ST**, we recommend using the following setting:

```r
spatial_expression_visualization(count_file, loc_file, gene_list, point_size = 6)
```

- For data profiled by **10x Visium**, we recommend using the following setting:

```r
spatial_expression_visualization(count_file, loc_file, , gene_list, point_size = 2)
```

- For **high-resolution data** profiled Slide-seq, Slide-seqV2, Stereo-seq, etc., we recommend using the following setting:

```r
spatial_expression_visualization(count_file, loc_file, , gene_list, point_size = 0.3)
```

**Example**

Example for the breast cancer data:

```r
count_file = 'st_bc2_count_matrix.txt'
loc_file = 'st_bc2_locations.txt'
spatial_expression_visualization(count_file, loc_file, c('SPINT2', 'CD74'), point_size = 6)
```

Spatial expression maps of SPINT2 and CD74 in the breast cancer data:

<img width="300" alt="example_spatial_expression" src="https://github.com/rwang-z/STModule/assets/57746198/464e3945-11dc-441f-9ef4-6599b43dfc3b">
<img width="300" alt="example_spatial_expression_cd74" src="https://github.com/rwang-z/STModule/assets/57746198/83c0922c-62fc-45a9-817d-99271b1f5edc">

&nbsp;
&nbsp;

## Applying the tissue modules to another tissue section<a id='apply-to-others'></a>

To apply the tissue modules identified from a tissue A to another section B, the same input data of section B is required, including the count matrix and spatial information in the same format as mentioned above. 

Spatial maps of the tissue modules for section B can be estimated using the function `run_spatial_map_estimation()`:

```r
run_spatial_map_estimation(res, count_file, loc_file, high_resolution = FALSE, file_sep = '\t', cell_thresh = 100, max_iter = 2000, version = 'cpu')
```

**Parameters**

- `res`: the result of function `run_STModule()` from tissue section A

- `count_file`: path and file name of the count matrix of tissue section B
  
- `loc_file`: path and file name of the spatial information of tissue section B
  
- `high_resolution`: a boolean value to indicate spatial resolution of the data
  - `high_resolution = TRUE`: used for data profiled by Slide-seq, Slide-seqV2, Stereo-seq, etc.
  - `high_resolution = FALSE`: used for data profiled by ST, 10x Visium, etc. (default)

- `file_sep`: the field separator character of the files for tissue section B, default '\t'
  
- `cell_thresh`: parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts (default 100). 
  - Used for tissue section B when `high_resolution = TRUE`

- `max_iter`: number of maximum iteration, default 2000

  For this parameter, we recommend using the default setting for data profiled by ST and 10x Visium, and setting `max_iter = 100` for high-resolution data for efficiency.

- `version`: to use CPU or GPU version of STModule
  - `version = 'cpu'`: used for ST data (default)
  - `version = 'gpu'`: used for 10x Visium and high-resolution data

The function will first load and pre-process the data of section B and then estimate the spatial maps of the tissue modules. The spatial maps can also be visualized using the function `spatial_map_visualization()`.

**Output**

- A list containing estimated results for variables in the model, as well as parameters, gene list, location list and location information

**Usage**

- For data profiled by **ST**, we recommend using the default setting:

```r
spatial_maps <- run_spatial_map_estimation(res, count_file , loc_file)
```

- For data profiled by **10x Visium**, we recommend using the GPU version:

```r
spatial_maps <- run_spatial_map_estimation(res, count_file , loc_file, version = 'gpu')
```

- For **high-resolution data** profiled Slide-seq, Slide-seqV2, Stereo-seq, etc., we recommend using the following setting:

```r
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution = TRUE, max_iter = 100, version = 'gpu')
```

**Example**

Example for the breast cancer data, estimating the spatial maps of the tissue modules on another section ‘layer 1’:

```r
# estimate the spatial maps
count_file = 'st_bc1_count_matrix.txt'
loc_file = 'st_bc1_locations.txt'
spatial_maps <- run_spatial_map_estimation(res, count_file , loc_file)

# plot the spatial maps
plots <- spatial_map_visualization(spatial_maps)
pdf('spatial_maps_st_bc1.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Spatial map of the tissue module representing DCIS for layer 1:

<img width="300" alt="example_generalization_bc1" src="https://github.com/rwang-z/STModule/assets/57746198/77e71cf3-e1a7-4700-8f6b-57bed5c5ab01">

&nbsp;
&nbsp;

## Example analysis: human dorsolateral prefrontal cortex profiled by 10x Visium<a id='example-dlpfc'></a>

```r
library('STModule')

# data pre-processing 
count_file = 'visium_dlpfc_151676_count_matrix.txt'
loc_file = 'visium_dlpfc_151676_locations.txt'
data <- data_preprocessing(count_file, loc_file)

# run STModule, using the GPU version
res = run_STModule(data, 10, version = 'gpu')

# visualize spatial maps
plots <- spatial_map_visualization(res, 'log', point_size = 2)
pdf('visium_dlpfc_151676_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()

# get and save assocaited genes
module_genes <- get_assocaited_genes(res)
write.table(module_genes, file = 'visium_dlpfc_151676_module_associated_genes.txt', sep = '\t', row.names = FALSE)

# visualize activities of assocaited genes of the modules
plots <- associated_gene_visualization(res)
pdf('visium_dlpfc_151676_associated_gene_activity.pdf')
for(p in plots){
    print(p)
}
dev.off()

# plot spatial expression of genes in the raw data
count_file = 'visium_dlpfc_151676_count_matrix.txt'
loc_file = 'visium_dlpfc_151676_locations.txt'
spatial_expression_visualization(count_file, loc_file, c('AQP4', 'PCP4'), point_size = 2)
```

Demonstration of a spatial map:

<img width="300" alt="example_dlpfc" src="https://github.com/rwang-z/STModule/assets/57746198/643ce62b-43cc-4ca2-9886-345e0d08ce5b">


&nbsp;
&nbsp;

## Example analysis: mouse olfactory bulb profiled by Slide-seqV2<a id='example-mob'></a>

```r
library('STModule')

# data pre-processing 
count_file = 'slide_seq_v2_mob_count_matrix.txt'
loc_file = 'slide_seq_v2_mob_locations.txt'
data <- data_preprocessing(count_file, loc_file, high_resolution = TRUE, gene_filtering = 50, cell_thresh = 100)

# run STModule, using the GPU version
res = run_STModule(data, 10, high_resolution = TRUE, max_iter = 100, version = 'gpu')

# visualize spatial maps
plots <- spatial_map_visualization(res, 'log', point_size = 0.3)
pdf('slide_seqv2_mob_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()

# get and save assocaited genes
module_genes <- get_assocaited_genes(res)
write.table(module_genes, file = 'slide_seqv2_mob_module_associated_genes.txt', sep = '\t', row.names = FALSE)

# visualize activities of assocaited genes of the modules
plots <- associated_gene_visualization(res)
pdf('slide_seqv2_mob_associated_gene_activity.pdf')
for(p in plots){
    print(p)
}
dev.off()

```

Demonstration of a spatial map:

<img width="500" alt="example_mob_slide_seqv2" src="https://github.com/rwang-z/STModule/assets/57746198/be1f58d9-86ae-4bcd-84ee-cd69adbef8a9">

&nbsp;
&nbsp;

## Example analysis: applying the tissue modules identified from ST MOB data to your own data<a id='example-mob-application'></a>

We provide the results of STModule for the ST MOB data *(<a href="https://www.science.org/doi/10.1126/science.aaf2403">Ståhl et al., 2016</a>)*, available <a href="https://drive.google.com/drive/folders/15jKtTqfeDMtPaXJgDYeg55aiD-HcvtQs">here</a>. 

You can load the result (STModule_res_st_mob.RData), apply the tissue modules to other MOB datasets (same data format required as mentioned above), and estimate the corresponding spatial maps as follows:

```r
library('STModule')
load('STModule_res_st_mob.RData')

# estimate spatial maps on the new data
count_file = 'path to the count matrix of your MOB data'
loc_file = 'path to the spatial information of your MOB data'
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution, file_sep, cell_thresh, max_iter, version)

# plot the spatial maps
plots <- spatial_map_visualization(spatial_maps, normalization, point_size)
pdf('spatial_maps_new_mob_data.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

&nbsp;

Example of applying the tissue modules to the Slide-seqV2 MOB data:

```r
library('STModule')
load('STModule_res_st_mob.RData')

count_file = 'slide_seq_v2_mob_count_matrix.txt'
loc_file = 'slide_seq_v2_mob_locations.txt'
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution = TRUE, max_iter = 100, version = 'gpu')

plots <- spatial_map_visualization(spatial_maps, 'log', point_size = 0.3)
pdf('spatial_maps_new_mob_data.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

&nbsp;

Demonstration of several tissue modules identified from ST MOB data and corresponding spatial maps for the Slide-seqV2 MOB data:

<img width="300" alt="mob_st_1" src="https://github.com/rwang-z/STModule/assets/57746198/fd84f490-08bb-4fd8-b3b0-ea77ed08e8d2">
<img width="300" alt="mob_st_2" src="https://github.com/rwang-z/STModule/assets/57746198/be0c502d-c77a-4269-88cf-74286005a379">
<img width="300" alt="mob_st_3" src="https://github.com/rwang-z/STModule/assets/57746198/b0ac93b9-99f7-4e8c-9d37-d27fe140f823">

<img width="300" alt="mob_st_to_slide_seqcv2_1" src="https://github.com/rwang-z/STModule/assets/57746198/db8d0023-890c-4321-94cd-cc56df48ca16">
<img width="300" alt="mob_st_to_slide_seqcv2_2" src="https://github.com/rwang-z/STModule/assets/57746198/e78d597b-b518-4016-97d3-bbb0f338eae1">
<img width="300" alt="mob_st_to_slide_seqcv2_3" src="https://github.com/rwang-z/STModule/assets/57746198/3676ddd8-9146-4d5c-bafe-5ceeabf5730b">
