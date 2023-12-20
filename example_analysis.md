# Example Analysis

&nbsp;

## Outline
[Souce STModule](#source-STModule)

[Data requirement](#data-requirement)

[Data pre-processing](#data-pre-processing)

[Run STModule](#run-STModule)

[Visualize spatial maps of the tissue modules](#visualize-spatial-maps)

[Associated genes of the tissue modules](#associated_genes)

[Visualize spatial expression interested genes](#visualize-spatial_expression)

[Applying the tissue modules to another tissue section](#apply-to-others)

[Example analysis of human dorsolateral prefrontal cortex profiled by 10x Visium](#example-dlpfc)

[Example analysis of mouse olfactory bulb profiled by Slide-seqV2](#example-mob)

[Example analysis of applying the tissue modules identified from ST MOB data to your own data](#example-mob-application)

&nbsp;

## Source STModule<a id='source-STModule'></a>

```r
> source('STModule.r')
```

&nbsp;

## Data requirement<a id='data-requirement'></a>

STModule requires the following input data of spatially resolved transcriptomics:

**1. count matrix: profiled gene expression at different sptos/cells**

- rows: id of spots/cells

- cols: gene names

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


**2. spatial information: spatial coordinates of the spots/cells**

- rows: id of spots/cells (same ids with the count matrix)

- cols: coordinates of the spots/cells, including two columns ‘x’ and ‘y’

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

## Data pre-processing<a id='data-pre-processing'></a>

Pre-processing and generating input data for STModule using the function **data_preprocessing**:

```r
data_preprocessing(count_file, loc_file, high_resolution = FALSE, gene_mode = 'HVG', top_hvg = 2000, gene_list = c(), 
                   file_sep = '\t', gene_filtering = 0.1, cell_thresh = 100)
```


**Parameters**:
- count_file: path and file name of the count matrix

- loc_file: path and file name of the spatial information

- high_resolution: to indicate spatial resolution of the data
  - TRUE: data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc.
  - FALSE: data profiled by 'ST', '10x Visium', etc. (default)
  
- gene_mode: how to select genes for the analysis
  - 'HVG': to use top highly variable genes selected by Seurat (default)
  - 'selected': to use user-selected genes included in the filtered data, provided by the 'gene_list' parameter
  - 'combined': to use the combination of top HVGs and user-selected genes

- top_hvg: the number of top HVGs to use in the analysis, default 2000

- gene_list: a vector of genes that the user want to use in the analysis, default: c()

- file_sep: the field separator character, default '\t'
  
- gene_filtering: parameter to filter the genes, removing genes expressed in less than gene_filtering locations/cells
  - when high_resolution is FALSE: gene_filtering indicates the threshold of percentage of locations (default 0.1)
  - when high_resolution is TRUE: gene_filtering indicates the threshold of number of cells
  
- cell_thresh: parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts (default 100)
  - used when high_resolution is TRUE

&nbsp;

Example for breast cancer layer 2 profiled by Spatial Transcriptomics:

```r
count_file = 'data/st_bc2_count_matrix.txt'
loc_file = 'data/st_bc2_locations.txt'
data <- data_preprocessing(count_file, loc_file)
```

&nbsp;

## Run STModule<a id='run-STModule'></a>

Run STModule on the generated data using the function **run_STModule**:

```r
run_STModule(data, num_modules, high_resolution = FALSE, max_iter = 2000)
```

**Parameters**:
- data: result of the function ***data_preprocessing***, generated data after data preprocessing
  
- num_modules: number of tissue modules to identify
  
- high_resolution: to indicate spatial resolution of the data
  - TRUE: data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc.
  - FALSE: data profiled by 'ST', '10x Visium', etc. (default)
  
- max_iter: maximum iteration, default 2000

&nbsp;

Example for the breast cancer data:

```r
res <- run_STModule(data, 10)
```

&nbsp;

## Visualize spatial maps of the tissue modules<a id='visualize-spatial-maps'></a>

Plotting spatial maps of the tissue modules using the function **spatial_map_visualization**:

```r
spatial_map_visualization(res, normalization = 'log', point_size = 6)
```

**Parameters**:
- res: result of function **run_STModule**
  
- normalization: normalization method of spatial maps
  - 'log': log-transformation by keeping the direction of activities (default)
  - 'scale': scale the spatial map of each module to a vector with zero-mean and unit-variance
    
- point_size: size of points in the plots
  recommended point size for some SRT technologies:
  - ST data: point_size = 6
  - 10x Visium data: point_size = 2
  - Slide-seqV2 and Stereo-seq data: point_size = 0.3

&nbsp;

Example for the breast cancer data:

```r
plots <- spatial_map_visualization(res, 'log', point_size = 6)
pdf('plots/st_bc2_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Spatial map of the tissue module indicating ductal carcinoma in situ (DCIS):

<img width="300" alt="example_spatial_map" src="https://github.com/rwang-z/STModule/assets/57746198/55977a89-7029-44bd-a27e-1076adf05540">


&nbsp;


## Associated genes of the tissue modules<a id='associated_genes'></a>

### Get associated genes and their activities of the tissue modules

Use the function **get_assocaited_genes** which returns a data frame with three columns ‘gene’, ‘activity’ and ‘module’:

```r
get_assocaited_genes(res)
```

- res: the result of 'run_STModule' running on a tissue section


&nbsp;

Example for the breast cancer data:

```r
module_genes <- get_assocaited_genes(res)
write.table(module_genes, file = 'results/st_bc2_module_associated_genes.txt', sep = '\t', row.names = FALSE)
```

&nbsp;

### Visualize the activities of associated genes:

Plot the activities of the associated genes of the tissue modules using the function **associated_gene_visualization** which returns a list of plots each representing the spatial map of a tissue module:

```r
associated_gene_visualization(res, quantile_thresh = 0.75)
```

**Parameters**:

- res: the result of 'run_STModule' running on a tissue section
  
- quantile_thresh: names will be demonstrated for genes with activity higher than the quantile threshold


&nbsp;

Example for the breast cancer data:

```r
plots <- associated_gene_visualization(res)
pdf('plots/st_pdac_a_st1_associated_gene_activity.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Gene activities of the tissue module indicating DCIS illustrated above:

<img width="500" alt="example_gene_activity" src="https://github.com/rwang-z/STModule/assets/57746198/fad94f9e-8b74-4f03-b43d-d14efb0c2dd8">


&nbsp;

## Visualize spatial expression interested genes<a id='visualize-spatial_expression'></a>

We also provide a function spatial_expression_visualization to plot the spatial expression of a list of interested genes:

```r
spatial_expression_visualization(count_file, loc_file, gene_list, file_sep = '\t', point_size = 6)
```

**Parameters**:

- count_file: path and file name of the count matrix

- loc_file: path and file name of the spatial information

- gene_list: genes to visualize

- file_sep: the field separator character, default '\t'

- point_size: size of points in the plots
  recommended point size for some SRT technologies:
  - ST data: point_size = 6
  - 10x Visium data: point_size = 2
  - Slide-seqV2 and Stereo-seq data: point_size = 0.3

The function plots the spatial expression of the genes into files in the folder ‘plots’.

&nbsp;

Example for the breast cancer data:

```r
count_file = 'data/st_bc2_count_matrix.txt'
loc_file = 'data/st_bc2_locations.txt'
spatial_expression_visualization(count_file, loc_file, c('SPINT2', 'CD74'))
```

<img width="300" alt="example_spatial_expression" src="https://github.com/rwang-z/STModule/assets/57746198/464e3945-11dc-441f-9ef4-6599b43dfc3b">
<img width="300" alt="example_spatial_expression_cd74" src="https://github.com/rwang-z/STModule/assets/57746198/83c0922c-62fc-45a9-817d-99271b1f5edc">

&nbsp;

## Applying the tissue modules to another tissue section<a id='apply-to-others'></a>

To apply the tissue modules identified from a tissue A to another section B, the same input data of section B is required, including the count matrix and spatial information in the same format as mentioned above. 

Spatial maps of the tissue modules on section B can be estimated using the function **run_spatial_map_estimation**:

```r
run_spatial_map_estimation(res, count_file, loc_file, high_resolution = FALSE, file_sep = '\t', cell_thresh = 100)
```

**Parameters**:

- res: the result of ***run_STModule*** from tissue section A

- count_file: path and file name of the count matrix of tissue section B
  
- loc_file: path and file name of the spatial information of tissue section B
  
- high_resolution: to indicate spatial resolution of the SRT data of tissue section B
      - TRUE: data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc.
      - FALSE: data profiled by 'ST', '10x Visium', etc. (default)

- file_sep: the field separator character of the files for tissue section B, default '\t'
  
- cell_thresh: parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts (default 100). 
      - Used for tissue section B when high_resolution is TRUE

The function will first load and pre-process the data of section B and estimate the spatial maps of the tissue modules. The spatial maps can be visualized using the function **spatial_map_visualization**.

&nbsp;

Example for the breast cancer data, estimating the spatial maps of the tissue modules on another section called ‘layer 1’:

```r
count_file = 'data/st_bc1_count_matrix.txt'
loc_file = 'data/st_bc1_locations.txt'
spatial_maps <- run_spatial_map_estimation(res, count_file , loc_file)

plots <- spatial_map_visualization(spatial_maps)
pdf('plots/spatial_maps_st_bc1.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Spatial map of the tissue module indicating DCIS for layer 1:

<img width="300" alt="example_generalization_bc1" src="https://github.com/rwang-z/STModule/assets/57746198/77e71cf3-e1a7-4700-8f6b-57bed5c5ab01">

&nbsp;

## Example analysis of human dorsolateral prefrontal cortex profiled by 10x Visium<a id='example-dlpfc'></a>

```r
source('STModule.r')

# data pre-processing 
count_file = 'data/visium_dlpfc_151676_count_matrix.txt'
loc_file = 'data/visium_dlpfc_151676_locations.txt'
data <- data_preprocessing(count_file, loc_file)

# run STModule
res = run_STModule(data, 10)

# visualize spatial maps
plots <- spatial_map_visualization(res, 'log', point_size = 2)
pdf('plots/visium_dlpfc_151676_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()

# get and save assocaited genes
module_genes <- get_assocaited_genes(res)
write.table(module_genes, file = 'results/visium_dlpfc_151676_module_associated_genes.txt', sep = '\t', row.names = FALSE)

# visualize activities of assocaited genes of the modules
plots <- associated_gene_visualization(res)
pdf('plots/visium_dlpfc_151676_associated_gene_activity.pdf')
for(p in plots){
    print(p)
}
dev.off()

# plot spatial expression of genes in the raw data
count_file = 'data/visium_dlpfc_151676_count_matrix.txt'
loc_file = 'data/visium_dlpfc_151676_locations.txt'
spatial_expression_visualization(count_file, loc_file, c('AQP4', 'PCP4'), point_size = 2)
```

Example of a spatial map:

<img width="300" alt="example_dlpfc" src="https://github.com/rwang-z/STModule/assets/57746198/643ce62b-43cc-4ca2-9886-345e0d08ce5b">


&nbsp;


## Example analysis of mouse olfactory bulb profiled by Slide-seqV2<a id='example-mob'></a>

```r
source('STModule.r')

# data pre-processing 
count_file = 'data/slide_seq_v2_mob_filtered_count_matrix.txt'
loc_file = 'data/slide_seq_v2_mob_filtered_locations.txt'
data <- data_preprocessing(count_file, loc_file, high_resolution = TRUE, gene_filtering = 50)

# run STModule
res = run_STModule(data, 10, high_resolution = TRUE, max_iter = 100)

# visualize spatial maps
plots <- spatial_map_visualization(res, 'log', point_size = 0.3)
pdf('plots/slide_seqv2_mob_spatial_maps.pdf')
for(p in plots){
    print(p)
}
dev.off()

```

Example of a spatial map:

<img width="500" alt="example_mob_slide_seqv2" src="https://github.com/rwang-z/STModule/assets/57746198/be1f58d9-86ae-4bcd-84ee-cd69adbef8a9">


&nbsp;

## Example analysis of applying the tissue modules identified from ST MOB data to your own data<a id='example-mob-application'></a>

We provide the results of STModule on the ST MOB data mentioned the study:

```r
load('results/STModule_res_st_mob.RData')
```

You can apply the tissue modules to other MOB datasets (same data format required as mentioned above) and estimate the corresponding spatial maps using the following code:

```r
# estimate spatial maps on the new data
count_file = 'path to the count matrix of your MOB data'
loc_file = 'path to the spatial information of your MOB data'
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution, file_sep, cell_thresh)

# plot the spatial maps
plots <- spatial_map_visualization(spatial_maps)
pdf('plots/spatial_maps_new_mob_data.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

&nbsp;

Example of applying the tissue modules to the Slide-seqV2 MOB data:

```r
source('STModule.r')

load('results/STModule_res_st_mob.RData')

count_file = 'data/slide_seq_v2_mob_filtered_count_matrix.txt'
loc_file = 'data/slide_seq_v2_mob_filtered_locations.txt'
spatial_maps <- run_spatial_map_estimation(res, count_file, loc_file, high_resolution = TRUE)

plots <- spatial_map_visualization(spatial_maps)
pdf('plots/spatial_maps_new_mob_data.pdf')
for(p in plots){
    print(p)
}
dev.off()
```

Spatial maps of tissue modules identified from ST MOB data:



Estimated spatial maps of the tissue modules for the Slide-seqV2 MOB data:












