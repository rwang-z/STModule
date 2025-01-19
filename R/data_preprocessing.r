#' Data preprocessing, to generate the input data for STModule from the user-provided data.
#' input data: count matrix & spatial information
#' - count matrix: locations * genes. rownames - id of spots/cells; colnames - gene names.
#' - spatial information: locations * coordinates. rownames - id of spots/cells; colnames - x, y.
#' 
#' @param count_file path and file name of the count matrix.
#' @param loc_file path and file name of the spatial information.
#' @param high_resolution a boolean value to indicate spatial resolution of the data. 
#' - TRUE for data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc. 
#' - FALSE for data profiled by 'ST' and '10x Visium'.
#' @param gene_mode the way to select genes for the analysis. 
#' - 'HVG': to use top highly variable genes selected by Seurat. 
#' - 'selected': to use user-selected genes included in the filtered data, provided by the 'gene_list' parameter.
#' - 'combined': to use the combination of top HVGs and user-selected genes.
#' @param top_hvg the number of top HVGs to use in the analysis, default 2000.
#' @param gene_list a vector of genes that the user want to use in the analysis.
#' @param file_sep the field separator character.
#' @param gene_filtering parameter to filter the genes, removing genes expressed in less than gene_filtering locations/cells.
#' - when 'high_resolution' is FALSE: gene_filtering indicates the threshold of percentage of locations.
#' - when 'high_resolution' is TRUE: gene_filtering indicates the threshold of number of cells.
#' @param cell_thresh parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts.
#' @import Seurat
#' @return a list containing the processed expression matrix, distance matrix and location information.
#' @export data_preprocessing

data_preprocessing = function(count_file, loc_file, high_resolution = FALSE, gene_mode = 'HVG', top_hvg = 2000, gene_list = c(), 
                                file_sep = '\t', gene_filtering = 0.1, cell_thresh = 100){
    # load count matrix
    print('Performing data pre-processing')
    print('Loading count matrix...')
    count_mat = read.table(count_file, header = TRUE, sep = file_sep, row.names = 1)
    print(paste0('Raw count mat: ', nrow(count_mat), ' spots, ', ncol(count_mat), ' genes'))

    # spatial information
    print('Loading spatial information')
    locations = read.table(loc_file, sep = file_sep, header = TRUE, row.names = 1)

    # data filtering, genes and locations/cells
    print('Data filtering...')
    filtered_count = data_filtering(count_mat, high_resolution, gene_filtering, cell_thresh)

    # normalization and HVG selection
    print('Normalizing data...')
    dataobj <- CreateSeuratObject(counts = t(filtered_count))
    dataobj <- NormalizeData(dataobj, normalization.method = "RC", scale.factor = 10000) # default setting
    dataobj <- FindVariableFeatures(dataobj, selection.method = 'vst', nfeatures = top_hvg)
    hvgs <- VariableFeatures(dataobj)
    # transformed_mat = GetAssayData(object = dataobj, slot = 'data')  # Seurat v4
    transformed_mat = dataobj[["RNA"]]$data     # LayerData(dataobj, assay = "RNA", layer = "data")
    transformed_mat = t(as.matrix(transformed_mat))
    if(gene_mode == 'HVG'){
        selected_genes = hvgs
    }else if(gene_mode == 'selected'){
        selected_genes = intersect(gene_list, colnames(transformed_mat))
    }else if(gene_mode == 'combined'){
        selected_genes = union(hvgs, intersect(gene_list, colnames(transformed_mat)))
    }
    transformed_mat = remove_zero_exp_locations(transformed_mat[,selected_genes])
    filtered_locations = locations[rownames(transformed_mat),]

    # calculate distance matrix of locations
    print('Calculating distance between locations/cells...')
    dist_mat = cal_distance_mat(filtered_locations)

    # result data
    data = list(exp_mat = transformed_mat, 
                dist_mat = dist_mat, 
                loc_info = filtered_locations)
    print('Data prepared!')
    return(data)
}


