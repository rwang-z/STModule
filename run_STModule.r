source('spatial_factorization.r')
source('spatial_map_estimation.r')

run_STModule = function(data, num_modules, high_resolution = FALSE,  max_iter = 2000){
    # Run STModule on spatially resolved transcriptomics data

    ###### parameters ######
    # data: result of the function data_preprocessing, generated data after data preprocessing
    # num_modules: number of tissue modules to identify
    # high_resolution: to indicate spatial resolution of the data
    #   - TRUE: used for data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc.
    #   - FALSE: used for data profiled by 'ST', '10x Visium', etc. (default)
    # max_iter: maximum iteration, default 2000

    exp_mat = as.matrix(data$exp_mat)
    dist_mat = data$dist_mat
    num_locs = nrow(exp_mat)
    num_genes = ncol(exp_mat)
    print(paste0('Number of locations: ', num_locs))
    print(paste0('Number of genes: ', num_genes))
    print(paste0('Number of factors: ', num_modules))

    # check data
    print('Checking data validity...')
    if(any(is.na(exp_mat))){
    stop('ERROR: the data contains missing data')
    }
    if(any(is.infinite(as.matrix(exp_mat)))){
    stop('ERROR: the data contains infinite data')
    }
    if(!is.numeric(exp_mat)){
    stop('ERROR: the data contains non-numeric data')
    }

    # run factorization
    print('Running STModule...')
    fixed_r = ifelse(num_locs > 1000, TRUE, FALSE)
    r_const = ifelse(high_resolution, 0.01, 1)
    params = list(N = num_locs, L = num_genes, C = num_modules, r_const = r_const, fixed_r = fixed_r, decrease_stop = FALSE,
                a = 1e-6, b = 1e6, c = 1e-6, d = 1e6, e = 1e-6, f = 1e6, g = 0, h = 0, u = 1e-6, v = 1e6, r = 1, z = 1)
    res <- spatial_factorization(params, exp_mat, dist_mat, max_iter)
    print("Estimation finished!")
    print(paste(res$maximumiteration,' Iterations were carried out.'))
    res$params = params
    res$gene_list = colnames(data$exp_mat)
    res$loc_list = rownames(data$exp_mat)
    res$locations = data$loc_info

    return(res)
}

run_spatial_map_estimation = function(res, count_file, loc_file, high_resolution = FALSE,  file_sep = '\t', cell_thresh = 100){
    # estimate spatial maps of tissue module identified from tissue A for another tissue section B

    ###### parameters ######
    # res: the result of 'run_STModule' from tissue section A
    # count_file: path and file name of the count matrix of tissue section B
    # loc_file: path and file name of the spatial information of tissue section B
    # high_resolution: to indicate spatial resolution of the data
    #   - TRUE: used for data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc.
    #   - FALSE: used for data profiled by 'ST', '10x Visium', etc. (default)
    # file_sep: the field separator character of the files for tissue section B, default '\t'
    # cell_thresh: parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts (default 100).
    #   - Used for tissue section B when high_resolution is TRUE

    # genes in tissue modules
    module_genes = res$gene_list
    num_modules = res$params$C
    
    # process new data, using genes included in tissue modules
    print('Loading new data...')
    data = data_preprocessing(count_file, loc_file, high_resolution, gene_mode = 'selected', gene_list = module_genes, 
                                file_sep = file_sep, gene_filtering = 0, cell_thresh = cell_thresh)
    exp_mat = as.matrix(data$exp_mat)
    dist_mat = data$dist_mat
    data_gene_list = colnames(exp_mat)
    num_locs = nrow(exp_mat)
    num_genes = ncol(exp_mat)
    print(paste0('Number of locations: ', num_locs))
    print(paste0('Number of genes: ', num_genes))

    # check new data
    print('Checking data validity...')
    if(any(is.na(exp_mat))){
    stop('ERROR: the data contains missing data')
    }
    if(any(is.infinite(as.matrix(exp_mat)))){
    stop('ERROR: the data contains infinite data')
    }
    if(!is.numeric(exp_mat)){
    stop('ERROR: the data contains non-numeric data')
    }

    # filter tissue module genes with overlapping genes
    module_loadings = as.matrix(res$X$mom1)
    colnames(module_loadings) = res$gene_list
    module_loadings = module_loadings[, data_gene_list]
    module_loadings = as.matrix(module_loadings)

    # estimate spatial maps
    print('Estimating spatial maps of new data ...')
    fixed_r = ifelse(num_locs > 1000, TRUE, FALSE)
    r_const = ifelse(high_resolution, 0.01, 1)
    params = list(N = num_locs, L = num_genes, C = num_modules, r_const = r_const, fixed_r = fixed_r, decrease_stop = FALSE,
                a = 1e-6, b = 1e6, c = 1e-6, d = 1e6, e = 1e-6, f = 1e6, g = 0, h = 0, u = 1e-6, v = 1e6, r = 1, z = 1)
    spatial_map_res <- spatial_map_estimation(params, exp_mat, dist_mat, module_loadings)
    print("Estimation finished!")
    print(paste(spatial_map_res$maximumiteration,' Iterations were carried out.'))
    spatial_map_res$params = params
    spatial_map_res$gene_list = colnames(data$exp_mat)
    spatial_map_res$loc_list = rownames(data$exp_mat)
    spatial_map_res$locations = data$loc_info

    return(spatial_map_res)
}


