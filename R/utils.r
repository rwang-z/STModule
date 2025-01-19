# utility functions

#' Calculate matrix inverse
#' 
#' @param mat matrix to calculate inverse.
#' @param params parameters of the estimation process. 
#' @return inverse matrix of the input matrix.
#' @export cal_mat_inv

cal_mat_inv = function(mat, params){
    if(params$version == 'cpu'){
        mat_inv = tryCatch(chol2inv(chol(mat)),
                            error = function(e){
                                tryCatch(return(solve(mat)),
                                        error = function(e){
                                            # add noise to cov_mat
                                            mat_noise = mat + diag(runif(dim(mat)[1], -1e-6, 0))
                                            return(solve(mat_noise))
                                        })
                            })
    }else if(params$version == 'gpu'){
        mat_inv = tryCatch(solve(mat),
                            error = function(e){return(ginv(mat))})
    }
    return(mat_inv)
}


#' Calculate inverse of covariance matrix
#' 
#' @param r variable r of the SE kernel.
#' @param dist_mat Euclidean distance matrix of the spots/cells.
#' @param params parameters of the estimation process. 
#' @return inverse matrix of the covariance matrix.
#' @export cal_cov_inv_mat

cal_cov_inv_mat = function(r, dist_mat, params){
    cov_mat = exp(-0.5 * r * dist_mat)
    cov_inv = cal_mat_inv(cov_mat, params)
    return(cov_inv)
}


#' Calculate the inverse of covariance matrices for all tissue modules.
#' 
#' @param r a list of variable r of the SE kernel for all tissue modules.
#' @param dist_mat Euclidean distance matrix of the spots/cells.
#' @param params parameters of the estimation process. 
#' @return a list of inverse matrices of the covariance matrices of all tissue moduels.
#' @export cal_cov_inv_list

cal_cov_inv_list = function(r, dist_mat, params){
    num_comp = dim(r)[1]
    cov_list = list()
    for(c in c(1:num_comp)){
        cov_list[[c]] = cal_cov_inv_mat(r[c], dist_mat, params)
    }
    return(cov_list)
}


#' Calculate the inverse of covariance matrices for all tissue modules, same for all modules.
#' 
#' @param r_const constant r of the SE kernel for all tissue modules.
#' @param num_comp number of tissue modules.
#' @param dist_mat Euclidean distance matrix of the spots/cells.
#' @param params parameters of the estimation process. 
#' @return a list of inverse matrices of the covariance matrices of all tissue moduels.
#' @export cal_cov_inv_list_same

cal_cov_inv_list_same = function(r_const, num_comp, dist_mat, params){
    cov_inv = cal_cov_inv_mat(r_const, dist_mat, params)
    cov_list = list()
    for(c in c(1:num_comp)){
        cov_list[[c]] = cov_inv
    }
    return(cov_list)
}


#' Calculate logarithm of determinant of a precision matrix.
#' 
#' @param prec_mat a precision matrix for spatial maps.
#' @return the logarithm of determinant of a precision matrix.
#' @export cal_precision_log_det

cal_precision_log_det = function(prec_mat){
    det_log = tryCatch(determinant(prec_mat)$modulus[1],
                        error = function(e){
                            return(2 * determinant(chol(prec_mat))$modulus[1])
                        })
    return(det_log)
}


#' Get diagonal of precision matrices.
#' 
#' @param prec_list a list of precision matrices of spatial maps.
#' @param dimension dimension of precision matrix.
#' @param num_comp number of tissue modules.
#' @return diagonal of precision matrices.
#' @export get_prec_mat_diag

get_prec_mat_diag = function(prec_list, dimension, num_comp){
    prec_diag_mat = matrix(0, dimension, num_comp)
    for(c in c(1:num_comp)){
        prec_diag_mat[, c] = diag(prec_list[[c]])
    }
    return(prec_diag_mat)
}


#' Calculate Euclidean distance matrix of spots/cells based on location information.
#' 
#' @param locations location information of spots/cells.
#' @param dimension dimension of the coordinates. '2D' or '3D'.
#' @return a Euclidean distance matrix of spots/cells.
#' @export cal_distance_mat

cal_distance_mat = function(locations, dimension = '2D'){
    num_loc = dim(locations)[1]
    if(dimension == '2D'){
        mat = locations[, c('x', 'y')]
    }else if(dimension == '3D'){
        print('calculating distance for 3D data')
        mat = locations[, c('x', 'y', 'z')]
    }
    dist_mat = as.matrix(dist(mat, method = "euclidean"))
    dist_mat = dist_mat^2
    return(dist_mat)
}


#' Convert results of GPU versions to matrices.
#' 
#' @param res the results of 'run_STModule' or 'run_spatial_map_estimation' running on a tissue section.
#' @param type type of estimation.
#' - 'all': for resutls of 'run_STModule'.
#' - 'estimate': for results of 'run_spatial_map_estimation'.
#' @return processed results.
#' @export convert_res_to_matrix

convert_res_to_matrix = function(res, type = 'all'){
    res$A = list(mu = as.matrix(res$A$mu))
    res$R = list(r = as.matrix(res$R$r))
    res$Delt = list(c = as.matrix(res$Delt$c), d = as.matrix(res$Delt$d))
    res$Lam = list(u = as.matrix(res$Lam$u), v = as.matrix(res$Lam$v))
    if(type == 'all'){
        res$Beta = list(e = as.matrix(res$Beta$e), f = as.matrix(res$Beta$f))
        res$Ph = as.matrix(res$Ph)
        res$Ps = as.matrix(res$Ps)
        res$Rho = as.matrix(res$Rho)
        res$X = list(gamma = as.matrix(res$X$gamma), sigma = as.matrix(res$X$sigma), m = as.matrix(res$X$m))
        res$X$mom1 = as.matrix(res$X$m * res$X$gamma)
    }else if(type == 'estimate'){
        res$X$mom1 = as.matrix(res$X$mom1)
    }
    return(res)
}


#' Filtering input data, removing genes and spots/cells by the criteria specified by the parameters
#' 
#' @param count_mat count matrix of input data
#' @param high_resolution a boolean value to indicate spatial resolution of the data. 
#' - TRUE for data profiled by 'Slide-seq', 'Slide-seqV2', 'Stereo-seq', etc. 
#' - FALSE for data profiled by 'ST' and '10x Visium'.
#' @param gene_filtering parameter to filter the genes, removing genes expressed in less than gene_filtering locations/cells.
#' - when 'high_resolution' is FALSE: gene_filtering indicates the threshold of percentage of locations.
#' - when 'high_resolution' is TRUE: gene_filtering indicates the threshold of number of cells.
#' @param cell_thresh parameter to filter cells for high-resolution data, removing cells with less than cell_thresh counts.
#' @return the count matrix after filtering
#' @export data_filtering

data_filtering = function(count_mat, high_resolution, gene_filtering, cell_thresh = 100){
    if(high_resolution){
        # gene_filtering: number of spots/cells
        if(gene_filtering > 0){
            gene_exp_cell = apply(count_mat, 2, function(x) length(which(x > 0))) # number of cells expressing each gene
            gene_remained = colnames(count_mat)[which(gene_exp_cell >= gene_filtering)]
        }else{
            gene_remained = colnames(count_mat)
        }
        cell_total_counts = rowSums(count_mat)
        cell_remained = rownames(count_mat)[which(cell_total_counts >= cell_thresh)]
        filtered_mat = count_mat[cell_remained, gene_remained]
    }else{
        # gene_filtering: percentage
        if(gene_filtering > 0){
            print(paste0('Removing genes expressed in less than ', gene_filtering * 100, '% locations'))
            gene_loc_exp = apply(count_mat, 2, function(x) length(which(x > 0)))
            remove_genes = which(gene_loc_exp < gene_filtering * nrow(count_mat))
            if(length(remove_genes) > 0){
                filtered_mat = count_mat[, -remove_genes]
            }else{
                filtered_mat = count_mat
            }
        }else{
            filtered_mat = count_mat
        }
        loc_total_counts = rowSums(filtered_mat)
        loc_remained = rownames(filtered_mat)[which(loc_total_counts > 0)]
        filtered_mat = filtered_mat[loc_remained, ]
    }
    print(paste0(ncol(filtered_mat), ' genes remained after filteration'))
    print(paste0(nrow(filtered_mat), ' locations/cells remained after filteration'))
    return(filtered_mat)
}


#' Removing spots/cells with zero expression
#' 
#' @param count_mat count matrix
#' @return the count matrix after filtering
#' @export remove_zero_exp_locations

remove_zero_exp_locations = function(count_mat){
    loc_sum = rowSums(count_mat)
    remove_locs = which(loc_sum == 0)
    if(length(remove_locs) > 0){
        count_mat = count_mat[-remove_locs,]
    }
    return(count_mat)
}

