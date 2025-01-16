# functions

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

cal_cov_inv_mat = function(r, dist_mat, params){
    cov_mat = exp(-0.5 * r * dist_mat)
    cov_inv = cal_mat_inv(cov_mat, params)
    return(cov_inv)
}

cal_cov_inv_list = function(r, dist_mat, params){
    num_comp = dim(r)[1]
    cov_list = list()
    for(c in c(1:num_comp)){
        cov_list[[c]] = cal_cov_inv_mat(r[c], dist_mat, params)
    }
    return(cov_list)
}

cal_cov_inv_list_same = function(r_const, num_comp, dist_mat, params){
    cov_inv = cal_cov_inv_mat(r_const, dist_mat, params)
    cov_list = list()
    for(c in c(1:num_comp)){
        cov_list[[c]] = cov_inv
    }
    return(cov_list)
}

cal_precision_log_det = function(prec_mat){
    det_log = tryCatch(determinant(prec_mat)$modulus[1],
                        error = function(e){
                            return(2 * determinant(chol(prec_mat))$modulus[1])
                        })
    return(det_log)
}

get_prec_mat_diag = function(prec_list, dimension, num_comp){
    prec_diag_mat = matrix(0, dimension, num_comp)
    for(c in c(1:num_comp)){
        prec_diag_mat[, c] = diag(prec_list[[c]])
    }
    return(prec_diag_mat)
}

cal_distance_mat = function(locations, dimension = '2D'){
    # locations: dataframe including columns 'x_pos' and 'y_pos' to calcualte the distances
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

data_filtering = function(count_mat, high_resolution, gene_filtering, cell_thresh = 100){
    # remove genes expressed in less than gene_filtering locations/cells
    # remove locations with zero expression, or cells with less than cell_thresh total counts
    if(high_resolution){
        # gene_filtering: number of locations/cells
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

remove_zero_exp_locations = function(count_mat){
    loc_sum = rowSums(count_mat)
    remove_locs = which(loc_sum == 0)
    if(length(remove_locs) > 0){
        count_mat = count_mat[-remove_locs,]
    }
    return(count_mat)
}

