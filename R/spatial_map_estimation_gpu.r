#' Estimating spatial maps of tissue modules on another tissue section, GPU version
#' 
#' @param params parameters of the estimation process.
#' @param profile pre-processed expression matrix.
#' @param dist_mat Euclidean distance matrix of spots/cells calculated from location information. 
#' @param module_loadings members of tissue modules.
#' @param maxiter number of maximum iteration.
#' @param track number of iterations to track the changes of tissue module members. 
#' @param debugging the way to output FE information for debugging
#' - 'no': not output FE information for debugging
#' - 'each_update': output FE information in each iteration
#' - 'iter_count': output FE information for each 'iter_count' iteration
#' @param iter_count number of iterations to check FE for debugging mode of 'iter_count'.
#' @param gpumatirx_type type of GPU matrix. 
#' - 'torch': using torch tensors.
#' - 'tensorflow': using tensorflow tensors.
#' @param pip_thresh parameter for stopping criteria. Number of changes of tissue module members.
#' @return a list of estimated results for the variables in the model.
#' @export spatial_map_estimation_gpu

spatial_map_estimation_gpu = function(params, profile, dist_mat, module_loadings, maxiter = 2000, track = 10, 
                                       debugging = 'iter_count', iter_count = 1, gpumatirx_type = 'torch', pip_thresh = 1){

    initialise_vars = function(params, dist_mat){
        list_of_vars = list()
        list_of_vars$Error = 0

        # A: N by C, first dimension, normal distribution
        list_of_vars$A = list(mu = gpu.matrix(rnorm(params$N * params$C), params$N, params$C, type = gpumatirx_type),
                               prec_term_2 = list())
        
        # R: C by 1, length scale
        if(params$fixed_r){
            list_of_vars$R = list(r = gpu.matrix(params$r_const, params$C, 1, type = gpumatirx_type),
                                   cov_inv = cal_cov_inv_mat(params$r_const, dist_mat, params))
        }else{
            r_initial = gpu.matrix(runif(params$C, 0.5, 2), params$C, 1, type = gpumatirx_type)
            list_of_vars$R = list(r = r_initial,
                                   cov_inv = cal_cov_inv_list(r_initial, dist_mat, params))
        }
        
        # Delta: C by 1, parameter of covariance matrix of MVN, gamma distribution
        list_of_vars$Delt = list(c = gpu.matrix(params$c, params$C, 1, type = gpumatirx_type),
                                  d = gpu.matrix(params$d, params$C ,1, type = gpumatirx_type))
        list_of_vars$Delt$mom1 = list_of_vars$Delt$c * list_of_vars$Delt$d

        # Lamda: N by 1, noise variance, gamma distribution. Expression in a spot share the same term.
        list_of_vars$Lam = list(u = gpu.matrix(params$u, params$N, 1, type = gpumatirx_type),
                                 v = gpu.matrix(params$v, params$N, 1, type = gpumatirx_type))
        list_of_vars$Lam$mom1 = list_of_vars$Lam$u * list_of_vars$Lam$v

        # X: C by L, second dimension, fixed using tissue modules identified from other tissues
        list_of_vars$X = list(mom1 = as.gpu.matrix(module_loadings, type = gpumatirx_type),
                               mom2 = as.gpu.matrix(module_loadings^2, type = gpumatirx_type))
        
        return(list_of_vars)
    }

    check_FE_decreasing = function(FE_res, params){
        FE_res$FEold = FE_res$FEcur
        FE_res$FEcur = Free_Energy(params)
        # if(FE_res$FEcur < FE_res$FEold){
        #     print("Decreased NFE in updating the variable")
        # }
        return(FE_res)
    }

    Free_Energy = function(params){
        #returns negative free energy 
        FE = 0

        ######## FE with respect to exp_mat
        # E[log(lamda)]
        FE = FE + 0.5 * params$L * sum(digamma(as.matrix(vars$Lam$u)) + log(vars$Lam$v))

        # E[log(p(y,params))]
        summation = (profile - vars$A$mu %*% vars$X$mom1)^2
        summation = summation + vars$A$mom2 %*% vars$X$mom2
        summation = summation - vars$A$mu^2 %*% vars$X$mom1^2
        FE = FE - 0.5 * sum(crossprod(vars$Lam$mom1, summation))

        ######## the terms from the prior and approx posteriors
        # FE with respect to A (P in the model)
        FE = FE + 0.5 * params$N * sum(digamma(as.matrix(vars$Delt$c)) + log(vars$Delt$d))
        temp_prod = 0
        det_prec = 0
        for(c in c(1 : params$C)){
            temp_prod = temp_prod + vars$Delt$mom1[c] * (crossprod(vars$A$mu[, c], vars$R$cov_inv) %*% vars$A$mu[, c])
            det_prec = det_prec + cal_precision_log_det(vars$A$precision[[c]])
        }
        FE = FE - 0.5 * temp_prod - 0.5 * det_prec

        # FE with respect to r
        FE = FE + (params$a - 1) * sum(log(vars$R$r)) - sum(vars$R$r) / params$b
        
        # FE with respect to Delta
        FE = FE + sum((params$c - vars$Delt$c) * digamma(as.matrix(vars$Delt$c)) + params$c * log(vars$Delt$d) + 
                      vars$Delt$c - vars$Delt$c * vars$Delt$d / params$d + lgamma(as.matrix(vars$Delt$c)))
        
        # FE with respect to Lamda
        FE = FE + sum((params$u - vars$Lam$u) * digamma(as.matrix(vars$Lam$u)) + params$u * log(vars$Lam$v) + 
                      vars$Lam$u - vars$Lam$u * vars$Lam$v / params$v + lgamma(as.matrix(vars$Lam$u)))
        
        return(as.numeric(FE))
    }

    updateA_precision = function(c){
        prec_mat = vars$Delt$mom1[c] * vars$R$cov_inv
        diag(prec_mat) = diag(prec_mat) + vars$Lam$mom1 * sum(vars$X$mom2[c,])

        return(prec_mat)
    }

    updateA_mu = function(params){
        a_mu = vars$A$mu
        mean_term1 = tcrossprod(as.vector(vars$Lam$mom1) * profile, vars$X$mom1)

        # the second term of vars$A$mu
        for(c in c(1 : params$C)){
            ind_list = c(1 : params$C)
            ind_list = ind_list[ind_list != c]
            x_product = tcrossprod(vars$X$mom1[c,], vars$X$mom1[ind_list,])
            mu_tmp = mean_term1[, c] - tcrossprod(a_mu[, ind_list], x_product) * vars$Lam$mom1
            prec_inv = cal_mat_inv(vars$A$precision[[c]], params)
            a_mu[,c] = prec_inv %*% mu_tmp
        }
       
        return(a_mu)
    }

    updateA_mom2 = function(params){
        a_mom2 = vars$A$mu ^ 2 + 1.0 / get_prec_mat_diag(vars$A$precision, params$N, params$C)

        return(a_mom2)
    }

    updateR = function(params){
        r_tildeF = function(r, c, cov_inv, mu_c){
            temp = crossprod(mu_c, cov_inv) %*% mu_c
            FE = - 0.5 * vars$Delt$mom1[c] * temp
            FE = FE + (params$a - 1) * log(r) - r / params$b
            return(FE)
        }

        # update r and cov_inv of vars$R
        R = vars$R
        max_iter = 10
        for(c in c(1 : params$C)){
            iter = 1
            r1 = R$r[c]
            r_old = R$r[c]
            cov_inv = R$cov_inv[[c]]
            mu_c = vars$A$mu[, c]
            FE_old = r_tildeF(r1, c, cov_inv, mu_c)
            while(iter <= max_iter && r1 > 0){
                r0 = r1
                kernel = exp(-0.5 * r0 * dist_mat)
                cov_dev = -0.5 * dist_mat * kernel
                prod_cov_inv_dev = cov_inv %*% cov_dev
                tmp_1 = prod_cov_inv_dev %*% cov_inv
                # r_dev
                tmp_2 = crossprod(mu_c,  - tmp_1) %*% mu_c
                r_dev = - 0.5 * vars$Delt$mom1[c] * tmp_2 + (params$a - 1) / r0 - 1.0 / params$b
                # r_dev_2
                cov_dev_2 = 0.25 * dist_mat^2 * kernel
                cov_inv_dev_2 = 2 * prod_cov_inv_dev %*% tmp_1
                cov_inv_dev_2 = cov_inv_dev_2 - cov_inv %*% cov_dev_2 %*% cov_inv
                tmp_3 = crossprod(mu_c, cov_inv_dev_2) %*% mu_c
                r_dev_2 = - 0.5 * vars$Delt$mom1[c] * tmp_3
                r_dev_2 = r_dev_2 - (params$a - 1) / (r0^2) 
                # updated r
                r1 = r0 - (r_dev[1] / r_dev_2[1])
                if(r1 <= 0){
                    break
                }
                cov_inv = cal_cov_inv_mat(r1, dist_mat, params)
                FE_new = r_tildeF(r1, c, cov_inv, mu_c)
                if(FE_new > FE_old){
                    R$r[c] = r1
                    R$cov_inv[[c]] = cov_inv
                    break
                }
                iter = iter + 1
            }
        }
        return(R)
    }

    updateDelt = function(params){
        Delt = vars$Delt
        Delt$c = gpu.matrix(params$c + 0.5 * params$N, params$C, 1, type = gpumatirx_type)

        # vars$Delt$d
        for(c in c(1 : params$C)){
            temp = sum(crossprod(vars$A$mu[, c], vars$R$cov_inv) %*% vars$A$mu[, c])
            Delt$d[c] = 1.0 / (1.0 / params$d + 0.5 * temp)
        }

        Delt$mom1 = Delt$c * Delt$d
        return(Delt)
    }

    updateLam=function(params){
        Lam = vars$Lam
        Lam$u = gpu.matrix(params$u + 0.5 * params$L, params$N, 1, type = gpumatirx_type)

        # vars$Lam$v
        component_1 = (profile - vars$A$mu %*% vars$X$mom1)^2
        component_2 = vars$A$mom2 %*% vars$X$mom2
        component_3 = vars$A$mu^2 %*% vars$X$mom1^2
        summation = component_1 + component_2 - component_3
        
        Lam$v = gpu.matrix(1.0 / (1.0/params$v + 0.5 * rowSums(summation)), params$N, 1, type = gpumatirx_type)
        Lam$mom1 = Lam$u * Lam$v
       
        return(Lam)
    }


    #### initialise variables ####
    print('Using STModule GPU version')
    rownames(profile) = c()
    colnames(profile) = c()
    rownames(dist_mat) = c()
    colnames(dist_mat) = c()
    dist_mat = as.gpu.matrix(dist_mat, type = gpumatirx_type, dtype = 'float32')
    profile = as.gpu.matrix(profile, type = gpumatirx_type, dtype = 'float32')

    #### initialise variables ####
    print('Initializing variables ...')
    vars = initialise_vars(params, dist_mat)
    for(c in c(1 : params$C)){
        prec_mat = vars$Delt$mom1[c] * vars$R$cov_inv
        diag(prec_mat) = diag(prec_mat) + vars$Lam$mom1 * sum(vars$X$mom2[c,])
        vars$A$precision[[c]] = as.gpu.matrix(prec_mat, type = gpumatirx_type, dtype = 'float32')
    }
    remove(prec_mat)
    vars$A$mom2 = vars$A$mu^2 + 1.0 / get_prec_mat_diag(vars$A$precision, params$N, params$C)
    continue = TRUE
    iteration = 1
    
    #### initialise FE ####
    FE_res = list(FEcur = -1e40, FEold = -1e50)
    trackingvec = rep(10 * track, track)

    # update variables in iterations
    print('Updating variables...')
    while(iteration <= maxiter & continue){
        print(paste0('Iteration ',iteration))

        # update A
        print('Updating A')
        for(c in c(1 : params$C)){
            vars$A$precision[[c]] = updateA_precision(c)
        }
        vars$A$mu = updateA_mu(params)
        vars$A$mom2 = updateA_mom2(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update r
        if(!params$fixed_r){
            print('Updating R')
            vars$R = updateR(params)
            if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}
        }
        
        # update Delta
        print('Updating Delta')
        vars$Delt = updateDelt(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}
        
        # update Lamda
        print('Updating Lamda')
        vars$Lam = updateLam(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # check FE each 10 iterations
        if(debugging == 'iter_count'){
            if(iteration %% iter_count == 0){
                FE_res = check_FE_decreasing(FE_res, params)
            }
        }

        # evaluate whether to stop...
        if(iteration > 1){
            if(abs(FE_res$FEcur - FE_res$FEold) < 1e-6){continue = FALSE} 
        }
        iteration = iteration + 1
    }
    FE_res$FEcur = Free_Energy(params)
    vars$maximumiteration = iteration - 1
    vars$last_FE = FE_res$FEcur

    # convert to matrices
    vars = convert_res_to_matrix(vars, 'estimate')

    return(vars)
}








