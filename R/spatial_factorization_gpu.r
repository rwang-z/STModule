# spatial factorization algorithm, GPU version for 10x Visium and high-resolution data
# based on packages torch & GPUmatrix

#' Identify tissue modules, GPU version
#' 
#' @param params parameters of the estimation process.
#' @param profile pre-processed expression matrix.
#' @param dist_mat Euclidean distance matrix of spots/cells calculated from location information. 
#' @param maxiter number of maximum iteration.
#' @param track number of iterations to track the changes of tissue module members. 
#' @param debugging the way to output FE information for debugging
#' - 'no': not output FE information for debugging
#' - 'each_update': output FE information in each iteration
#' - '10_iter': output FE information for each 10 iteration
#' @param gpumatirx_type type of GPU matrix. 
#' - 'torch': using torch tensors.
#' - 'tensorflow': using tensorflow tensors.
#' @param pip_thresh parameter for stopping criteria. Number of changes of tissue module members.
#' @return a list of estimated results for the variables in the model.
#' @export spatial_factorization_gpu

spatial_factorization_gpu = function(params, profile, dist_mat, maxiter = 2000, track=10, debugging = 'no', 
                                      gpumatirx_type = 'torch', pip_thresh = 1){

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

        # Beta: C by 1, variance of W, gamma distribution
        list_of_vars$Beta = list(e = gpu.matrix(params$e, params$C, 1, type = gpumatirx_type),
                                  f = gpu.matrix(params$f, params$C, 1, type = gpumatirx_type),
                                  mom1 = gpu.matrix(params$e * params$f, params$C, 1, type = gpumatirx_type))

        # print('Initializing Ph, Ps, Rho')
        list_of_vars$Ph = gpu.matrix(0.5, params$C, params$L, type = gpumatirx_type)
        list_of_vars$Ps = gpu.matrix(0.5, params$C, params$L, type = gpumatirx_type)
        list_of_vars$Rho = gpu.matrix(rbeta(params$C, params$r, params$z), 1, params$C, type = gpumatirx_type)

        # X: C by L, second dimension, sparse factor matrix, spike and slab prior, normal-bernoulli distribution
        list_of_vars$X = list(gamma = gpu.matrix(0.5, params$C, params$L, type = gpumatirx_type),  
                                mom1 = gpu.matrix(0, params$C, params$L, type = gpumatirx_type),
                                mom2 = gpu.matrix(0, params$C, params$L, type = gpumatirx_type),
                                sigma = gpu.matrix(100, params$C, params$L, type = gpumatirx_type),
                                m = gpu.matrix(as.vector(list_of_vars$Beta$e * list_of_vars$Beta$f), params$C, params$L, type = gpumatirx_type))
        list_of_vars$X$mom1 = list_of_vars$X$m * list_of_vars$X$gamma
        list_of_vars$X$mom2 = (list_of_vars$X$m ^ 2 + 1 / list_of_vars$X$sigma) * list_of_vars$X$gamma
        
        return(list_of_vars)
    }

    check_FE_decreasing = function(FE_res, params){
        FE_res$FEold = FE_res$FEcur
        FE_res$FEcur = Free_Energy(params)
        print(paste0("Current FE: ", FE_res$FEcur))
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
                      vars$Delt$c - vars$Delt$c * vars$Delt$d/params$d + lgamma(as.matrix(vars$Delt$c)))
        
        # FE with respect to w
        Wmom2 = vars$X$gamma * ((1/vars$X$sigma) + (vars$X$m) ** 2) + (1-vars$X$gamma) * matrix(1 / vars$Beta$mom1, params$C, params$L)
        FE = FE + 0.5 * sum(
            matrix(digamma(as.matrix(vars$Beta$e)) + log(vars$Beta$f), params$C, params$L) - matrix(vars$Beta$mom1, params$C, params$L) * Wmom2
        )

        FE = FE + sum(
            -0.5 * (
                vars$X$gamma * log(vars$X$sigma) + (1 - vars$X$gamma) * log(matrix(vars$Beta$mom1, params$C, params$L))
            )
        )

        # FE with respect to Rho
        FE = FE + sum(
            (params$r - 1) * log(vars$Rho) + (params$z - 1) * log(1 - vars$Rho)
        )
        
        # FE with respect to Psi
        FE = FE + sum(
            (params$g - 1) * log(vars$Ps) + (params$h - 1) * log(1 - vars$Ps)
        )
        
        # FE with respect to Phi
        FE = FE + sum(
            vars$Ph * matrix(log(vars$Rho), params$C, params$L) + (1 - vars$Ph) * matrix(log(1 - vars$Rho), params$C, params$L)
        )
        
        # FE with respect to s
        FE = FE + sum(
            vars$X$gamma * log(vars$Ph * vars$Ps) + (1 - vars$X$gamma) * log(1 - vars$Ph * vars$Ps)
        )
        
        Xtmp = (- (1 - vars$X$gamma) * log(1 - vars$X$gamma) - vars$X$gamma * log(vars$X$gamma))
        if(any(vars$X$gamma == 0 | vars$X$gamma == 1)){
            Xtmp[vars$X$gamma == 0 | vars$X$gamma == 1] = 0
        }
        FE = FE + sum(Xtmp)
        
        # FE with respect to Beta
        FE = FE + sum(
            params$e * log(abs(vars$Beta$f)) + (params$e - vars$Beta$e) * digamma(as.matrix(vars$Beta$e)) +
                vars$Beta$e - vars$Beta$mom1 / params$f + lgamma(as.matrix(vars$Beta$e))
        )
        
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
            mu_tmp = mean_term1[, c] - tcrossprod(a_mu[,ind_list], x_product) * vars$Lam$mom1
            prec_inv = cal_mat_inv(vars$A$precision[[c]], params)
            a_mu[, c] = prec_inv %*% mu_tmp
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
                cov_dev_2 = 0.25 * dist_mat ^ 2 * kernel
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

    updateX = function(params){
        X = vars$X

        # sigma term 2
        sigma_tmp = vars$Beta$mom1 + crossprod(vars$A$mom2, vars$Lam$mom1)
        X$sigma = gpu.matrix(as.vector(sigma_tmp), params$C, params$L, type = gpumatirx_type)

        # mean term 1
        m_term1 = crossprod(vars$Lam$mom1[,1] * vars$A$mu, profile)

        for(c in c(1 : params$C)){
            ind_list = c(1 : params$C)
            ind_list = ind_list[ind_list != c]
            a_product = crossprod(vars$Lam$mom1 * vars$A$mu[, c], vars$A$mu[,ind_list])
            tmp = a_product %*% X$mom1[ind_list,]
            X$m[c,] = (m_term1[c,] - tmp) / X$sigma[c,]
            temp_prod = as.vector(vars$Ps[c,] * vars$Ph[c,])
            u = -0.5 * log(X$sigma[c,]) + 0.5 * (X$m[c,] ^ 2) * X$sigma[c,]
                + log(temp_prod) - log(1 - temp_prod)  + 0.5 * log(vars$Beta$mom1[c])
            X$gamma[c,] = 1 / (1 + exp(-u))
            X$mom1[c,] = X$m[c,] * X$gamma[c,]
            X$mom2[c,] = (1 / X$sigma[c,] + X$m[c,] ^ 2) * X$gamma[c,]
        }

        return(X)
    }

    updateLam = function(params){
        Lam = vars$Lam
        Lam$u = gpu.matrix(params$u + 0.5 * params$L, params$N, 1, type = gpumatirx_type)

        # vars$Lam$v
        component_1 = (profile - vars$A$mu %*% vars$X$mom1)^2
        component_2 = vars$A$mom2 %*% vars$X$mom2
        component_3 = vars$A$mu^2 %*% vars$X$mom1^2
        summation = component_1 + component_2 - component_3
        
        Lam$v = gpu.matrix(1.0 / (1.0 / params$v + 0.5 * rowSums(summation)), params$N, 1, type = gpumatirx_type)
        Lam$mom1 = Lam$u * Lam$v
       
        return(Lam)
    }

    updateBeta = function(params){
        Beta = vars$Beta
        Wmom2 = vars$X$gamma * (1 / vars$X$sigma  + vars$X$m ^ 2) +
                (1 - vars$X$gamma) * (matrix(1 / Beta$mom1, params$C, params$L))

        Beta$e = (params$e + params$L/2) * gpu.matrix(1, params$C, 1, type = gpumatirx_type)

        for (c in c(1 : params$C)){
            Beta$f[c] = 1 / (1 / params$f + 0.5 * sum(Wmom2[c,]))
        }
        Beta$mom1 = Beta$e * Beta$f
        return(Beta)
    }

    updateRho = function(params){
        Rho = vars$Rho
        for(c in c(1:params$C)){
            Rho[c] = (params$r - 1 + sum(vars$Ph[c,])) / (params$L + params$r + params$z -2)
        }
        
        return(Rho)
    }

    updatePhiPsi = function(params){

        Grad = function(x, y, c, l, x_gamma, rho_c){
            v1 = x_gamma / x - (1 - x_gamma) / (1 / y - x) + log(rho_c) - log(1 - rho_c)
            v2 = x_gamma / y - (1 - x_gamma)/(1 / x - y) + (params$g - 1) / y - (params$h - 1)/(1 - y)
            vec = c(v1, v2)
            return(vec)
        }

        Hess = function(X, c, l, x_gamma){
            mat = matrix(0, 2, 2)
            mat[1, 1] = - x_gamma / (X[1] ^ 2) -  (1 - x_gamma) / ((1 / X[2] - X[1]) ^ 2)
            mat[1, 2] = - (1 - x_gamma) / (1 - X[1] * X[2]) ^ 2
            mat[2, 1] = mat[1, 2]
            mat[2, 2] = - x_gamma / (X[2] ^ 2) - (1 - x_gamma) / ((1 / X[1] - X[2]) ^ 2) - (params$g - 1) / (X[2] ^ 2) - (params$h - 1) / ((1 - X[2]) ^ 2)
            return(mat)
        }

        tildeF = function(X, c, l, x_gamma, rho_c){
            FE = 0
            FE = FE + (x_gamma) * log(X[1] * X[2]) + (1 - x_gamma) * log(1 - X[1] * X[2]) +
                (params$g - 1) * log(X[2]) + (params$h - 1) * log(1 - X[2]) + X[1] * log(rho_c) + (1 - X[1]) * log(1 - rho_c)
            return(FE)
        }

        Phi = as.matrix(vars$Ph)
        Psi = as.matrix(vars$Ps)
        x_gamma = as.matrix(vars$X$gamma)
        rho = as.matrix(vars$Rho)

        #use Newton's method for finding Ph, Ps (c,l)
        Xtol = 1e-6
        ftol = 1e-17
        for(c in c(1:params$C)){
            for (l in c(1:params$L)){
                tmpY = c(Phi[c, l], Psi[c, l])
                X = tmpY
                g = Grad(tmpY[1], tmpY[2], c, l, x_gamma[c, l], rho[c])
                H = Hess(tmpY, c, l, x_gamma[c, l])
                dH = H[1, 1] * H[2, 2] - H[1, 2] * H[2, 1]
                Direction = (matrix(c(H[2, 2], - H[2, 1], - H[1, 2], H[1, 1]), 2, 2)) %*% g * (1 / dH)
                alpha = 0.1
                i = 0
                current = tildeF(X, c, l, x_gamma[c, l], rho[c])
                while(alpha ^ i > 1e-10){
                    tmpY = c(X - (alpha ^ i) * Direction)
                    if(all(tmpY > 1e-10) && all(tmpY < 1 - 1e-10)){
                        if(tildeF(tmpY, c, l, x_gamma[c, l], rho[c]) > current){
                            Phi[c, l] = tmpY[1]
                            Psi[c, l] = tmpY[2]
                            break
                        }
                    }
                    i = i + 1
                }
            }
        }
        return(list(Phi = as.gpu.matrix(Phi, type = gpumatirx_type), Psi = as.gpu.matrix(Psi, type = gpumatirx_type)))
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
        vars$A$precision[[c]] <- as.gpu.matrix(prec_mat, type = gpumatirx_type, dtype = 'float32')
    }
    remove(prec_mat)
    vars$A$mom2 = vars$A$mu ^ 2 + 1.0 / get_prec_mat_diag(vars$A$precision, params$N, params$C)
    continue = TRUE
    iteration = 1
    
    #### initialise FE ####
    FE_res = list(FEcur = -1e50, FEold = -1e50)
    trackingvec = rep(10 * track, track)

    # update variables in iterations
    print('Updating variables...')
    while(iteration <= maxiter & continue){
        print(paste0('Iteration ',iteration))

        # update Beta
        print('Updating Beta')
        vars$Beta = updateBeta(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update X
        print('Updating X')
        vars$X = updateX(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update Rho
        print('Updating Rho')
        vars$Rho = updateRho(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update Phi and Psi 
        print('Updating Phi and Psi')
        tmp = updatePhiPsi(params)
        vars$Ph = tmp$Phi
        vars$Ps = tmp$Psi
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

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

        if(debugging == '10_iter'){
            if(iteration %% 10 == 0){
                FE_res = check_FE_decreasing(FE_res, params)
            }
        }

        # evaluate whether to stop...
        PIP = round(vars$X$gamma)
        if(iteration > 1){
            indexingvar = iteration %% track
            if(indexingvar == 0){
                indexingvar = track
            }
            trackingvec[indexingvar] = sum(abs(PIP - PIP_old))
            if(mean(trackingvec) < pip_thresh){
                continue = FALSE
            } 
        }
        PIP_old = PIP
        iteration = iteration + 1
    }


    if(debugging != 'each_update'){
        FE_res$FEcur = Free_Energy(params)
    }

    vars$maximumiteration = iteration - 1
    vars$last_FE = FE_res$FEcur

    # convert to matrices
    vars = convert_res_to_matrix(vars, 'all')

    return(vars)
}








