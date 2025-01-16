# spatial factorization algorithm

spatial_factorization <- function(params, profile, dist_mat, maxiter = 2000, track=10, 
                                  debugging = 'no', pip_thresh = 1){

    initialise_vars <- function(params, dist_mat){
        list_of_vars <- list()
        list_of_vars$Error=0

        # A: N by C, first dimension, normal distribution
        list_of_vars$A <- list(mu = matrix(rnorm(params$N * params$C),params$N,params$C), 
                               prec_term_2 = list())

        # R: C by 1, length scale
        if(params$fixed_r){
            list_of_vars$R <- list(r = matrix(params$r_const, params$C, 1),
                                    cov_inv = cal_cov_inv_list_same(params$r_const, params$C, dist_mat, params))
        }else{
            r_initial = matrix(runif(params$C, 0.5, 2), params$C, 1)
            list_of_vars$R <- list(r = r_initial,
                                    cov_inv = cal_cov_inv_list(r_initial, dist_mat, params))
        }
        
        # Delta: C by 1, parameter of covariance matrix of MVN, gamma distribution
        list_of_vars$Delt <- list(c = matrix(params$c,params$C,1),
                                  d = matrix(params$d,params$C,1))
        list_of_vars$Delt$mom1 = list_of_vars$Delt$c * list_of_vars$Delt$d

        # Lamda: N by 1, noise variance, gamma distribution. Expression in a spot share the same term.
        list_of_vars$Lam <- list(u = matrix(params$u,params$N,1),
                                 v = matrix(params$v,params$N,1))
        list_of_vars$Lam$mom1 = list_of_vars$Lam$u * list_of_vars$Lam$v

        # Beta: C by 1, variance of W, gamma distribution
        list_of_vars$Beta <- list(e = matrix(params$e,params$C,1),
                                  f = matrix(params$f,params$C,1),
                                  mom1 = matrix(params$e * params$f,params$C,1))

        list_of_vars$Ph <- matrix(0.5,params$C,params$L)
        list_of_vars$Ps <- matrix(0.5,params$C,params$L)
        list_of_vars$Rho <- matrix(rbeta(params$C,params$r,params$z),1,params$C)

        # X: C by L, second dimension, sparse factor matrix, spike and slab prior, normal-bernoulli distribution
        list_of_vars$X <- list(gamma = matrix(0.5,params$C,params$L),  
                                mom1 = matrix(0,params$C,params$L),
                                mom2 = matrix(0,params$C,params$L),
                                sigma = matrix(c(100),params$C,params$L),
                                m = matrix(list_of_vars$Beta$e *list_of_vars$Beta$f,params$C,params$L))
        list_of_vars$X$mom1 = list_of_vars$X$m * list_of_vars$X$gamma
        list_of_vars$X$mom2 = (list_of_vars$X$m**2 + 1/list_of_vars$X$sigma) * list_of_vars$X$gamma
        
        return(list_of_vars)
    }

    check_FE_decreasing = function(FE_res, params){
        FE_res$FEold = FE_res$FEcur
        FE_res$FEcur = Free_Energy(params)
        print(paste0("Current FE: ", FE_res$FEcur))
        if(FE_res$FEcur < FE_res$FEold){
            print("Decreased NFE in updating the variable")
            if(params$decrease_stop){
                break
            }
        }
        return(FE_res)
    }

    Free_Energy<-function(params){
        #returns negative free energy 
        FE = 0

        ######## FE with respect to exp_mat
        # E[log(lamda)]
        FE = FE + 0.5 * params$L * sum(digamma(vars$Lam$u)+log(vars$Lam$v))

        # E[log(p(y,params))]
        component_1 = (profile - vars$A$mu %*% vars$X$mom1)^2
        component_2 = vars$A$mom2 %*% vars$X$mom2
        component_3 = vars$A$mu^2 %*% vars$X$mom1^2
        summation = component_1 + component_2 - component_3 
        FE = FE - 0.5 * sum(matrix(vars$Lam$mom1, params$N, params$L) * summation) 

        ######## the terms from the prior and approx posteriors
        # FE with respect to A (P in the model)
        FE = FE + 0.5 * params$N * sum(digamma(vars$Delt$c)+log(vars$Delt$d))
        temp_prod = 0
        det_prec = 0
        for(c in c(1:params$C)){
            temp_prod = temp_prod + vars$Delt$mom1[c] * (crossprod(vars$A$mu[, c], vars$R$cov_inv[[c]]) %*% vars$A$mu[, c])
            det_prec = det_prec + cal_precision_log_det(vars$A$precision[[c]])
        }
        FE = FE - 0.5 * temp_prod - 0.5 * det_prec

        # FE with respect to r
        FE = FE + (params$a - 1) * sum(log(vars$R$r)) - sum(vars$R$r) / params$b
       
        # FE with respect to Delta
        FE = FE + sum((params$c - vars$Delt$c)*digamma(vars$Delt$c) + params$c*log(vars$Delt$d) + 
                      vars$Delt$c - vars$Delt$c*vars$Delt$d/params$d + lgamma(vars$Delt$c))
        
        # FE with respect to w
        Wmom2 = vars$X$gamma * ((1/vars$X$sigma) + (vars$X$m)**2) + (1-vars$X$gamma) * matrix(1/vars$Beta$mom1,params$C,params$L)
        FE = FE + 0.5*sum(
            matrix(digamma(vars$Beta$e) + log(vars$Beta$f),params$C,params$L) - matrix(vars$Beta$mom1,params$C,params$L)*Wmom2
        )
        FE = FE + sum(
            -0.5*(
                vars$X$gamma * log(vars$X$sigma) + (1-vars$X$gamma)*log(matrix(vars$Beta$mom1,params$C,params$L))
            )
        )
        
        # FE with respect to Rho
        FE = FE + sum(
            (params$r - 1)*log(vars$Rho) + (params$z - 1)*log(1-vars$Rho)
        )

        # FE with respect to Psi
        FE = FE + sum(
            (params$g - 1) * log(vars$Ps) + (params$h - 1)*log(1 - vars$Ps)
        )
       
        # FE with respect to Phi
        FE = FE + sum(
            vars$Ph * matrix(log(vars$Rho),params$C,params$L) + (1-vars$Ph)*matrix(log(1-vars$Rho),params$C,params$L)
        )
        
        # FE with respect to s
        FE = FE + sum(
            vars$X$gamma*log(vars$Ph*vars$Ps) + (1-vars$X$gamma)*log(1-vars$Ph*vars$Ps)
        )
        Xtmp<- (-(1-vars$X$gamma)*log(1-vars$X$gamma) - vars$X$gamma*log(vars$X$gamma))
        if(any(vars$X$gamma==0 | vars$X$gamma==1)){
            Xtmp[vars$X$gamma==0 | vars$X$gamma==1]=0
        }
        FE = FE + sum(Xtmp)
        
        # FE with respect to Beta
        FE = FE + sum(
            params$e * log(abs(vars$Beta$f)) + (params$e - vars$Beta$e)*digamma(vars$Beta$e) +
                vars$Beta$e - vars$Beta$mom1/params$f + lgamma(vars$Beta$e)
        )
       
        # FE with respect to Lamda
        FE = FE + sum((params$u - vars$Lam$u)*digamma(vars$Lam$u) + params$u*log(vars$Lam$v) + 
                      vars$Lam$u - vars$Lam$u*vars$Lam$v/params$v + lgamma(vars$Lam$u))
        
        return(FE)
    }

    updateA=function(params){
        A=vars$A 

        # precision
        prec_inv = list()
        for(c in 1:params$C){
            prec_term_1 = matrix(0, params$N, params$N)
            diag(prec_term_1) = vars$Lam$mom1 * sum(vars$X$mom2[c,])
            prec_term_2 = vars$Delt$mom1[c] * vars$R$cov_inv[[c]]
            prec_mat = prec_term_1 + prec_term_2
            A$precision[[c]] = prec_mat
            prec_inv[[c]] = cal_mat_inv(prec_mat, params)
        }

        # the first term of vars$A$mu
        mean_term1 = tcrossprod(as.vector(vars$Lam$mom1) * profile, vars$X$mom1)
       
        # the second term of vars$A$mu
        for(c in 1:params$C){
            x_product = colSums(vars$X$mom1[c,] * t(vars$X$mom1[-c,]))
            mu_tmp = mean_term1[,c] - (A$mu[,-c] %*% x_product) * vars$Lam$mom1
            A$mu[,c] = prec_inv[[c]] %*% mu_tmp
        }

        A$mom2 = A$mu^2 + 1.0 / get_prec_mat_diag(A$precision, params$N, params$C)

        return(A)
    }

    updateR = function(params){
        cal_r_dev = function(r, c, cov_inv){
            cov_dev = -0.5 * dist_mat * exp(-0.5 * r * dist_mat)
            cov_inv_dev = - cov_inv %*% cov_dev %*% cov_inv
            temp = crossprod(vars$A$mu[, c],  cov_inv_dev) %*% vars$A$mu[, c]
            dev = - 0.5 * vars$Delt$mom1[c] * temp
            dev = dev + (params$a - 1) / r - 1.0 / params$b
            return(dev[1])
        }

        cal_r_dev_2 = function(r, c, cov_inv){
            cov_dev = -0.5 * dist_mat * exp(-0.5 * r * dist_mat)
            cov_dev_2 = 0.25 * dist_mat^2 * exp(-0.5 * r * dist_mat)
            cov_inv_dev_2 = 2 * cov_inv %*% cov_dev %*% cov_inv %*% cov_dev %*% cov_inv
            cov_inv_dev_2 = cov_inv_dev_2 - cov_inv %*% cov_dev_2 %*% cov_inv
            temp = crossprod(vars$A$mu[, c], cov_inv_dev_2) %*% vars$A$mu[, c]
            dev = - 0.5 * vars$Delt$mom1[c] * temp
            dev = dev - (params$a - 1) / (r^2) 
            return(dev[1])
        }

        tildeF = function(r, c, cov_inv){
            temp = crossprod(vars$A$mu[, c], cov_inv) %*% vars$A$mu[, c]
            FE = - 0.5 * vars$Delt$mom1[c] * temp
            FE = FE + (params$a - 1) * log(r) - r / params$b
            return(FE)
        }

        # update r and cov_inv of vars$R
        R = vars$R
        max_iter = 10
        for(c in c(1:params$C)){
            iter = 1
            r1 = R$r[c]
            r_old = R$r[c]
            cov_inv = R$cov_inv[[c]]
            FE_old = tildeF(r1, c, cov_inv)
            while(iter <= max_iter && r1 > 0){
                r0 = r1
                r1 = r0 - (cal_r_dev(r0, c, cov_inv) / cal_r_dev_2(r0, c, cov_inv))
                if(r1 <= 0){
                    break
                }
                cov_inv = cal_cov_inv_mat(r1, dist_mat, params)
                FE_new = tildeF(r1, c, cov_inv)
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
        Delt=vars$Delt
        Delt$c = params$c + 0.5 * params$N

        # vars$Delt$d
        for(c in c(1:params$C)){
            temp = sum(crossprod(vars$A$mu[, c], vars$R$cov_inv[[c]]) %*% vars$A$mu[, c])
            Delt$d[c] = 1.0 / (1.0/params$d + 0.5 * temp)
        }
        Delt$mom1 = Delt$c * Delt$d

        return(Delt)
    }

    updateX=function(params){
        X=vars$X

        # sigma term 2
        sigma_tmp = vars$Beta$mom1 + crossprod(vars$A$mom2, vars$Lam$mom1)
        X$sigma = matrix(sigma_tmp, params$C, params$L)

        # mean term 1
        m_term1 = crossprod(vars$A$mu, (as.vector(vars$Lam$mom1) * profile))

        for(c in 1:params$C){
            a_product = colSums(as.vector(vars$Lam$mom1) * (vars$A$mu[,c] * vars$A$mu[,-c])) 
            tmp = a_product %*% X$mom1[-c,] 
            X$m[c,] = (m_term1[c,] - tmp) / X$sigma[c,]  

            u = -0.5*log(X$sigma[c,]) + 0.5 * (X$m[c,]^2) * X$sigma[c,] 
                + log(vars$Ps[c,]*vars$Ph[c,]) + 0.5*log(vars$Beta$mom1[c]) - log(1-vars$Ps[c,]*vars$Ph[c,])

            X$gamma[c,] = 1/(1+exp(-u))
            X$mom1[c,] = X$m[c,] * X$gamma[c,]
            X$mom2[c,] = (1 / X$sigma[c,] + X$m[c,]^2) * X$gamma[c,]
        } 

        return(X)
    }

    updateLam=function(params){
        Lam=vars$Lam
        Lam$u = params$u + 0.5 * params$L

        # vars$Lam$v
        component_1 = (profile - vars$A$mu %*% vars$X$mom1)^2
        component_2 = vars$A$mom2 %*% vars$X$mom2
        component_3 = vars$A$mu^2 %*% vars$X$mom1^2
        summation = component_1 + component_2 - component_3
        
        Lam$v = 1.0 / (1.0/params$v + 0.5 * rowSums(summation))
        Lam$mom1 = Lam$u * Lam$v

        return(Lam)
    }

    updateBeta=function(params){
        Beta=vars$Beta
        Wmom2 = vars$X$gamma * (1/vars$X$sigma  + vars$X$m^2) +
            (1-vars$X$gamma)*(matrix(1/Beta$mom1,params$C,params$L))
        Beta$e = (params$e + params$L/2)*matrix(1,params$C,1)
        for (c in 1:params$C){
            Beta$f[c] = 1/(1/params$f + 0.5*sum(Wmom2[c,]))
        }
        Beta$mom1 = Beta$e*Beta$f

        return(Beta)
    }

    updateRho=function(params){
        Rho=vars$Rho
        for(c in 1:params$C){
            Rho[c] = (params$r - 1 + sum(vars$Ph[c,]))/(params$L + params$r + params$z -2)
        }
        return(Rho)
    }

    updatePhiPsi=function(params){
        Grad=function(x,y,c,l){
            v1 = vars$X$gamma[c,l]/x - (1-vars$X$gamma[c,l])/(1/y -x) + log(vars$Rho[c]) - log(1-vars$Rho[c])
            v2 = vars$X$gamma[c,l]/y - (1-vars$X$gamma[c,l])/(1/x - y) + (params$g - 1)/y - (params$h - 1)/(1-y)
            vec=c(v1,v2)
            return(vec)
        }
        Hess = function(X,c,l){
            mat=matrix(0,2,2)
            mat[1,1] = -vars$X$gamma[c,l]/(X[1]^2) -  (1-vars$X$gamma[c,l])/((1/X[2] - X[1])^2)
            mat[1,2] = -(1-vars$X$gamma[c,l])/(1-X[1]*X[2])^2
            mat[2,1] = mat[1,2]
            mat[2,2] = -vars$X$gamma[c,l]/(X[2]^2) - (1 - vars$X$gamma[c,l])/((1/X[1] - X[2])^2) - (params$g - 1)/(X[2]^2) - (params$h - 1)/((1-X[2])^2)
            return(mat)
        }
        tildeF = function(X,c,l){
            FE =0
            FE = FE + (vars$X$gamma[c,l])*log(X[1]*X[2]) + (1-vars$X$gamma[c,l])*log(1-X[1]*X[2]) +
                (params$g - 1)*log(X[2]) + (params$h -1)*log(1-X[2]) + X[1]*log(vars$Rho[c]) +(1-X[1])*log(1-vars$Rho[c])
            return(FE)
        }
        Phi=vars$Ph
        Psi=vars$Ps
        #use Newton's method to update Ph, Ps
        Xtol = 1e-6
        ftol = 1e-17
        for(c in 1:params$C){
            for (l in 1:params$L){
                tmpY = c(Phi[c,l],Psi[c,l])
                X=tmpY
                g = Grad(tmpY[1],tmpY[2],c,l)
                H=Hess(tmpY,c,l)
                dH = H[1,1]*H[2,2] - H[1,2]*H[2,1]
                Direction = (matrix(c(H[2,2],-H[2,1],-H[1,2],H[1,1]),2,2))%*%g * (1/dH)
                alpha = 0.1
                i=0
                current=tildeF(X,c,l)
                while(alpha^i>1e-10){
                    tmpY=c(X - (alpha^i)*Direction)
                    if(all(tmpY>1e-10) && all(tmpY<1-1e-10)){
                        if(tildeF(tmpY,c,l)>current){
                            Phi[c,l]=tmpY[1]
                            Psi[c,l]=tmpY[2]
                            break
                        }
                    }
                    i=i+1
                }
            }

        }
        return(list(Phi=Phi,Psi=Psi))
    }

    #### initialise variables ####
    print('Initializing...')
    vars<-initialise_vars(params, dist_mat)
    for(c in c(1:params$C)){
        term_1 = matrix(0, params$N, params$N)
        diag(term_1) = vars$Lam$mom1 * sum(vars$X$mom2[c,])
        term_2 = vars$Delt$mom1[c] * vars$R$cov_inv[[c]]
        prec_mat = term_1 + term_2
        vars$A$precision[[c]] = prec_mat
    }
    vars$A$mom2 = vars$A$mu^2 + 1.0 / get_prec_mat_diag(vars$A$precision, params$N, params$C)
    vars$dist_tr = sum(diag(dist_mat))
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
        vars$Beta=updateBeta(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update X
        print('Updating X')
        vars$X=updateX(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update Rho
        print('Updating Rho')
        vars$Rho=updateRho(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update Phi and Psi 
        print('Updating Phi and Psi')
        tmp = updatePhiPsi(params)
        vars$Ph = tmp$Phi
        vars$Ps = tmp$Psi
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update A
        print('Updating A')
        vars$A=updateA(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # update r
        if(!params$fixed_r){
            print('Updating R')
            vars$R=updateR(params)
            if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}
        }

        # update Delta
        print('Updating Delta')
        vars$Delt=updateDelt(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}
        
        # update Lamda
        print('Updating Lamda')
        vars$Lam=updateLam(params)
        if(debugging == 'each_update'){FE_res = check_FE_decreasing(FE_res, params)}

        # check FE each 10 iterations
        if(debugging == '10_iter'){
            if(iteration %% 10 == 0){
                FE_res = check_FE_decreasing(FE_res, params)
            }
        }

        # evaluate whether to stop...
        PIP = round(vars$X$gamma)
        if(iteration > 1){
            indexingvar = iteration%%track
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

    return(vars)
}
