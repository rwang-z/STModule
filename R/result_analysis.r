#' Visiualization of spatial maps
#' 
#' @param res the results of 'run_STModule' or 'run_spatial_map_estimation' running on a tissue section.
#' @param normalization normalization method of spatial maps.
#' - 'log': log-transformation by keeping the direction of activities, default.
#' - 'scale': scale the spatial map of each module to a vector with zero-mean and unit-variance.
#' @param point_size size of points in the spatial map plot
#' - ST data: point_size = 6
#' - 10x Visium data: point_size = 2
#' - Slide-seqV2 and Stereo-seq data: point_size = 0.3
#' @import ggplot2
#' @return a list of plots each for a tissue module
#' @export spatial_map_visualization

spatial_map_visualization = function(res, normalization = 'log', point_size = 6){
    # result processing
	loc_loadings = res$A$mu
    locations = res$locations
    num_module = ncol(loc_loadings)
	
	# pattern loadings
	if(normalization == 'log'){
		pattern_loadings = sign(loc_loadings) * log2(abs(loc_loadings) + 1)
	}else if(normalization == 'scale'){
		pattern_loadings = apply(loc_loadings, 2, scale)
	}

	spatial_maps = list()
	for(g in c(1:num_module)){
		loading_vec = pattern_loadings[,g]
		pattern_df = locations
		pattern_df$activity = loading_vec
		scatter <- ggplot(pattern_df, aes(x = y, y = x, color = activity)) +
                          geom_point(size = point_size) +
                          scale_color_distiller(palette = 'RdBu', direction = -1) +
                          ylim(max(locations$x), min(locations$x)) +
                          ggtitle(paste0('Module ', g)) +
                          theme_void() +
                          theme(legend.position = "right")
		spatial_maps[[g]] = scatter
	}
	return(spatial_maps)
}


#' Get associated genes of tissue modules
#' 
#' @param res the results of 'run_STModule' running on a tissue section.
#' @return a dataframe of member information of the tissue modules. Columns: 'gene', 'activity', 'module'.
#' @export get_assocaited_genes

get_assocaited_genes = function(res){
    module_loadings = t(as.matrix(res$X$mom1))
    num_module = ncol(module_loadings)
    module_res = c()
    for(module in c(1:num_module)){
        member_loading = module_loadings[,module]
        member_df = data.frame(gene = res$gene_list, 
                                activity = member_loading,
                                act_abs = abs(member_loading), 
                                act_sign = sign(member_loading))
        member_df = member_df[member_df$activity != 0,]
        member_df = member_df[order(member_df$act_abs, decreasing = TRUE), ]
        member_df$module = module
        module_res = rbind(module_res, member_df[, c('gene', 'activity', 'module')])
    }
    return(module_res)
}


#' Visulization of associated genes of tissue modules
#' 
#' @param res the results of 'run_STModule' running on a tissue section.
#' @param quantile_thresh names will be demonstrated for genes with activity higher than the quantile threshold.
#' @import ggplot2
#' @import ggrepel
#' @return a list of plots each for a tissue module
#' @export associated_gene_visualization

associated_gene_visualization = function(res, quantile_thresh = 0.75){
    module_loadings = t(as.matrix(res$X$mom1))
    num_module = ncol(module_loadings)
    module_plots = list()
    for(module in c(1:num_module)){
        member_loading = module_loadings[,module]
        member_df = data.frame(gene = res$gene_list, 
                                activity = member_loading,
                                act_abs = abs(member_loading), 
                                act_sign = sign(member_loading))
        member_df = member_df[member_df$activity != 0,]
        quantile_value = quantile(member_df$act_abs, quantile_thresh)
        p <- ggplot(member_df, aes(x = gene, y = act_abs, color = activity)) +
            geom_point(aes(size = act_abs)) + 
            geom_text_repel(aes(label = ifelse(act_abs > quantile_value, as.character(gene),''))) +
            scale_color_gradient2(low = '#0f388a', mid = '#C0C0C0', high = '#a11212') +
            ggtitle(paste0('Module ', module)) +
            theme_classic() + 
            theme(legend.position = "none",
                  axis.title.x = element_blank(), 
                  axis.title.y = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
        module_plots[[module]] = p
    }
    return(module_plots)
}


#' Plot spatial expression of interested genes
#' 
#' @param count_file path and file name of the count matrix
#' @param loc_file path and file name of the spatial information
#' @param gene_list genes to visualize
#' @param file_sep the field separator character.
#' @param point_size size of points in the spatial expression plot
#' - ST data: point_size = 6
#' - 10x Visium data: point_size = 2
#' - Slide-seqV2 and Stereo-seq data: point_size = 0.3
#' @import ggplot2
#' @import viridis
#' @return pdf files, each for the spatial expression of a query gene
#' @export spatial_expression_visualization

spatial_expression_visualization = function(count_file, loc_file, gene_list, file_sep = '\t', point_size = 6){
    # load data
    print('Loading data...')
    count_mat = read.table(count_file, header = TRUE, sep = file_sep, row.names = 1)
    locations = read.table(loc_file, sep = file_sep, header = TRUE, row.names = 1)
    count_mat = count_mat[rownames(locations),]

    # spatial expression plots
    plots = list()
    for(gene in gene_list){
        if(gene %in% colnames(count_mat)){
            gene_exp = count_mat[, gene]
            if(sum(gene_exp) > 0){
                print(paste0('Plotting for gene ', gene))
                plot_file = paste0('spatial_expression_', gene, '.pdf')
                pdf(plot_file)
                pattern_df = locations
                pattern_df$expression = log2(gene_exp + 1)
                scatter <- ggplot(pattern_df, aes(x = y, y = x, color = expression)) +
                            geom_point(size = point_size) +
                            scale_color_viridis(option = 'D', direction = 1) +
                            ylim(max(locations$x), min(locations$x)) +
                            ggtitle(gene) +
                            theme_void() +
                            theme(legend.position = 'right',
                                  plot.title = element_text(size = 30))
                print(scatter)
                dev.off()
            }else{
                print(paste0(gene, ' is not expressed in any spot'))
            }
        }else{
            print(paste0(gene, ' is not detected in this dataset'))
        }
    }
}






