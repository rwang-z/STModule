# functions for result analysis

library(ggplot2)
library(ggrepel)
library(viridis)

spatial_map_visualization = function(res, normalization = 'log', point_size = 6){
    # return plots of spatial maps of tissue modules

    ###### parameters ######
    # res: result of function run_STModule
    # normalization: normalization method of spatial maps
    #   - 'log': log-transformation while keeping the direction of activities, default
    #   - 'scale': scale the spatial map of each module to a vector with zero-mean and unit-variance
    # point_size: size of points in the plots. Recommended point size for some SRT technologies:
    #   - ST data: point_size = 6
    #   - 10x Visium data: point_size = 2
    #   - Slide-seqV2 and Stereo-seq data: point_size = 0.3

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

get_assocaited_genes = function(res){
    # return a dataframe with columns 'gene', 'activity', 'module'

    ###### parameters ######
    # res: the result of 'run_STModule' running on a tissue section

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

associated_gene_visualization = function(res, quantile_thresh = 0.75){
    # return plots demonstraing gene activities of each tissue module

    ###### parameters ######
    # res: the result of 'run_STModule' running on a tissue section
    # quantile_thresh: names will be demonstrated for genes with activity higher than the quantile threshold

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

spatial_expression_visualization = function(count_file, loc_file, gene_list, file_sep = '\t', point_size = 6){
    # visualize spatial expression of the selected genes in the raw data
    # plot files to the folder 'plots'

    ###### parameters ######
    # count_file: path and file name of the count matrix
    # loc_file: path and file name of the spatial information
    # gene_list: genes to visualize
    # file_sep: the field separator character, default '\t'
    # point_size: size of points in the plots. Recommended point size for some SRT technologies:
    #   - ST data: point_size = 6
    #   - 10x Visium data: point_size = 2
    #   - Slide-seqV2 and Stereo-seq data: point_size = 0.3

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
                plot_file = paste0('plots/spatial_expression_', gene, '.pdf')
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






