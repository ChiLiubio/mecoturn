#' @title Analyze the 'turnover' of microbial communities.
#' 
#' @description
#' Analyze the 'turnover' of microbial communities, i.e. beta-diversity along a gradient <doi:10.1111/j.1461-0248.2010.01552.x>.
#' The workflow consists of the steps of dissimilarity matrix generation, matrix conversion, differential test and visualization.
#'
#' @export
betaturn <- R6::R6Class(classname = "betaturn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param measure default "bray"; beta diversity dissimilarity metric; 
		#' 	 must be one of \code{c("bray", "jaccard", "wei_unifrac", "unwei_unifrac", "betaMPD", "betaMNTD", "betaNRI", "betaNTI", "ses_UniFrac", "RCbray")}
		#' 	 or other options in parameter \code{method} of \code{vegan::vegdist} function.
		#' 	 If the distance matrix has been in the beta_diversity list of microtable object, 
		#' 	 the function can ignore this step. Otherwise, the function can generate the corresponding beta diversity distance matrix in the microtable object.
		#' 	 bray: Bray-Curtis; RCbray: Raup–Crick based Bray-Curtis; wei_unifrac: weighted UniFrac; ses_UniFrac: standardized deviation of UniFrac.
		#' @param filter_thres default 0; the relative abundance threshold used to filter features with low abundance.
		#' @param abundance.weighted default TRUE; whether use abundance-weighted method for the phylogenetic metrics.
		#' @param null.model default NULL; one of \code{c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")},
		#'   in which "taxa.labels" can only be used for phylogenetic analysis. 
		#'   See \code{null.model} parameter of \code{ses.mntd} function in \code{picante} package for the algorithm details.
		#' @param runs default 1000; simulation number of times for null model.
		#' @param iterations default 1000; iteration number for part null models to perform; see iterations parameter of \code{picante::randomizeMatrix} function.
		#' @param ... parameters passed to \code{cal_betadiv} function of \code{microtable} class when provided measure is not in the current vector;
		#'   parameters passed to \code{cal_betamntd} function of \code{trans_nullmodel} class when \code{measure = "betaMNTD"};
		#'   parameters passed to \code{cal_ses_betamntd} function of \code{trans_nullmodel} class when \code{measure = "betaNTI"}.
		#' @return \code{dataset}, stored in the object. The new dataset has a beta_diversity list and the calculated distance matrix in the list.
		#' @examples
		#' data(wheat_16S)
		#' b1 <- betaturn$new(wheat_16S, measure = "bray")
		initialize = function(dataset, measure = "bray", filter_thres = 0, abundance.weighted = TRUE, null.model = NULL, runs = 1000, iterations = 1000, ...){
			# check input format
			if(!inherits(dataset, "microtable")){
				stop("Input dataset must be microtable class!")
			}
			tmp <- clone(dataset)
			if(filter_thres > 0){
				tmp$otu_table %<>% .[apply(., 1, sum)/sum(.) > filter_thres, ]
				tmp$tidy_dataset()
			}
			if(is.null(measure)){
				stop("Input measure should not be NULL!")
			}
			if(measure %in% names(tmp$beta_diversity)){
				message("The input measure has been in the beta_diversity list of dataset ...")
			}else{
				if(is.null(null.model)){
					if(measure == "RCbray"){
						null.model <- "independentswap"
						message("Use independentswap as the null model for ", measure, " ...")
					}else{
						if(measure %in% c("betaNRI", "betaNTI", "ses_UniFrac")){
							null.model <- "taxa.labels"
							message("Use taxa.labels as the null model for ", measure, " ...")
						}
					}
				}else{
					# check the input null model
					if(! null.model %in% c('taxa.labels', 'richness', 'frequency', 'sample.pool', 'phylogeny.pool', 'independentswap', 'trialswap')){
						stop("The input null.model must be one of c('taxa.labels', 'richness', 'frequency', 'sample.pool', 'phylogeny.pool', 'independentswap', 'trialswap')!")
					}
				}
				if(! measure %in% c("betaMPD", "betaMNTD", "betaNRI", "betaNTI", "ses_UniFrac", "RCbray")){
					if(measure %in% c("bray", "jaccard")){
						suppressMessages(tmp$cal_betadiv(method = NULL, unifrac = FALSE))
					}else{
						if(measure %in% c("wei_unifrac", "unwei_unifrac")){
							suppressMessages(tmp$cal_betadiv(unifrac = TRUE))
						}else{
							suppressMessages(tmp$cal_betadiv(unifrac = FALSE, method = measure, ...))
						}
					}
				}else{
					trans_nullmodel_object <- trans_nullmodel$new(tmp, filter_thres = 0)
					if(measure == "betaMPD"){
						trans_nullmodel_object$cal_betampd(abundance.weighted = abundance.weighted)
						res_dis <- trans_nullmodel_object$res_betampd
					}
					if(measure == "betaMNTD"){
						trans_nullmodel_object$cal_betamntd(abundance.weighted = abundance.weighted, ...)
						res_dis <- trans_nullmodel_object$res_betamntd
					}
					if(measure == "betaNRI"){
						trans_nullmodel_object$cal_ses_betampd(abundance.weighted = abundance.weighted, runs = runs, null.model = null.model, iterations = iterations)
						res_dis <- trans_nullmodel_object$res_ses_betampd
					}
					if(measure == "betaNTI"){
						trans_nullmodel_object$cal_ses_betamntd(abundance.weighted = abundance.weighted, runs = runs, null.model = null.model, iterations = iterations, ...)
						res_dis <- trans_nullmodel_object$res_ses_betamntd
					}
					if(measure == "RCbray"){
						trans_nullmodel_object$cal_rcbray(runs = runs, null.model = null.model)
						res_dis <- trans_nullmodel_object$res_rcbray
					}
					if(measure == "ses_UniFrac"){
						feature_table <- t(tmp$otu_table)
						res_dis <- private$ses_unifrac(feature_table, tmp$phylo_tree, runs = runs, weighted = abundance.weighted, null.model = null.model, 
							iterations = iterations)
					}
					tmp$beta_diversity[[measure]] <- res_dis
				}
			}
			self$measure <- measure
			self$dataset <- tmp
			message('The distance matrix is stored in object$beta_diversity list ...')
		},
		#' @description
		#' Convert sample distances within groups or between groups.
		#'
		#' @param group one colname of sample_table in \code{microtable} object used for group distance convertion.
		#' @param within_group default TRUE; whether transform sample distance within groups? If FALSE, transform sample distances between any two groups.
		#' @param by_group default NULL; one colname of sample_table in \code{microtable} object.
		#'   If provided, convert distances according to the provided by_group parameter. This is especially useful for ordering and filtering values further.
		#'   When \code{within_group = TRUE}, the result of by_group parameter is the format of paired groups.
		#'   When \code{within_group = FALSE}, the result of by_group parameter is the format same with the group information in \code{sample_table}.
		#' @param ordered_group default NULL; a vector representing the ordered elements of \code{group} parameter; only useful when within_group = FALSE.
		#' @param sep default TRUE; a character string to separate the group names after merging them into a new name.
		#' @param add_cols default NULL; add several columns of sample_table to the final \code{res_group_distance} table according to the \code{by_group} column; 
		#'   invoked only when \code{within_group = FALSE}.
		#' @return \code{res_group_distance} stored in object.
		#' @examples
		#' b1$cal_group_distance(group = "Type", within_group = FALSE, by_group = "Plant_ID")
		cal_group_distance = function(group, within_group = TRUE, by_group = NULL, ordered_group = NULL, sep = " vs ", add_cols = NULL){
			dataset <- self$dataset
			measure <- self$measure
			trans_beta_object <- trans_beta$new(dataset = dataset, group = group, measure = measure)
			suppressMessages(trans_beta_object$cal_group_distance(within_group = within_group, by_group = by_group, ordered_group = ordered_group, sep = sep))
			if(!is.null(add_cols) & within_group == FALSE){
				if(!all(add_cols %in% colnames(dataset$sample_table))){
					stop("Please provide correct add_cols parameter!")
				}
				trans_beta_object$res_group_distance %<>% dplyr::left_join(., unique(dataset$sample_table[, c(by_group, add_cols)]), by = by_group)
			}
			self$res_group_distance <- trans_beta_object$res_group_distance
			message('The result is stored in object$res_group_distance ...')
			self$tmp_trans_beta <- trans_beta_object
		},
		#' @description
		#' Differential test of distances among groups.
		#'
		#' @param ... parameters passed to \code{cal_group_distance_diff} function of \code{\link{trans_beta}} class.
		#' @return \code{res_group_distance_diff} stored in object.
		#' @examples
		#' b1$cal_group_distance_diff(method = "wilcox")
		cal_group_distance_diff = function(...){
			self$tmp_trans_beta$res_group_distance <- self$res_group_distance
			suppressMessages(self$tmp_trans_beta$cal_group_distance_diff(...))
			self$res_group_distance_diff <- self$tmp_trans_beta$res_group_distance_diff
			message('The result is stored in object$res_group_distance_diff ...')
		},
		#' @description
		#' Plot the distance between samples within or between groups.
		#'
		#' @param ... parameters passed to \code{plot_group_distance} function of \code{\link{trans_beta}} class.
		#' @return \code{ggplot}.
		#' @examples
		#' b1$plot_group_distance()
		plot_group_distance = function(...){
			self$tmp_trans_beta$plot_group_distance(...)
		}
	),
	private = list(
		# calculate standardized deviation of unifrac
		ses_unifrac = function(samp, tre, null.model, runs, weighted, iterations) {
			# samp: OTU table, rownames are sample; tre: phylonetic tree
			if (sum(!(colnames(samp) %in% tre$tip.label)) != 0) {
				stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\tin the OTU table and the tree should match!")
			}
			message("Run the observed UniFrac ...")
			betaobs <- private$cal_unifrac(samp, tre, weighted)
			N <- nrow(samp)
			sesbeta_matrix <- matrix(nrow = N, ncol = N)
			rownames(sesbeta_matrix) <- colnames(sesbeta_matrix) <- rownames(samp)
			betaobs_vec <- as.vector(betaobs)
			message("Start the randomization ...")
			if(null.model %in% c("taxa.labels", "phylogeny.pool")){
				if(null.model == "taxa.labels"){
					betauni_rand <- replicate(runs, private$cal_unifrac(samp, picante::tipShuffle(tre), weighted))
				}else{
					betauni_rand <- replicate(runs, private$cal_unifrac(picante::randomizeMatrix(samp, null.model = "richness"), tre, weighted))
				}
			}else{
				if(null.model == "sample.pool"){
					betauni_rand <- replicate(runs, private$cal_unifrac(picante::randomizeMatrix(samp, null.model = "richness", iterations = iterations), tre, weighted))
				}else{
					betauni_rand <- replicate(runs, private$cal_unifrac(picante::randomizeMatrix(samp, null.model = null.model, iterations = iterations), tre, weighted))
				}
			}
			message("Calculate the deviation ...")
			betauni_rand_mean <- apply(X = betauni_rand, MARGIN = 1, FUN = mean, na.rm = TRUE)
			betauni_rand_sd <- apply(X = betauni_rand, MARGIN = 1, FUN = sd, na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - betauni_rand_mean)/betauni_rand_sd
			for (i in 1:(N - 1)) {
				x <- 1
				for (j in (i + 1):N) {
					sesbeta_matrix[j, i] <- beta_obs_z[x]
					x <- x + 1
				}
				beta_obs_z <- beta_obs_z[-c(1:(N - i))]
			}
			as.matrix(as.dist(sesbeta_matrix))
		},
		cal_unifrac = function(eco_table, phylo_tree, weighted){
			unifrac1 <- GUniFrac::GUniFrac(eco_table, phylo_tree, alpha = c(0, 0.5, 1))
			unifrac2 <- unifrac1$unifracs
			if(weighted){
				unifrac2[,, "d_1"]
			}else{
				unifrac2[,, "d_UW"]
			}
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
