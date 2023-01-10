#' @title Analyze the 'turnover' of beta-diversity.
#' 
#' @description
#' Analyze the 'turnover' of beta-diversity.
#'
#' @export
betaturn <- R6::R6Class(classname = "betaturn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param measure default NULL; beta diversity dissimilarity metric; 
		#' 	 must be one of \code{c("bray", "jaccard", "wei_unifrac", "unwei_unifrac", "betaMPD", "betaMNTD", "betaNRI", "betaNTI", "ses_UniFrac", "RCbray")}
		#' 	 or other options in parameter \code{method} of \code{\link{vegan::vegdist}} function.
		#' @param filter_thres default 0; the relative abundance threshold used to filter features.
		#' @param abundance.weighted default TRUE; whether use abundance-weighted method for the phylogenetic metrics.
		#' @param null.model default NULL; one of \code{c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap")},
		#'   in which "taxa.labels" can only used for phylogenetic analysis.
		#' @param runs default 1000; simulation runs of null model.
		#' @param iterations default 1000; iteration number for part null models to perform; see iterations parameter of \code{picante::randomizeMatrix} function.
		#' @param ... parameters passed to \code{cal_betadiv} function of \code{microtable} class when provided measure is not in the current vector;
		#'   parameters passed to \code{cal_betamntd} function of \code{trans_nullmodel} class when \code{measure = "betaMNTD"};
		#'   parameters passed to \code{cal_ses_betamntd} function of \code{trans_nullmodel} class when \code{measure = "betaNTI"}.
		#' @return \code{dataset}, stored in the object. The new dataset has a beta_diversity list and the calculated distance matrix in the list.
		#' @examples
		#' data(wheat_16S)
		#' b1 <- betaturn$new(wheat_16S, measure = "bray")
		initialize = function(dataset, measure = NULL, filter_thres = 0, abundance.weighted = TRUE, null.model = NULL, runs = 1000, iterations = 1000, ...){
			set.seed(123) # set the random seed
			tmp <- clone(dataset)
			if(filter_thres > 0){
				tmp$otu_table %<>% .[apply(., 1, sum)/sum(.) > filter_thres, ]
				tmp$tidy_dataset()
			}
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
			}
			if(!is.null(measure)){
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
						res_dis <- private$ses_unifrac(tmp$otu_table, tmp$phylo_tree, runs = runs, weighted = abundance.weighted, null.model = null.model, iterations = iterations)
					}
					tmp$beta_diversity[[measure]] <- res_dis
				}
			}
			self$measure <- measure
			self$dataset <- tmp
			message('The distance matrix is stored in object$beta_diversity list ...')
		},
		#' @description
		#' Transform sample distances within groups or between groups.
		#'
		#' @param group one colname of sample_table in \code{microtable} object used for group distance convertion.
		#' @param within_group default TRUE; whether transform sample distance within groups, if FALSE, transform sample distance between any two groups.
		#' @param by_group default NULL; one colname of sample_table in \code{microtable} object.
		#'   If provided, convert distances according to the provided by_group parameter. This is especially useful for ordering and filtering values further.
		#'   When \code{within_group = TRUE}, the result of by_group parameter is the format of paired groups.
		#'   When \code{within_group = FALSE}, the result of by_group parameter is the format same with the group information in \code{sample_table}.
		#' @param ordered_group default NULL; a vector representing the ordered elements of \code{group} parameter; only useful when within_group = FALSE.
		#' @param sep default TRUE; a character string to separate the group names after merging them into a new name.
		#' @return \code{res_group_distance} stored in object.
		#' @examples
		#' \donttest{
		#' b1$cal_group_distance(group = "Type", within_group = FALSE, by_group = "Plant_ID")
		#' }
		cal_group_distance = function(group, within_group = TRUE, by_group = NULL, ordered_group = NULL, sep = " vs "){
			dataset <- self$dataset
			measure <- self$measure
			trans_beta_object <- trans_beta$new(dataset = dataset, group = group, measure = measure)
			suppressMessages(trans_beta_object$cal_group_distance(within_group = within_group, by_group = by_group, ordered_group = ordered_group, sep = sep))
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
		#' \donttest{
		#' b1$cal_group_distance_diff(method = "wilcox")
		#' }
		cal_group_distance_diff = function(...){
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
		#' \donttest{
		#' b1$plot_group_distance()
		#' }
		plot_group_distance = function(...){
			self$tmp_trans_beta$plot_group_distance(...)
		}
	),
	private = list(
		# This function is used for calculating ses.unifrac
		ses_unifrac = function(samp, tre, null.model = "taxa.labels", runs = 1000, weighted = TRUE, iterations = 1000) {
			# samp: OTU table, rownames are OTU id; tre: phylonetic tree, .tre file
			library(phyloseq) # download this package from bioconductor
			if (sum(!(rownames(samp) %in% tre$tip.label)) != 0) {
				stop("The OTU table contains unknown OTUs! OTU names\n\t\t\t\tin the OTU table and the tree should match!")
			}
			samp_table <- otu_table(samp, taxa_are_rows=TRUE)
			message("Run the observed UniFrac ...")
			betaobs <- UniFrac(phyloseq(samp_table, tre), weighted)
			N <- dim(samp)[2]
			sesbeta_matrix <- matrix(nrow = N, ncol = N)
			rownames(sesbeta_matrix) <- colnames(sesbeta_matrix) <- colnames(samp)
			betaobs_vec <- as.vector(betaobs)
			message("Start the randomization ...")
			if(null.model %in% c("taxa.labels", "phylogeny.pool")){
				if(null.model == "taxa.labels"){
					betauni_rand <- replicate(runs, UniFrac(phyloseq(samp_table, picante::tipShuffle(tre)), weighted))
				}else{
					betauni_rand <- replicate(runs, UniFrac(phyloseq(picante::randomizeMatrix(samp_table, null.model = "richness"), tre), weighted))
				}
			}else{
				if(null.model == "sample.pool"){
					betauni_rand <- replicate(runs, UniFrac(phyloseq(picante::randomizeMatrix(samp_table, null.model = "richness", iterations = iterations), tre), weighted))
				}else{
					betauni_rand <- replicate(runs, UniFrac(phyloseq(picante::randomizeMatrix(samp_table, null.model = null.model, iterations = iterations), tre), weighted))
				}
			}
			message("Calculate the deviation ...")
			betauni_rand_mean <- apply(X = betauni_rand, MARGIN = 1, FUN = mean, 
									 na.rm = TRUE)
			betauni_rand_sd <- apply(X = betauni_rand, MARGIN = 1, FUN = sd, 
								   na.rm = TRUE)
			beta_obs_z <- (betaobs_vec - betauni_rand_mean)/betauni_rand_sd
			for (i in 1:(N - 1)) {
				x <- 1
				for (j in (i + 1):N) {
					sesbeta_matrix[j,i] <- beta_obs_z[x]
					x <- x+1
				}
				beta_obs_z <- beta_obs_z[-c(1:(N-i))]
			}
			return(as.matrix(as.dist(sesbeta_matrix)))
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
