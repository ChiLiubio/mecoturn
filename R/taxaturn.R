#' @title Analyze the 'turnover' of taxa.
#' 
#' @description
#' Analyze the 'turnover' of taxa along a defined gradient.
#'
#' @export
taxaturn <- R6::R6Class(classname = "taxaturn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param taxa_level default "Phylum"; taxonomic rank name, such as "Genus".
		#' 	  If the provided taxonomic name is not a colname in tax_table of input dataset, 
		#' 	  the function will use the features in input \code{microtable$otu_table} automatically.
		#' @param group sample group used for the selection; a colname of input \code{microtable$sample_table}.
		#' @param ordered_group a vector representing the elements of \code{group} parameter.
		#' @param by_group default NULL; a colname of sample_table of input dataset used to for the analysis for each element in the \code{by_group} to get the change result.
		#' 	  If \code{by_group} is not the smallest unit and has repetitions in the elements, mean and sd will be calculated.
		#' @param by_group_final default NULL; NULL or other colname of sample_table of input dataset; 
		#' 	 NULL represents the output is the consistent change across all the elements in \code{by_group};
		#'   a colname of sample_table of input dataset means the consistent change is obtained by those groups instead of all the elements in \code{by_group};
		#'   Note that the by_group_final can be same with by_group, in which the final change is the result of each element in \code{by_group};'
		#'   This parameter is applied only when \code{by_group} is not \code{NULL}.
		#' @return \code{res_change} and \code{res_abund}.
		#' @examples
		#' \donttest{
		#' data(wheat_16S)
		#' t1 <- taxaturn$new(wheat_16S, taxa_level = "Phylum", group = "Type", ordered_group = c("S", "RS", "R"))
		#' t1 <- taxaturn$new(wheat_16S, taxa_level = "Phylum", group = "Type", ordered_group = c("S", "RS", "R"), by_group = "Plant_ID")
		#' }
		initialize = function(dataset, taxa_level = "Phylum", group, ordered_group, by_group = NULL, by_group_final = NULL){
			tmp_dataset <- clone(dataset)
			# first check groups
			if(!group %in% colnames(tmp_dataset$sample_table)){
				stop("Provided group parameter is not the colname of sample_table!")
			}
			if(!all(ordered_group %in% tmp_dataset$sample_table[, group])){
				stop("Part of elements in ordered_group are not found in sample_table[, ", group, "]!")
			}
			if(!is.null(by_group)){
				if(!by_group %in% colnames(tmp_dataset$sample_table)){
					stop("Provided by_group parameter is not the colname of sample_table!")
				}
				if(!is.null(by_group_final)){
					if(!by_group_final %in% colnames(tmp_dataset$sample_table)){
						stop("Provided by_group_final parameter is not the colname of sample_table!")
					}
				}
			}
			# generate abudance table
			if(is.null(tmp_dataset$taxa_abund)){
				message("No taxa_abund found. First calculate it with cal_abund function ...")
				suppressMessages(tmp_dataset$cal_abund())
			}
			# make sure taxa_level can be extracted from taxa_abund
			if(! taxa_level %in% names(tmp_dataset$taxa_abund)){
				# recalculate taxa_abund with rownames as features in otu_table
				message("Provided taxa_level: ", taxa_level, " not in tax_table of dataset; use features in otu_table ...")
				tmp_dataset$add_rownames2taxonomy(use_name = taxa_level)
				suppressMessages(tmp_dataset$cal_abund(rel = TRUE))
			}
			abund_table <- tmp_dataset$taxa_abund[[taxa_level]]
			abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]}

			# all the abudance
			res_abund <- private$abund_change(abund_table = abund_table, sampleinfo = tmp_dataset$sample_table, group = group, by_group = by_group)
			res_change <- data.frame(Taxa = unique(res_abund$Taxa))
			res_change$Direction <- paste0(ordered_group, collapse = " -> ")
			if(is.null(by_group)){
				res_change$Change <- private$chang_func(abund_table = res_abund, group = group, ordered_group = ordered_group)
			}else{
				if(is.null(by_group_final)){
					tmp <- data.frame(Taxa = unique(res_abund$Taxa))
					for(i in unique(res_abund[, by_group])){
						res_abund_bygroup <- res_abund[res_abund[, by_group] == i, ]
						tmp[, paste0(i, " | Change")] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
					}
					res_change$Change <- apply(tmp[, -1, drop = FALSE], 1, function(x){ifelse(length(unique(x)) == 1, unique(x), "")})
				}else{
					if(by_group_final == by_group){
						for(i in unique(res_abund[, by_group])){
							res_abund_bygroup <- res_abund[res_abund[, by_group] == i, ]
							res_change[, paste0(i, " | Change")] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
						}
					}else{
						map_table <- unique(tmp_dataset$sample_table[, c(by_group, by_group_final)])
						rownames(map_table) <- map_table[, by_group]
						res_abund[, by_group_final] <- map_table[res_abund[, by_group], by_group_final]
						for(j in unique(res_abund[, by_group_final])){
							tmp <- data.frame(Taxa = unique(res_abund$Taxa))
							for(i in unique(res_abund[res_abund[, by_group_final] == j, by_group])){
								res_abund_bygroup <- res_abund[res_abund[, by_group] == i, ]
								tmp[, i] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
							}
							res_change[, paste0(j, " | Change")] <- apply(tmp[, -1, drop = FALSE], 1, function(x){ifelse(length(unique(x)) == 1, unique(x), "")})
						}
					}
				}
			}
			
			self$res_abund <- res_abund
			message('Taxa abundance data is stored in object$res_abund ...')
			self$res_change <- res_change
			message('Abundance change data is stored in object$res_change ...')
			self$tmp_dataset <- tmp_dataset
			self$group <- group
			self$ordered_group <- ordered_group
			self$by_group <- by_group
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Differential test of taxonomic abudance across groups
		#'
		#' @param method default "wilcox"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'   }
		#' @param by_group default NULL; a column of sample_table used to perform the differential test 
		#'   among groups (\code{group} parameter) for each group (\code{by_group} parameter). 
		#'   So \code{by_group} has a larger scale than \code{group} parameter in the previous function.
		#' @param by_ID default NULL; a column of sample_table used to perform paired t test or paired wilcox test for the paired data,
		#'   such as the data of plant compartments for different plant species (ID). 
		#'   So \code{by_ID} in sample_table should be the smallest unit of sample collection without any repetition in it.
		#' @param ... parameters (except measure) passed to \code{plot_alpha} function of \code{\link{trans_alpha}} class.
		#' @return \code{res_change}, updated in the object.
		#' @examples
		#' \donttest{
		#' t1$cal_diff()
		#' }
		cal_diff = function(method = c("wilcox", "t.test")[1], by_group = NULL, by_ID = NULL, ...){
			ordered_group <- self$ordered_group
			group <- self$group
			res_change <- self$res_change
			taxa_level <- self$taxa_level
			tmp_dataset <- self$tmp_dataset
			res_table <- data.frame()

			method <- match.arg(method, c("wilcox", "t.test"))
			
			for(i in seq_len(length(ordered_group) - 1)){
				use_ordered_group <- ordered_group[i:(i + 1)]
				tmp <- clone(tmp_dataset)
				tmp$sample_table %<>% base::subset(.[, group] %in% use_ordered_group)
				tmp$tidy_dataset()
				t1 <- suppressMessages(trans_diff$new(dataset = tmp, method = method, group = group, by_group = by_group, by_ID = by_ID, taxa_level = taxa_level, ...))
				t2 <- t1$res_diff
				if(is.null(by_group)){
					rownames(t2) <- t2$Taxa
					sig_value <- t2[res_change$Taxa, "Significance"]
					sig_value[is.na(sig_value)] <- "ns"
					col_name <- paste0(paste0(use_ordered_group, collapse = " -> "), " | Significance")
					res_change[, col_name] <- sig_value
				}else{
					all_bygroups <- tmp_dataset$sample_table %>% dropallfactors %>% .[, by_group] %>% unique
					for(j in all_bygroups){
						t2_tmp <- t2[t2$by_group == j, ]
						rownames(t2_tmp) <- t2_tmp$Taxa
						sig_value <- t2_tmp[res_change$Taxa, "Significance"]
						sig_value[is.na(sig_value)] <- "ns"
						col_name <- paste0(paste0(use_ordered_group, collapse = " -> "), " | ", j, " | Significance")
						res_change[, col_name] <- sig_value
					}
				}
				res_table %<>% rbind(., t1$res_diff)
			}
			res_change$diff_method <- unique(res_table$Test_method)
			self$res_diff_raw <- res_table
			message('Raw differential test results are stored in object$res_diff_raw ...')
			self$res_change <- res_change
			message('Differential test results have been added in object$res_change ...')
		}
	),
	private = list(
		# This function is used for calculating ses.unifrac
		abund_change = function(abund_table, sampleinfo, group, by_group){
			res_abund <- reshape2::melt(tibble::rownames_to_column(abund_table, "Taxa"), id.vars = "Taxa") %>%
				`colnames<-`(c("Taxa", "Sample", "Abund"))
			res_abund <- suppressWarnings(dplyr::left_join(res_abund, tibble::rownames_to_column(sampleinfo), by = c("Sample" = "rowname"))) %>%
				microeco:::summarySE_inter(., measurevar = "Abund", groupvars = c("Taxa", group, by_group))
			res_abund
		},
		chang_func = function(abund_table, group, ordered_group){
			lapply(unique(abund_table$Taxa), function(x){
				abund_sub <- abund_table[abund_table$Taxa == x, ]
				abund_sub <- abund_sub[match(abund_sub[, group], ordered_group), ]
				lagged_abund <- diff(abund_sub$Mean, lag = 1)
				if(all(lagged_abund < 0)){
					"Decrease"
				}else{
					if(all(lagged_abund > 0)){
						"Increase"
					}else{
						""
					}
				}
			}) %>% unlist
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)

