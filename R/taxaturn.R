#' @title Analyze the 'turnover' of taxa.
#' 
#' @description
#' Analyze the 'turnover' of taxa along a defined gradient.
#' The workflow of taxaturn class includes the taxonomic abundance calculation, abundance transformation, abundance change summary, statistical analysis and visualization.
#'
#' @export
taxaturn <- R6::R6Class(classname = "taxaturn",
	public = list(
		#' @param dataset the object of \code{\link{microtable}} class.
		#' @param taxa_level default "Phylum"; taxonomic rank name, such as "Genus". An integer is also acceptable.
		#' 	  If the provided taxa_level is not found in \code{taxa_abund} list,
		#' 	  the function will invoke the \code{cal_abund} function to obtain the relative abudance automatically.
		#' @param group sample group used for the selection; a colname of input \code{microtable$sample_table}.
		#' @param ordered_group a vector representing the ordered elements of \code{group} parameter.
		#' @param by_ID default NULL; a column of sample_table used to obtain the consistent change along provided elements.
		#'   So by_ID can be ID (unique repetition) or even group (with repetitions). 
		#'   If it denotes unique ID, consistent change can be performed across each ID.
		#'   It is also especially useful for the paired wilcox test (or paired t test) in the following analysis.
		#'   If it does not represent unique ID, the mean of each group will be calculated, and consistent change across groups will be obtained.
		#' @param by_group default NULL; NULL or other colname of sample_table of input dataset used to show the result for different groups; 
		#' 	 NULL represents the output is the default consistent change across all the elements in \code{by_ID};
		#'   a colname of sample_table of input dataset means the consistent change is obtained for each group instead of all the elements in \code{by_group};
		#'   Note that the by_group can be same with by_ID, in which the final change is the result of each element in \code{by_group}.
		#'   So generally \code{by_group} has a larger scale than \code{by_ID} parameter in terms of the sample numbers in each element.
		#' @param filter_thres default 0; the mean abundance threshold used to filter features with low abudance.
		#' @return \code{res_change} and \code{res_abund}.
		#' @examples
		#' data(wheat_16S)
		#' t1 <- taxaturn$new(wheat_16S, taxa_level = "Phylum", group = "Type", 
		#'	 ordered_group = c("S", "RS", "R"), by_ID = "Plant_ID")
		initialize = function(dataset, taxa_level = "Phylum", group, ordered_group, by_ID = NULL, by_group = NULL, filter_thres = 0){
			if(!inherits(dataset, "microtable")){
				stop("Input dataset must be microtable class!")
			}
			tmp_dataset <- clone(dataset)
			# first check groups
			if(!group %in% colnames(tmp_dataset$sample_table)){
				stop("Provided group parameter is not the colname of sample_table!")
			}
			if(!all(ordered_group %in% tmp_dataset$sample_table[, group])){
				stop("Part of elements in ordered_group are not found in sample_table[, ", group, "]!")
			}
			if(!is.null(by_ID)){
				if(!by_ID %in% colnames(tmp_dataset$sample_table)){
					stop("Provided by_ID parameter is not the colname of sample_table!")
				}
			}
			if(!is.null(by_group)){
				if(!by_group %in% colnames(tmp_dataset$sample_table)){
					stop("Provided by_group parameter is not the colname of sample_table!")
				}
			}
			if(is.null(taxa_level)){
				stop("The input taxa_level should not be NULL!")
			}else{
				if(length(taxa_level) > 1){
					stop("The input taxa_level must only have one!")
				}else{
					if(is.numeric(taxa_level)){
						if(taxa_level > ncol(tmp_dataset$tax_table)){
							stop("The input taxa_level is larger than the column number of tax_table!")
						}else{
							taxa_level <- colnames(tmp_dataset$tax_table)[taxa_level]
							message("Rename taxa_level: ", taxa_level, " ...")
						}
					}else{
						if(!taxa_level %in% colnames(tmp_dataset$tax_table)){
							if(!taxa_level %in% names(tmp_dataset$taxa_abund)){
								stop("The input taxa_level is not found in the colnames of tax_table and names of taxa_abund list!")
							}
						}
					}
				}
			}
			# generate abudance table
			if(is.null(tmp_dataset$taxa_abund)){
				message("No taxa_abund found. First calculate it with cal_abund function ...")
				suppressMessages(tmp_dataset$cal_abund(select_cols = taxa_level))
			}else{
				if(!taxa_level %in% names(tmp_dataset$taxa_abund)){
					suppressMessages(tmp_dataset$cal_abund(select_cols = taxa_level))
				}
			}
			abund_table <- tmp_dataset$taxa_abund[[taxa_level]]
			if(filter_thres > 0){
				abund_table %<>% .[apply(., 1, mean) > filter_thres, ]
				if(nrow(abund_table) == 0){
					stop("Please check the filter_thres parameter! No feature remained after filtering!")
				}
			}
			abund_table %<>% {.[!grepl("__$|uncultured$|Incertae..edis$|_sp$", rownames(.), ignore.case = TRUE), ]}

			# transform the abundance along by_ID or across by_group
			if(is.null(by_ID)){
				res_abund <- private$abund_change(abund_table = abund_table, sampleinfo = tmp_dataset$sample_table, group = group, by_group = by_group)
			}else{
				res_abund <- private$abund_change(abund_table = abund_table, sampleinfo = tmp_dataset$sample_table, group = group, by_group = by_ID)
			}
			res_change <- data.frame(Taxa = unique(res_abund$Taxa))
			res_change$Direction <- paste0(ordered_group, collapse = " -> ")
			
			if(is.null(by_ID)){
				if(is.null(by_group)){
					res_change$Change <- private$chang_func(abund_table = res_abund, group = group, ordered_group = ordered_group)
				}else{
					for(i in unique(res_abund[, by_group])){
						res_abund_bygroup <- res_abund[res_abund[, by_group] == i, ]
						res_change[, paste0(i, " | Change")] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
					}
				}
			}else{
				if(is.null(by_group)){
					tmp <- data.frame(Taxa = unique(res_abund$Taxa))
					for(i in unique(res_abund[, by_ID])){
						res_abund_bygroup <- res_abund[res_abund[, by_ID] == i, ]
						tmp[, paste0(i, " | Change")] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
					}
					res_change$Change <- apply(tmp[, -1, drop = FALSE], 1, function(x){ifelse(length(unique(x)) == 1, unique(x), "")})
				}else{
					if(by_group == by_ID){
						for(i in unique(res_abund[, by_ID])){
							res_abund_bygroup <- res_abund[res_abund[, by_ID] == i, ]
							res_change[, paste0(i, " | Change")] <- private$chang_func(abund_table = res_abund_bygroup, group = group, ordered_group = ordered_group)
						}
					}else{
						map_table <- unique(tmp_dataset$sample_table[, c(by_ID, by_group)])
						rownames(map_table) <- map_table[, by_ID]
						res_abund[, by_group] <- map_table[res_abund[, by_ID], by_group]
						for(j in unique(res_abund[, by_group])){
							tmp <- data.frame(Taxa = unique(res_abund$Taxa))
							for(i in unique(res_abund[res_abund[, by_group] == j, by_ID])){
								res_abund_bygroup <- res_abund[res_abund[, by_ID] == i, ]
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
			self$by_ID <- by_ID
			self$by_group <- by_group
			self$taxa_level <- taxa_level
		},
		#' @description
		#' Differential test of taxonomic abundance across groups
		#'
		#' @param method default "wilcox"; see the following available options:
		#'   \describe{
		#'     \item{\strong{'wilcox'}}{Wilcoxon Rank Sum and Signed Rank Tests for all paired groups }
		#'     \item{\strong{'t.test'}}{Student's t-Test for all paired groups}
		#'     \item{\strong{'betareg'}}{Beta Regression based on the \code{betareg} package}
		#'     \item{\strong{'lme'}}{lme: Linear Mixed Effect Model based on the \code{lmerTest} package}
		#'     \item{\strong{'anova'}}{one-way or multi-way anova}
		#'   }
		#' @param ... parameters passed to \code{trans_diff$new}.
		#' @return \code{res_change} or \code{res_diff}.
		#' @examples
		#' t1$cal_diff(method = "wilcox")
		cal_diff = function(method = c("wilcox", "t.test", "betareg", "lme", "anova")[1], ...){
			ordered_group <- self$ordered_group
			group <- self$group
			by_ID <- self$by_ID
			by_group <- self$by_group
			res_change <- self$res_change
			taxa_level <- self$taxa_level
			tmp_dataset <- self$tmp_dataset
			method <- match.arg(method, c("wilcox", "t.test", "betareg", "lme", "anova"))
			
			if(method %in% c("wilcox", "t.test")){
				res_table <- data.frame()
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
				self$res_diff_raw <- res_table
				message('Raw differential test results are stored in object$res_diff_raw ...')
				self$res_change <- res_change
				message('Differential test results have been added into object$res_change ...')
			}
			if(method == "betareg"){
				tmp <- clone(tmp_dataset)
				message("Convert ", group, " to numeric class for the beta regression ...")
				tmp$sample_table[, group] %<>% factor(levels = ordered_group) %>% as.numeric
				if(is.null(by_group)){
					t1 <- suppressMessages(trans_diff$new(dataset = tmp, method = method, taxa_level = taxa_level, ...))
					res_table <- t1$res_diff
				}else{
					res_table <- data.frame()
					all_bygroups <- tmp_dataset$sample_table %>% dropallfactors %>% .[, by_group] %>% unique
					for(j in all_bygroups){
						tmp2 <- clone(tmp)
						tmp2$sample_table %<>% .[.[, by_group] %in% j, ]
						tmp2$tidy_dataset()
						t1 <- suppressMessages(trans_diff$new(dataset = tmp2, method = method, taxa_level = taxa_level, ...))
						tmp_output <- data.frame(by_group = j, t1$res_diff)
						res_table %<>% rbind(., tmp_output)
					}
				}
				self$res_diff <- res_table
				message('The results are stored in object$res_diff ...')
			}
			if(method %in% c("lme", "anova")){
				tmp <- clone(tmp_dataset)
				t1 <- suppressMessages(trans_diff$new(dataset = tmp, method = method, taxa_level = taxa_level, ...))
				self$res_diff <- t1$res_diff
				message('The results are stored in object$res_diff ...')
			}
		},
		#' @description
		#' Plot the line chart.
		#'
		#' @param select_taxon default NULL; a taxon name.
		#'   Note that if \code{delete_prefix} is TRUE, the provided select_taxon should be taxa names without long prefix (those before |);
		#'   if \code{delete_prefix} is FALSE, the select_taxon should be full names same with those in the \code{res_abund} of the object.
		#' @param color_values default \code{RColorBrewer::brewer.pal}(8, "Dark2"); colors palette for the plotting.
		#' @param delete_prefix default TRUE; whether delete the prefix in the taxa names.
		#' @param plot_type default c("point", "line", "errorbar", "smooth")[1:3]; a vector of visualization types. Multiple elements are available. 
		#'   'smooth' denotes the fitting with \code{geom_smooth} function of ggplot2 package.
		#' @param errorbar_SE default TRUE; TRUE: plot the errorbar with mean ± se; FALSE: plot the errorbar with mean ± sd.
		#' @param rect_fill default TRUE; Whether fill color in each rectangular area.
		#' @param rect_color default c("grey70", "grey90"); the colors used to fill different rectangular area.
		#' @param rect_alpha default 0.2; the fill color transparency in rectangular area.
		#' @param position default position_dodge(0.1); Position adjustment for the points and lines, either as a string (such as "identity"), 
		#'   or the result of a call to a position adjustment function.
		#' @param errorbar_size default 1; errorbar size.
		#' @param errorbar_width default 0.1; errorbar width.
		#' @param point_size default 3; point size for taxa.
		#' @param point_alpha default 0.8; point transparency.
		#' @param line_size default 0.8; line size.
		#' @param line_alpha default 0.8; line transparency.
		#' @param line_type default 1; an integer; line type.
		#' @param ... parameters passed to \code{geom_smooth} when 'smooth' is in plot_type parameter.
		#' @return ggplot2 plot. 
		#' @examples
		#' t1$plot()
		plot = function(
			select_taxon = NULL,
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			delete_prefix = TRUE,
			plot_type = c("point", "line", "errorbar", "smooth")[1:3],
			errorbar_SE = TRUE,
			rect_fill = TRUE,
			rect_color = c("grey70", "grey90"),
			rect_alpha = 0.2,
			position = position_dodge(0.1),
			errorbar_size = 1,
			errorbar_width = 0.1,
			point_size = 3,
			point_alpha = 0.8,
			line_size = 0.8,
			line_alpha = 0.8,
			line_type = 1,
			...
			){
			by_ID <- self$by_ID
			by_group <- self$by_group
			plot_data <- self$res_abund
			group <- self$group
			ordered_group <- self$ordered_group
			plot_data[, group] %<>% factor(levels = ordered_group)
			if(delete_prefix){
				plot_data$Taxa %<>% gsub(".*\\|", "", .)
			}
			if(is.null(select_taxon)){
				# sort abudance of taxa for the selection
				plot_data$total_abund <- plot_data$N * plot_data$Mean
				taxa_abund <- tapply(plot_data$total_abund, plot_data$Taxa, sum) %>% sort(decreasing = TRUE)
				use_taxa_name <- names(taxa_abund[1])
				message("No select_taxon provided! Select the taxon with the highest abudance: ", use_taxa_name, " ...")
			}else{
				if(length(select_taxon) > 1){
					message("Provided select_taxon has multiple elements! Select the first one ...")
				}
				use_taxa_name <- select_taxon[1]
				if(! use_taxa_name %in% plot_data$Taxa){
					stop("Please check the input select_taxon parameter! It is not correct!")
				}
			}
			plot_data %<>% .[.$Taxa %in% use_taxa_name, ]
			plot_data$Taxa %<>% factor(levels = use_taxa_name)
			
			if(is.null(by_ID)){
				if(is.null(by_group)){
					p <- ggplot(plot_data, aes(!!as.symbol(group), Mean, group = Taxa))
				}else{
					p <- ggplot(plot_data, aes(!!as.symbol(group), Mean, group = !!as.symbol(by_group), color = !!as.symbol(by_group)))
				}
			}else{
				if(is.null(by_group)){
					p <- ggplot(plot_data, aes(!!as.symbol(group), Mean, group = !!as.symbol(by_ID)))
				}else{
					p <- ggplot(plot_data, aes(!!as.symbol(group), Mean, group = !!as.symbol(by_ID), color = !!as.symbol(by_group)))
				}
			}
			# + facet_grid(Taxa ~ ., drop = TRUE, scale = "free", space = "fixed")
			p <- p + theme_bw()
			if(rect_fill){
				for(i in seq_len(length(ordered_group) - 1)){
					if(i %% 2 == 1){
						p <- p + geom_rect(xmin = i, xmax = i + 1, ymin = -Inf, ymax = Inf, alpha = rect_alpha, fill = rect_color[1], colour = rect_color[1])
					}else{
						p <- p + geom_rect(xmin = i, xmax = i + 1, ymin = -Inf, ymax = Inf, alpha = rect_alpha, fill = rect_color[2], colour = rect_color[2])
					}
				}
			}
			if("errorbar" %in% plot_type){
				if(("SE" %in% colnames(plot_data)) & errorbar_SE){
					p <- p + geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = errorbar_width, position = position, linewidth = errorbar_size)
				}else{
					if(("SD" %in% colnames(plot_data)) & errorbar_SE){
						p <- p + geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = errorbar_width, position = position, linewidth = errorbar_size)
					}
				}
			}
			if("point" %in% plot_type){
				p <- p + geom_point(size = point_size, alpha = point_alpha, position = position)
			}
			if("line" %in% plot_type){
				p <- p + geom_line(linewidth = line_size, alpha = line_alpha, linetype = line_type, position = position)
			}
			if("smooth" %in% plot_type){
				p <- p + geom_smooth(position = position, ...)
			}
			p <- p + theme(strip.background = element_rect(fill = "grey95"), strip.text.y = element_text(angle = 360)) +
				ylab("Relative abundance") + 
				xlab("") +
				scale_color_manual(values = color_values)

			p
		}
	),
	private = list(
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
						paste0(sapply(lagged_abund, private$get_symbol), collapse = "|")
					}
				}
			}) %>% unlist
		},
		get_symbol = function(x){
			if(!is.numeric(x)){
				stop("The input must be numeric!")
			}
			if(x > 0){
				"+"
			}else{
				if(x < 0){
					"-"
				}else{
					"0"
				}
			}
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)

