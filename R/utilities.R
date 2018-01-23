library(magrittr)
#' Plot each row or column of a matrix
#' 
#' @param sample_mat The matrix to plot
#' @param plot_dim Which dimension to plot - must be one of "columns" or "rows"
#' @param add Should the plot be added to an existing. False by default.
#' 
#' 
plot_matrix <- function(target_matrix, plot_dim =c("columns","rows"), add=F, ...){
  plot_dim <- match.arg(plot_dim)
  if (plot_dim == "rows"){
    target_matrix <- t(target_matrix)
  }
  
  iters <- 1:ncol(target_matrix)
  if (!add){
    plot(target_matrix[,1], ylim=range(target_matrix), type="l", ...)  
  }
  sapply(iters, function(j) points(target_matrix[,j], type="l", ...))
}


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Take from:
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



#' Get empirical log-rates dataframe
#' @param prefix Gives the path to the top level of the repository.
#' @return Dataframe containing empirical log rates in tidy format
load_empirical_rates <- function(prefix=""){
  expos <-  readRDS(paste0(prefix, "data/hmd/exposures_hmd.rds"))
  deaths <- readRDS(paste0(prefix, "data/hmd/deaths_hmd.rds"))
  expos %<>% dplyr::select(-OpenInterval) %>% 
    tidyr::gather(Sex, Exposure, -Year, -Age)
  deaths %<>% dplyr::select(-OpenInterval) %>% 
    tidyr::gather(Sex, Deaths, -Year, -Age)
  emp_rates_df <- dplyr::left_join(expos, deaths) %>%
    dplyr::mutate(Rate= Deaths/Exposure, Log_rate = log(Rate)) %>% 
    as_tibble()
  return(emp_rates_df)
}


#' Capitalise words
#' 
#' Taken from help(tolower) examples
#' 
#' @param character vector to be capitialised
#' @param strict Should already in all caps be decapitalised?
capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  return(sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s))))
}

#' Find everything in the given path that is not a directory
#' @param path The path to scanned
#' @return A vector of the names of the files at the given path.
get_non_directories <- function(path="."){
  return(setdiff(list.files(path), 
                 list.dirs(recursive = FALSE, full.names = FALSE)))
  
}

#' Tidy up ONS variant NPP projections
#'
tidy_ons_data <- function(variant){
  var_df <- gather(variant, key=Year_mid, value=qx,-Sex,-Age)
  var_df <- var_df %>% mutate(qx = qx/100000,
                              Year_mid =get_ons_year(Year_mid))
  
  var_df <- var_df %>% mutate( log_qx =log(qx))
  var_df <- var_df %>% as_tibble %>% mutate(Age = ifelse(Age=="Birth", -1,Age))
  var_df <- var_df %>% mutate(Age=as.numeric(Age))
  return(var_df)
}

#' Extract year from ons column label.
#'
get_ons_year<-function(year_string){
  year_dots  <-  gsub("X", "",year_string) 
  year <-  sapply(year_dots, function(x) strsplit(x,"\\.")[[1]][1])
  return(as.numeric(year))
}