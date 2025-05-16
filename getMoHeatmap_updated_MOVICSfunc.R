

getMoHeatmap <- function(data             = NULL,
                         is.binary        = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                         row.title        = c("Data1","Data2","Data3","Data4","Data5","Data6"),
                         legend.name      = c("Data1","Data2","Data3","Data4","Data5","Data6"),
                         clust.res        = NULL,
                         clust.dend       = NULL,
                         show.col.dend    = TRUE,
                         show.colnames    = FALSE,
                         show.row.dend    = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                         show.rownames    = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                         clust.dist.row   = c("pearson","pearson","pearson","pearson","pearson","pearson"),
                         clust.method.row = c("ward.D","ward.D","ward.D","ward.D","ward.D","ward.D"),
                         clust.col        = c("black","grey","#FF9F1C","#BDD5EA","#FFA5AB","#011627","#023E8A","#9D4EDD"),
                         color            = rep(list(c("#00FF00", "#000000", "#FF0000")),length(data)),
                         annCol           = NULL,
                         annColors        = NULL,
                         annRow           = NULL,
                         width            = 6,
                         height           = 4,
                         fig.path         = getwd(),
                         fig.name         = "moheatmap") {
  
  ht_opt$message = FALSE
  defaultW <- getOption("warn")
  options(warn = -1)
  
  # check data
  if(is.null(names(data))){
    names(data) <- sprintf("dat%s", 1:length(data))
  }
  
  n_dat <- length(data)
  if(n_dat > 6){
    stop('current verision of MOVICS can support up to 6 datasets.')
  }
  if(n_dat < 2){
    stop('current verision of MOVICS needs at least 2 omics data.')
  }
  
  colvec <- clust.col[1:length(unique(clust.res$clust))]
  names(colvec) <- paste0(#"CS",
                          unique(clust.res$clust))
  
  if(!is.null(annCol) & !is.null(annColors)) {
    
    annCol <- annCol[colnames(data[[1]]), , drop = FALSE]
    annCol$COCA <- paste0("CS",clust.res[colnames(data[[1]]),"clust"])
    annColors[["COCA"]] <- colvec
    
    if(is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust),]
      annCol <- annCol[clust.res$samID, , drop = FALSE]
    }
    
    ha <- ComplexHeatmap::HeatmapAnnotation(df     = annCol,
                                            col    = annColors,
                                            border = FALSE)
  } else {
    annCol <- data.frame("BCC" = paste0(# "CS",
                                         clust.res[colnames(data[[1]]),"clust"]), # Change name here----
                         row.names = colnames(data[[1]]),
                         stringsAsFactors = FALSE)
    annColors <- list("BCC" = colvec) # Change name here----
    
    if(is.null(clust.dend)) {
      clust.res <- clust.res[order(clust.res$clust),]
      annCol <- annCol[clust.res$samID,,drop = FALSE]
    }
    
    ha <- ComplexHeatmap::HeatmapAnnotation(df     = annCol,
                                            col    = annColors,
                                            border = FALSE)
  }
  
  if(!is.null(annRow)) {
    if(!is.list(annRow)) {stop("argument of annRow should be a list!")}
  }
  
  ht <- list()
  for (i in 1:n_dat) {
    
    hcg <- hclust(ClassDiscovery::distanceMatrix(as.matrix(t(data[[i]])), clust.dist.row[i]), clust.method.row[i])
    
    if(is.null(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    } else if (is.na(annRow[[i]][1])) {
      rowlab <- ""
      rowlab.index <- 0
    } else {
      rowlab <- intersect(rownames(data[[i]]),annRow[[i]])
      rowlab.index <- match(rowlab, rownames(data[[i]]))
    }
    
    if(is.null(clust.dend)) {
      data <- lapply(data, function(x) x[,clust.res$samID])
      
      if(!is.binary[i]) {
        print(ha)
        ha@anno_list[["BCC"]]@show_legend <- FALSE # Remove legend of main cluster - need update name of cluster
        
        ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(data[[i]]),
                                            row_title            = row.title[i],
                                            name                 = legend.name[i],
                                            cluster_columns      = FALSE,
                                            cluster_rows         = hcg,
                                            show_column_dend     = FALSE,
                                            show_column_names    = show.colnames,
                                            show_row_dend        = show.row.dend[i],
                                            show_row_names       = show.rownames[i],
                                            col                  = grDevices::colorRampPalette(color[[i]])(64),
                                            top_annotation       = switch((i == 1) + 1, NULL, ha),
                                            width                = grid::unit(width, "cm"),
                                            height               = grid::unit(height, "cm"),
                                            heatmap_legend_param = list(at     = pretty(range(data[[i]])),
                                                                        labels = pretty(range(data[[i]]))),
                                            show_heatmap_legend = FALSE, # Remove L1, Alu, LTR legends if FALSE
                                            right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                                   anno_mark(at         = rowlab.index,
                                                                                                             labels     = rowlab,
                                                                                                             which      = "row",
                                                                                                             lines_gp   = grid::gpar(fontsize = 5),
                                                                                                             link_width = grid::unit(3, "mm"),
                                                                                                             padding    = grid::unit(0.8, "mm"),
                                                                                                             labels_gp  = grid::gpar(fontsize = 7))))
      } else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        
        ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(data[[i]]),
                                            row_title            = row.title[i],
                                            name                 = legend.name[i],
                                            cluster_columns      = FALSE,
                                            cluster_rows         = hcg,
                                            show_column_dend     = FALSE,
                                            show_column_names    = show.colnames,
                                            show_row_dend        = show.row.dend[i],
                                            show_row_names       = show.rownames[i],
                                            col                  = color[[i]],
                                            top_annotation       = switch((i == 1) + 1, NULL, ha),
                                            width                = grid::unit(width, "cm"),
                                            height               = grid::unit(height, "cm"),
                                            heatmap_legend_param = list(at        = c(0, 1),
                                                                        legend_gp = grid::gpar(fill = col_fun(c(0,1))),
                                                                        labels    = c("0", "1")),
                                            right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                                   anno_mark(at         = rowlab.index,
                                                                                                             labels     = rowlab,
                                                                                                             which      = "row",
                                                                                                             lines_gp   = grid::gpar(fontsize = 5),
                                                                                                             link_width = grid::unit(3, "mm"),
                                                                                                             padding    = grid::unit(0.8, "mm"),
                                                                                                             labels_gp  = grid::gpar(fontsize = 7))))
      }
      
    } else {
      if(!is.binary[i]) {
        ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(data[[i]]),
                                            row_title            = row.title[i],
                                            name                 = legend.name[i],
                                            cluster_columns      = clust.dend,
                                            cluster_rows         = hcg,
                                            show_column_dend     = show.col.dend,
                                            show_column_names    = show.colnames,
                                            show_row_dend        = show.row.dend[i],
                                            show_row_names       = show.rownames[i],
                                            col                  = grDevices::colorRampPalette(color[[i]])(64),
                                            top_annotation       = switch((i == 1) + 1, NULL, ha),
                                            width                = grid::unit(width, "cm"),
                                            height               = grid::unit(height, "cm"),
                                            heatmap_legend_param = list(at     = pretty(range(data[[i]])),
                                                                        labels = pretty(range(data[[i]]))),
                                            right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                                   anno_mark(at         = rowlab.index,
                                                                                                             labels     = rowlab,
                                                                                                             which      = "row",
                                                                                                             lines_gp   = grid::gpar(fontsize = 5),
                                                                                                             link_width = grid::unit(3, "mm"),
                                                                                                             padding    = grid::unit(0.8, "mm"),
                                                                                                             labels_gp  = grid::gpar(fontsize = 7))))
      } else {
        col_fun = circlize::colorRamp2(c(0, 1), color[[i]])
        
        ht[[i]] <-  ComplexHeatmap::Heatmap(matrix               = as.matrix(data[[i]]),
                                            row_title            = row.title[i],
                                            name                 = legend.name[i],
                                            cluster_columns      = clust.dend,
                                            cluster_rows         = hcg,
                                            show_column_dend     = show.col.dend,
                                            show_column_names    = show.colnames,
                                            show_row_dend        = show.row.dend[i],
                                            show_row_names       = show.rownames[i],
                                            col                  = color[[i]],
                                            top_annotation       = switch((i == 1) + 1, NULL, ha),
                                            width                = grid::unit(width, "cm"),
                                            height               = grid::unit(height, "cm"),
                                            heatmap_legend_param = list(at        = c(0, 1),
                                                                        legend_gp = grid::gpar(fill = col_fun(c(0,1))),
                                                                        labels    = c("0", "1")),
                                            right_annotation     = ComplexHeatmap::rowAnnotation(link =
                                                                                                   anno_mark(at         = rowlab.index,
                                                                                                             labels     = rowlab,
                                                                                                             which      = "row",
                                                                                                             lines_gp   = grid::gpar(fontsize = 5),
                                                                                                             link_width = grid::unit(3, "mm"),
                                                                                                             padding    = grid::unit(0.8, "mm"),
                                                                                                             labels_gp  = grid::gpar(fontsize = 7))))
      }
    }
  }
  
  if(n_dat == 1){
    ht_list <- ht[[1]]
  }
  if(n_dat == 2){
    ht_list <- ht[[1]] %v% ht[[2]]
  }
  if(n_dat == 3){
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]]
  }
  if(n_dat == 4){
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]]
  }
  if(n_dat == 5){
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% ht[[5]]
  }
  if(n_dat == 6){
    ht_list <- ht[[1]] %v% ht[[2]] %v% ht[[3]] %v% ht[[4]] %v% ht[[5]] %v% ht[[6]]
  }
  
  outFile <- file.path(fig.path,paste0(fig.name,".pdf"))
  if(is.null(annCol)) {
    pdf(outFile, width = width, height = height * n_dat/2)
  } else {
    pdf(outFile, width = width, height = height * n_dat/1.5)
  }
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right") # output to pdf
  invisible(dev.off())
  
  draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "right") # output to screen
  
  options(warn = defaultW)
}
