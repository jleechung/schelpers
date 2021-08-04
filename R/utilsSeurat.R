# ExploreClusterResolutions 
# CombinedPlot
# CondByIdent
# FindMarkersBaseline


#' Visualise clustering parameters
#' @param seu seurat object
#' @param out.dir out directory
#' @param umap.pt umap point size
#' @param tsne.pt tsne point size
#' @param dotplot.cex dotplot axis text size
#' @param n.markers number of markers
#' @param comb.width combined plot width
#' @param comb.height combined plot height 
#' @param umap.width umap plot width
#' @param umap.height umap plot height
#' @param tsne.width tsne plot width
#' @param tsne.height tsne plot height
#' @importFrom Seurat Idents DimPlot FindAllMarkers DotPlot
#' @importFrom ggplot2 ggtitle theme scale_x_discrete
#' @importFrom gridExtra grid.arrange
#' @importFrom dplyr group_by top_n `%>%`
#' @return list of umaps, tsnes and markers, along with plots in out.dir
#' @export
ExploreClusterResolutions <- function(seu, 
                                      out.dir = 'plots',
                                      umap.pt = 2,
                                      tsne.pt = 2,
                                      dotplot.cex = 10,
                                      n.markers = 10,
                                      comb.width = 20,
                                      comb.height = 16,
                                      umap.width = 15,
                                      umap.height = 15,
                                      tsne.width = 15,
                                      tsne.height = 15) {
    
    if (!dir.exists(out.dir)) {
        message(paste0(out.dir, ' does not exist. Creating ', out.dir))
    }
    
    mdNames <- names(seu@meta.data)
    clustNames <- mdNames[grep('res', mdNames)]
    if (length(clustNames) == 0) stop('No clustering output')
    nclust <- length(clustNames)
    
    UMAPS <- vector('list', nclust)
    TSNES <- vector('list', nclust)
    MARKERS <- vector('list', nclust)
    
    for (i in seq_len(nclust)) {
        
        currClustName <- clustNames[i]
        message(paste0('Processing ', currClustName))
        currClust <- seu@meta.data[[clustNames[i]]]
        
        Idents(seu) <- currClust
        
        ## Dimplots
        umap <- DimPlot(seu, reduction = 'umap', cols = scpalette(), pt.size = umap.pt) +
            ggtitle(clustNames[i])
        tsne <- DimPlot(seu, reduction = 'tsne', cols = scpalette(), pt.size = tsne.pt) +
            ggtitle(clustNames[i])
        
        UMAPS[[i]] <- umap
        TSNES[[i]] <- tsne
        
        markers <- FindAllMarkers(seu)
        MARKERS[[i]] <- markers
        ## Get top n for each cluster
        features <- markers %>% group_by(cluster) %>% top_n(n = -n.markers, 
                                                            wt = p_val_adj)
        markerGenes <- unique(features$gene)
        geneSym <- ifelse(test = !grepl('NA', markerGenes), 
                          yes = sub('.*?-', '', markerGenes),
                          no = sub('-.*', '', markerGenes))
        
        ## Dotplot
        dot <- DotPlot(seu, features = unique(features$gene), 
                       cols = 'RdBu', cluster.idents = T) + 
            theme(axis.text.x = element_text(angle = 45, size = dotplot.cex, 
                                             vjust = 0.8, hjust = 0.8)) + 
            scale_x_discrete(labels= geneSym)
        
        layout <- rbind(c(1,2),
                        c(3,3))
        pdf(paste0(out.dir, '/combined_', clustNames[i] ,'.pdf'), 
            width = comb.width, height = comb.height)
        grid.arrange(umap, tsne, dot, layout_matrix = layout)
        dev.off()
    }
    
    pdf(paste0(out.dir, '/umaps.pdf'), height = umap.height, width = umap.width)
    do.call('grid.arrange', UMAPS)
    dev.off()
    
    pdf(paste0(out.dir, '/tsnes.pdf'), height = tsne.height, width = tsne.width)
    do.call('grid.arrange', TSNES)
    dev.off()
    
    names(UMAPS) <- clustNames
    names(TSNES) <- clustNames
    names(MARKERS) <- clustNames
    
    return(list(umaps = UMAPS, tsnes = TSNES, markers = MARKERS))
    
}


#' Combined plots
#' @param seu Seurat object
#' @param features FindAllMarkers output
#' @param out.dir out directory
#' @param plot.name file name for plot
#' @param umap.pt umap point size
#' @param tsne.pt tsne point size
#' @param dotplot.cex dotplot axis test size
#' @param comb.height plot height
#' @param comb.width plot width
#' @importFrom Seurat DimPlot DotPlot
#' @importFrom ggplot2 theme scale_x_discrete
#' @importFrom gridExtra grid.arrange
#' @return plots in out.dir
#' @export
CombinedPlot <- function(seu, 
                         features,
                         out.dir = 'plots',
                         plot.name = 'combined',
                         umap.pt = 2,
                         tsne.pt = 2,
                         dotplot.cex = 10,
                         comb.height = 16, 
                         comb.width = 20) {
    
    pdf(paste0(out.dir, '/', plot.name, '.pdf'), 
        height = comb.height, width = comb.width)
    
    umap <- DimPlot(seu, reduction = 'umap', cols = scpalette(), 
                    label = TRUE, repel = TRUE, pt.size = umap.pt)
    tsne <- DimPlot(seu, reduction = 'tsne', cols = scpalette(), 
                    label = TRUE, repel = TRUE, pt.size = tsne.pt) 
    dot <- DotPlot(seu, features = unique(features$gene), cols = 'RdBu', 
                   cluster.idents = T) + 
        theme(axis.text.x = element_text(angle = 45, size = dotplot.cex, 
                                         vjust = 0.8, hjust = 0.8)) + 
        scale_x_discrete(labels = features$sym)
    layout <- rbind(c(1,2),
                    c(3,3))
    grid.arrange(umap, tsne, dot, layout_matrix = layout)
    dev.off()
    
}

#' Compare conditions for each identity class
#' @param seu Seurat Object
#' @param condition (metadata column which contains levels of conditions
#' @param xlsx.out.dir marker table dir
#' @param xlsx.out.name marker table name
#' @param plot.out.dir plot dir
#' @param plot.name plot name
#' @param plot.height plot height
#' @param plot.width plot width 
#' @param dotplot.cex dotplot axis text size
#' @importFrom openxlsx write.xlsx
#' @importFrom Seurat Idents FindAllMarkers DotPlot
#' @importFrom data.table data.table
#' @importFrom ggplot2 theme element_text scale_x_discrete ggtitle
#' @importFrom gridExtra grid.arrange
#' @return plots and table of markers
#' @export
CondByIdent <- function(seu, 
                        condition = 'SHAM.MI',
                        xlsx.out.dir = 'out',
                        xlsx.out.name = 'condition-by-ident',
                        plot.out.dir = 'plots',
                        plot.name = 'condition-by-ident',
                        plot.height = 12,
                        plot.width = 16,
                        dotplot.cex = 10) {
    
    if (!(condition %in% names(seu@meta.data))) {
        stop('Name a valid condition to split idents by. Should be in metadata')
    }
    
    idents <- unique(Idents(seu))
    out <- DOT <- vector('list', length = length(idents))
    names(out) <- idents
    
    for (i in 1:length(idents)) {
        
        message(paste0('Finding markers for ', idents[i]))
        seuCT <- subset(seu, idents = idents[i])
        Idents(seuCT) <- seuCT@meta.data[[condition]]
        
        if (any(tabulate(Idents(seuCT)) < 2)) {
            message(paste0('Skipping ', idents[i]))
            next
        }
        
        markersCT <- FindAllMarkers(seuCT)
        markersCT <- data.table(markersCT)
        markersCT$sym = ifelse(test = !grepl('NA', markersCT$gene), 
                               yes = sub('.*?-', '', markersCT$gene), 
                               no = sub('-.*', '', markersCT$gene))
        out[[i]] <- markersCT
        p_val_adj <- NULL
        markersCT <- markersCT[p_val_adj < 0.5]
        
        if (nrow(markersCT) > 1) {
            message(paste0('Generating plot for ', idents[i]))
            DOT[[i]] <- DotPlot(seuCT, 
                                features = unique(markersCT$gene), 
                                cols = 'RdBu') + 
                theme(axis.text.x = element_text(angle = 45, 
                                                 size = dotplot.cex, 
                                                 vjust = 0.8, hjust = 0.8)) +
                scale_x_discrete(labels = markersCT$sym) +
                ggtitle(idents[i])
        }
    }
    
    keep <- unlist(lapply(DOT, is.null))
    pdf(paste0(plot.out.dir, "/", plot.name, ".pdf"), 
        height = plot.height, width = plot.width)
    do.call(grid.arrange, DOT[!keep])
    dev.off()
    
    xlsx.out <- paste0(xlsx.out.dir, '/', xlsx.out.name, '.xlsx')
    write.xlsx(out, file = xlsx.out)
    
    names(out) <- names(DOT) <- paste0('cluster ', idents)
    return(list(markers = out, plots = DOT))
}


#' Find Markers for all clusters against a baseline cluster
#' @param seu Seurat object
#' @param baseline cluster to compare all clusters against
#' @param only.pos return up-regulated markers against the baseline only
#' @param n.markers number of markers to draw in boxplot
#' @param dotplot.cex dotplot axis text size
#' @importFrom data.table data.table
#' @importFrom Seurat DotPlot FindMarkers
#' @importFrom ggplot2 theme element_text scale_x_discrete
#' @return plots and table of markers
#' @export
FindMarkersBaseline <- function(seu, baseline, only.pos = TRUE, 
                                n.markers = 30, dotplot.cex = 10) {
    
    if (!(baseline %in% Idents(seu))) stop('Specify a valid baseline 
                                           cluster to compare to')
    
    others <- setdiff(as.numeric(unique(Idents(seu))), baseline)
    others <- sort(others)
    OUT <- DOT <- vector('list', length = length(others))
    
    for (i in 1:length(others)) {
        
        message(paste0('Processing cluster ', others[i]))
        markers <- FindMarkers(seu, ident.1 = others[i], ident.2 = baseline)
        
        if (nrow(markers) == 0) next 
        
        features <- data.table(markers, keep.rownames = 'gene')
        features <- features[p_val_adj < 0.05]
        features <- features[order(avg_log2FC, decreasing = TRUE)]
        features$sym <- ifelse(test = !grepl("NA", features$gene), 
                               yes = sub(".*?-", "", features$gene), 
                               no = sub("-.*", "", features$gene))
        
        if (only.pos) features <- features[avg_log2FC > 0]
        
        OUT[[i]] <- features
        features <- features[1:min(n.markers, nrow(markers))]
        
        if (nrow(features) == 0) next
        
        message(paste0('Generating plot for cluster ', others[i]))
        markerGenes <- unique(features$gene)
        DOT[[i]] <- DotPlot(seu, 
                            features = markerGenes, 
                            idents = c(others[i], baseline),
                            cols = "RdBu", cluster.idents = T) + 
            theme(axis.text.x = element_text(angle = 45,  size = dotplot.cex, 
                                             vjust = 0.8, hjust = 0.8)) + 
            scale_x_discrete(labels = features$sym[match(features$gene, 
                                                         markerGenes)])
        
    }
    
    
    names(OUT) <- names(DOT) <- paste0('cluster ', others)
    return(list(markers = OUT, plots = DOT))
}
