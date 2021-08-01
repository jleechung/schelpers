
#' Peak into objects
#' @param x a data.frame or matrix
#' @return preview
#' @export
see <- function(x) {
    x[1:5,1:5]
}

#' Unique length
#' @param x a vector
#' @return unique length
#' @export
ulen <- function(x) {
    length(unique(x))
}

#' Kelly palette
#' @param n number of colors
#' @return colors
#' @export
kellyPal <- function(n=21) {
    
    col <- c(
        rgb(0, 0, 0, maxColorValue = 255, alpha = 200),
        rgb(255, 179, 0, maxColorValue = 255, alpha = 200),
        rgb(128, 62, 117, maxColorValue = 255, alpha = 200),
        rgb(255, 104, 0, maxColorValue = 255, alpha = 200),
        
        rgb(166, 189, 215, maxColorValue = 255, alpha = 200),
        rgb(193, 0, 32, maxColorValue = 255, alpha = 200),
        rgb(206, 162, 98, maxColorValue = 255, alpha = 200),
        rgb(0, 125, 52, maxColorValue = 255, alpha = 200),
        
        rgb(246, 118, 142, maxColorValue = 255, alpha = 200),
        rgb(0, 83, 138, maxColorValue = 255, alpha = 200),
        rgb(129, 112, 102, maxColorValue = 255, alpha = 200),
        rgb(179, 40, 81, maxColorValue = 255, alpha = 200),
        
        rgb(255, 122, 92, maxColorValue = 255, alpha = 200),
        rgb(83, 55, 122, maxColorValue = 255, alpha = 200),
        rgb(255, 142, 0, maxColorValue = 255, alpha = 200),
        rgb(244, 200, 0, maxColorValue = 255, alpha = 200),
        
        rgb(127, 24, 13, maxColorValue = 255, alpha = 200),
        rgb(147, 170, 0, maxColorValue = 255, alpha = 200),
        rgb(89, 51, 21, maxColorValue = 255, alpha = 200),
        rgb(241, 58, 19, maxColorValue = 255, alpha = 200),
        rgb(35, 44, 22, maxColorValue = 255, alpha = 200)
        
    )
    
    if (n > length(col)) stop('Not enough colors.')
    
    return(col[seq_len(n)])
}

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
        umap <- DimPlot(seu, reduction = 'umap', cols = kellyPal(), pt.size = umap.pt) +
            ggtitle(clustNames[i])
        tsne <- DimPlot(seu, reduction = 'tsne', cols = kellyPal(), pt.size = tsne.pt) +
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

#' Map ensembl IDs to gene names
#' @param gcm gene-cell matrix with rownames as genes
#' @param map map with ensembl gene id, transcript id and gene name as columns
#' @return gcm with ensemblID-gene as rownames
#' @export
ens2gene <- function(gcm, map) {
    ens <- gsub('\\..*', '', rownames(gcm))
    key <- gsub('\\..*', '', map$V1)
    id <- match(ens, key)
    val <- map$V3[id]
    rownames(gcm) <- paste0(rownames(gcm), '-', val)
    return(gcm)
}

#' Map FindAllMarkers output to a reference
#' @param features FindAllMarkers output
#' @param ref reference panel 
#' @importFrom dplyr as_tibble filter `%>%`
#' @return list of major cell class and tabulated hits
#' @export
returnHits <- function(features, ref) {
    
    upgenes <- features %>% dplyr::filter(avg_log2FC > 0)
    if (length(upgenes) == 0) {
        stop('No upregulated genes.')
    } else {
        markers <- upgenes$sym
        refHit <- as_tibble(ref) %>% dplyr::filter(Gene %in% markers, 
                                                   Log2FoldChange > 0)
        hitTab <- table(refHit$Ident)
        hitName <- names(which.max(hitTab))
    }
    
    return(list(MajorCT = hitName, AllHits = hitTab))
}

#' Parse RSEM genes.results to a GCM
#' @param path path of directory
#' @param field field (FPKM / TPM)
#' @return gcm
#' @export
parseRSEM <- function(path, field = 'FPKM') {
    files <- list.files(path, full.names = TRUE)
    quant <- lapply(files, function(x) {
        name <- gsub('.*/|.genes.results', '', x)
        message(name)
        x <- read.delim(x)
        x <- x[, c('gene_id', field)]
        names(x) <- c('gene_id', name)
        x})
    gcm <- do.call(cbind.data.frame, quant)
    rownames(gcm) <- gcm$gene_id
    gcm <- gcm[, -grep('gene_id', names(gcm))]
    gcm <- as.matrix(gcm)
    return(gcm)
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
combPlot <- function(seu, 
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
    
    umap <- DimPlot(seu, reduction = 'umap', cols = kellyPal(), 
                    label = TRUE, repel = TRUE, pt.size = umap.pt)
    tsne <- DimPlot(seu, reduction = 'tsne', cols = kellyPal(), 
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
condByIdent <- function(seu, 
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
    
}

