
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

#' Kelly and Alphabet palette
#' @param n number of colors
#' @return colors
#' @export
scpalette <- function(n=45, alpha = 200) {
    
    maxCV <- 255
    col <- c(
    
        ## Kelly palette
        rgb(0, 0, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 179, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(128, 62, 117, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 104, 0, maxColorValue = maxCV, alpha = alpha),
        
        rgb(166, 189, 215, maxColorValue = maxCV, alpha = alpha),
        rgb(193, 0, 32, maxColorValue = maxCV, alpha = alpha),
        rgb(206, 162, 98, maxColorValue = maxCV, alpha = alpha),
        rgb(0, 125, 52, maxColorValue = maxCV, alpha = alpha),
        
        rgb(246, 118, 142, maxColorValue = maxCV, alpha = alpha),
        rgb(0, 83, 138, maxColorValue = maxCV, alpha = alpha),
        rgb(129, 112, 102, maxColorValue = maxCV, alpha = alpha),
        rgb(179, 40, 81, maxColorValue = maxCV, alpha = alpha),
        
        rgb(255, 122, 92, maxColorValue = maxCV, alpha = alpha),
        rgb(83, 55, 122, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 142, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(244, 200, 0, maxColorValue = maxCV, alpha = alpha),
        
        rgb(127, 24, 13, maxColorValue = maxCV, alpha = alpha),
        rgb(147, 170, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(89, 51, 21, maxColorValue = maxCV, alpha = alpha),
        rgb(241, 58, 19, maxColorValue = maxCV, alpha = alpha),
        rgb(35, 44, 22, maxColorValue = maxCV, alpha = alpha),
        
        ## Alphabet palette
        rgb(240, 163, 255, maxColorValue = maxCV, alpha = alpha),
        rgb(0, 117, 220, maxColorValue = maxCV, alpha = alpha),
        rgb(153, 63, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(76, 0, 92, maxColorValue = maxCV, alpha = alpha),

        rgb(0, 92, 49, maxColorValue = maxCV, alpha = alpha),
        rgb(43, 206, 72, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 204, 153, maxColorValue = maxCV, alpha = alpha),
        rgb(148, 255, 181, maxColorValue = maxCV, alpha = alpha),
        
        rgb(143, 124, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(157, 204, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(194, 0, 136, maxColorValue = maxCV, alpha = alpha),
        rgb(0, 51, 128, maxColorValue = maxCV, alpha = alpha),
        
        rgb(255, 164, 5, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 168, 187, maxColorValue = maxCV, alpha = alpha),
        rgb(66, 102, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 0, 0, maxColorValue = maxCV, alpha = alpha),
        
        rgb(94, 241, 242, maxColorValue = maxCV, alpha = alpha),
        rgb(0, 153, 143, maxColorValue = maxCV, alpha = alpha),
        rgb(224, 255, 102, maxColorValue = maxCV, alpha = alpha),
        rgb(116, 10, 255, maxColorValue = maxCV, alpha = alpha),
        
        rgb(153, 0, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 255, 128, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 225, 0, maxColorValue = maxCV, alpha = alpha),
        rgb(255, 80, 5, maxColorValue = maxCV, alpha = alpha)
        
    )
    
    if (n > length(col)) stop('Not enough colors.')
    
    return(col[seq_len(n)])
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
#' @importFrom data.table fread
#' @return gcm
#' @export
parseRSEM <- function(path, field = "FPKM") {
    files <- list.files(path, full.names = TRUE)
    quant <- lapply(files, function(x) {
        name <- gsub(".*/|.genes.results", "", x)
        message(name)
        x <- fread(x, select = c('gene_id', field))
        names(x) <- c("gene_id", name)
        x
    })
    message('Merging data')
    gcm <- do.call(cbind.data.frame, quant)
    rownames(gcm) <- gcm$gene_id
    gcm <- gcm[, -grep("gene_id", names(gcm))]
    gcm <- as.matrix(gcm)
    return(gcm)
}

#' Parse STAR <Sample>.Log.final.out files prefixed by sample name
#' @param path path of directory
#' @importFrom data.table fread
#' @return mapping statistics
#' @export
parseSTARlog <- function(path) {
    files <- list.files(path, full.names = TRUE)
    message('Parsing files')
    quant <- lapply(files, function(x) {
        name <- gsub(".*/|.Log.final.out", "", x)
        x <- fread(x, sep = '\t', fill = TRUE, 
                   header = FALSE, blank.lines.skip = TRUE)
        x$V1 <- sub(' \\|', '', x$V1)
        x$V2 <- sub('%', '', x$V2)
        x <- x[V2 != '']
        names(x) <- c('Stat', name)
        x
    })
    message('Merging ')
    stat <- do.call(cbind.data.frame, stat)
    fields <- stat[,1]
    out <- data.frame(out[,-grep('Stat', colnames(out))])
    rownames(out) <- fields
    return(out)
}
