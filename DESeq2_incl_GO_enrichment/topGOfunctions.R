getTopGO_exp_matching <-
    function(DEgeneSet,
             dds, nR = 10) {
        
        
        df <- data.frame(sign = as.numeric(rownames(dds) %in% as.character(DEgeneSet)),
                         geneExp = rowMeans(DESeq2::counts(dds, T)))
        
        match_res <-
            MatchIt::matchit(
                sign ~ geneExp,
                df,
                method = "nearest",
                distance = "mahalanobis",
                ratio = nR
            )
        
        background <- as.vector(match_res$match.matrix)
        background
    }

plotMatchingRes <- function(DEgeneSet, matchRes, dds) {
    meanExp <- rowMeans(DESeq2::counts(dds, T))
    dplyr::bind_rows(
        "all" = data.frame(exp = meanExp[rownames(dds)]),
        "DE" = data.frame(exp = meanExp[DEgeneSet]),
        "Background" = data.frame(exp = meanExp[matchRes]),
        .id = "group"
    ) %>%
        ggplot2::ggplot() +
        ggplot2::geom_density(ggplot2::aes(x = exp, color = group)) +
        ggplot2::scale_x_log10(name = "log10(Mean expression)")
}

makeTopGOobject <- function(DEgeneSet, matchRes, dds, nodeSize = 5){
    # df <- df[which(!rownames(df) %in% background), ]
    # background <- unique(na.exclude(background))
    # ## total number of matched de_genes
    # length(background)
    
    geneIDs <- rownames(dds)
    inUniverse <- geneIDs %in% c(DEgeneSet,  matchRes)
    inSelection <-  geneIDs %in% DEgeneSet
    
    alg <- factor(as.integer(inSelection[inUniverse]))
    names(alg) <- geneIDs[inUniverse]
    
    tgd <-
        new(
            "topGOdata",
            ontology = "BP",
            allGenes = alg,
            nodeSize = nodeSize,
            annot = topGO::annFUN.gene2GO,
            gene2GO = GOList
        )
    tgd
}

exportGOtable <-
    function(tgd_object,
             orderTest = "Fisher.weight01",
             n = 25,
             joinfun = "intersect") {
        resultTopGO.weight01 <-
            topGO::runTest(tgd_object, algorithm = "weight01", statistic = "Fisher")
        #resultTopGO.elim <-
        #runTest(tgd_object, algorithm = "elim", statistic = "Fisher")
        #resultTopGO.classic <-
        #runTest(tgd_object, algorithm = "classic", statistic = "Fisher")
        #resultTopGO.parentchild <-
        #runTest(
        #tgd_object,
        #algorithm = "parentchild",
        #statistic = "Fisher",
        #joinFun = joinfun
        #)
        topGO::GenTable(
            tgd_object,
            #Fisher.elim = resultTopGO.elim,
            #Fisher.classic = resultTopGO.classic,
            Fisher.weight01 = resultTopGO.weight01,
            #Fisher.pc = resultTopGO.parentchild,
            orderBy = orderTest,
            topNodes = min(resultTopGO.weight01@geneData[4], n)
        )
    }

extractSigGenes<-function(Term, tgd, DEset){
    allgenes<-genesInTerm(tgd)[[Term]]
    DEset$geneid[DEset$geneid %in% allgenes]
}


runMatchedTopGoAnalysis <- function(DEgeneSet,
                                    dds,
                                    GOdb,
                                    onts = c("MF", "BP", "CC"),
                                    nR = 10,
                                    nodeSize = 5,
                                    plot = FALSE) {
    match <- getTopGO_exp_matching(DEgeneSet, dds, nR = nR)
    if(plot == TRUE) {print(plotMatchingRes(DEgeneSet, matchRes = match, dds = dds))}
    tgd <- makeTopGOobject(DEgeneSet = DEgeneSet, matchRes = match, dds = dds, nodeSize = nodeSize) 
    tgd
}

runTopGoAnalysis <- function(DEgeneSet,
                             dds,
                             GOdb,
                             onts = c("MF", "BP", "CC"),
                             nR = 10,
                             nodeSize = 5)
{
    match <- rownames(dds)
    tgd <- makeTopGOobject(DEgeneSet = DEgeneSet, matchRes = match, dds = dds, nodeSize = nodeSize) 
    tgd
}
