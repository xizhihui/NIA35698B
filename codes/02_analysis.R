library(DESeq2)
library(ggplot2)
library(dplyr)

geneinfo <- read.delim("GRCm38_STAR_2.7.10b/geneInfo.tab", header = F, skip = 1)
colnames(geneinfo) <- c("geneid", "genename", "genetype")
samples <- list.files("02_alignment")
names(samples) <- gsub("^.*_", "", samples)
cutoff_pval <- 0.05
cutoff_fc <- 1
dir.create("03_DESeq2")

## load counts ----
count_files <- file.path("02_alignment", samples, "ReadsPerGene.out.tab")
names(count_files) <- names(samples)
counts <- lapply(count_files, function(cfile) {
    count <- read.delim(cfile, skip = 4, header = FALSE)
    colnames(count) <- c("geneid", "unstrand", "strand_forward", "strand_reverse")
    rownames(count) <- count$geneid
    count <- count[, "unstrand", drop = FALSE]
    colnames(count) <- names(count_files)[count_files == cfile]
    return(count)
}) %>% do.call(what = cbind)


## DESeq2 ----
coldata <- data.frame(
    sample = names(samples),
    group = factor(gsub("\\d$", "", names(samples)), c("A","B","C","D")),
    rownames = names(samples)
)
dds <- DESeqDataSetFromMatrix(
    as.matrix(counts[, coldata$sample]),
    coldata,
    design = ~ group
)
dds <- dds[rowSums(countds(dds)) >= 100, ]
dds <- DESeq(dds)
dds_trans <- rlog(dds)

comparelist <- combn(levels(coldata$group), m = 2, simplify = FALSE)
names(comparelist) <- sapply(comparelist, function(x) paste(rev(x), collapse = " vs "))

plotPCA(dds_trans, intgroup = c("group", "sample")) + 
    coord_cartesian() + 
    cowplot::theme_cowplot() +
    geom_label(aes(label = sample, fill = gsub(":.*$", "", group)), color = "white") +
    theme(legend.position = "none")
ggsave("03_DESeq2/pca_plot.png", width = 6, height = 6)


markers <- lapply(comparelist, function(comp) {
    df <- results(dds, contrast = c("group", comp[2], comp[1]))
    df <- as.data.frame(df)
    df$compare <- paste(rev(comp), collpase = " vs ")
    df$geneid <- rownames(df)
    
    df$updown <- "others"
    df$updown[df$log2FoldChange >= cutoff_fc & df$padj <= cutoff_pval] <- "Sig_UP"
    df$updown[df$log2FoldChange >= cutoff_fc & df$padj > cutoff_pval] <- "UP"
    df$updown[df$log2FoldChange <= -cutoff_fc & df$padj <= cutoff_pval] <- "Sig_DOWN"
    df$updown[df$log2FoldChange <= -cutoff_fc & df$padj > cutoff_pval] <- "DOWN"

    df <- merge(df, geneinfo, by = "geneid")
    df
})

## save counts and transformed counts ----
sapply(markers, function(df) {
    outfile <- sprintf("03_DESeq2/DEGs_%s.csv", gsub(" ", "_", df$compare[1])))
    write.csv(df, file = outfile)
})
write.csv(counts(dds, normalized = TRUE), file = "03_DESeq2/DESeq2_normalized_counts.csv")
counts_tpm <- calcTPM(counts, "Mus_musculus.GRCm38.93.gtf")
write.csv(counts_tpm, file = "03_DESeq2/counts_TPM.csv")


## plot volcanos ----
volcanoplot <- function (
    markers, fc.col = "log2FoldChange", p.col = "padj",
    p.maxlimit = NULL, p.inflimit = FALSE, 
    gene.col = "genename", label.gene = NULL, label.n = 20,
    label.by = c("p", "fc"), sigcolors = NULL)
{
    sigcolors <- c(
        others = "gray",
        DOWN = "steelblue", Sig_DOWN = "blue", 
        UP = "#C44221", Sig_UP = "red"
    )
    
    label.by <- match.arg(label.by)
    markers[, c("FC", "PADJ", "GENE")] <- markers[, c(fc.col, p.col, gene.col)]
    markers$logpadj <- -log10(markers[, p.col])
    if (p.inflimit) {
        markers$logpadj <- -log10(markers[, p.col])
        logpadj_inf <- which(markers$logpadj == Inf)
        logpadj_max <- max(markers$logpadj[-logpadj_inf])
        markers$logpadj[logpadj_inf] <- logpadj_max + 1:length(logpadj_max)
    }

    if (!is.null(p.maxlimit)) {
        markers$logpadj <- ifelse(
            markers$logpadj > p.maxlimit, 
            (markers$logpadj - p.maxlimit) * 0.1 + p.maxlimit,
            markers$logpadj
        )
    }
    if (!is.null(label.gene)) {
        tolabel <- subset(markers, subset = GENE %in% label.gene)
    } else if (label.by == "p") {
        tolabel <- top_n(markers, n = label.n, wt = logpadj)
    } else {
        tolabel <- top_n(markers, n = label.n, wt = abs(FC))
    }
    breaks <- seq(0, max(markers$FC), by = 2)
    breaks <- unique(c(-breaks, breaks, -1, 1, 0.5, -0.5))
    ggplot(markers, aes(x = FC, y = logpadj, color = updown)) + 
        geom_point() +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
        geom_vline(xintercept = c(-1, -0.5, 0.5, 1), linetype = "dashed", color = "gray") + 
        theme_bw() + scale_color_manual(values = sigcolors) + 
        scale_x_continuous(breaks = breaks) +
        ggrepel::geom_text_repel(aes(label = GENE), data = tolabel, show.legend = F) + 
        labs(y = "-log10(p.adjust)", x = fc.col, color = "Significance")
}
add_arrow <- function(plot, x, xend, y, label, both = FALSE, hjust = 0.5, vjust = 0.5) {
    if (both) {
        x <- c(x, -x)
        xend <- c(xend, -xend)
        stopifnot(length(label) == 2)
    }
    plot + annotate(
        "segment", x = x, xend = xend, y = y, yend = y,
        color = "gray50", arrow = arrow(type = "open")
    ) + annotate(
        "text", x = (x + xend) / 2, y = y, 
        label = label, hjust = hjust, vjust = vjust, size = 12, family = "Arial"
    )
}
# just plot DEGs for B vs A here
tolabel <- subset(markers$B_vs_A, subset = updown %in% c("Sig_UP", "Sig_DOWN"))$genename
volcanoplot(markers$B_vs_A, label.gene = tolabel) %>%
    add_arrow(x = -2, xend = -6, y = 0, label = c("A", "B"), both = T, vjust = 0)
ggsave("03_DESeq2/volcano_B_vs_A.png", width = 8, height = 8)

## plot heatmaps ----
library(ComplexHeatmap)
normcount <- "03_DESeq2/DESeq2_normalized_counts.csv" %>%
    read.csv(row.names = 1) %>% 
    as.matrix()
aimgenes <- markers$B_vs_A %>%
    dplyr::filter(updown %in% c("Sig_UP", "Sig_DOWN") | genename == "Siglece") %>%
    dplyr::select(geneid, genename)
mat <- pheatmap:::scale_mat(normcount[aimgenes$geneid, ], "row") %>%
    magrittr::set_rownames(make.unique(aimgenes$genename))
hm <- Heatmap(
    mat, name = " ",
    show_row_names = F,
    cluster_columns = F,
    cluster_rows = T,
    column_split = plyr::mapvalues(
        dds$group,
        from = LETTERS[1:4],
        to = c(
          "WT BMDM\nXBB.1", "WT BMDM\nXBB.1(S375)",
          "Siglece-/- BMDM\nXBB.1", "Siglece-/- BMDM\nXBB.1(S375)"
        )
    )
)
pdf(file = "03_DESeq2/heatmap.pdf", width = 2.5, height = 2.5, family = "Arial")
    pushViewport(viewport(gp = gpar(fontsize = 12, fontfamily = "Arial")))
    draw(hm, newpage = FALSE)
    popViewport()
dev.off()

## enrichments ----
siggene <- lapply(markers, function(marker) {
    split(make.unique(marker$genename), marker$updown)[c("Sig_UP", "Sig_DOWN")]
})
en <- recur_enrich(siggene, species = "mouse", use_internal = T)
endf <- unlist_enrichlist(en)
str_wrap <- function(x) stringr::str_wrap(x)

comp_current <- "C_vs_A"
tempdf <- endf %>%
    dplyr::mutate(newgroup = group2) %>%
    tidyr::separate(
      col = "newgroup",
      into = c("compare_g1", "compare_g2", "updown"),
      sep = "\\.|_vs_"
    ) %>% dplyr::mutate(groupname = ifelse(updown == "Sig_UP", compare_g1, compare_g2)) %>%
    dplyr::filter(grepl(comp_current, group))
tempdf %>%
    dplyr::select(-group,-group2,-compare_g1, -compare_g2) %>%
    write.csv(file = sprintf("03_DESeq2/%s/enrichment_hyper.csv", comp_current))
tempdf %>%
    dplyr::filter(Count >= 3) %>%
    dotplot_enrich() +
    scale_y_discrete(label = str_wrap) +
    facet_grid(category ~ groupname, space = "free", scales = "free")
ggsave(sprintf("03_DESeq2/%s/enrichment_hyper.png", comp_current), width = 9, height = 9)
