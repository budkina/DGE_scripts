library(optparse)
library(edgeR)
library(DESeq2)

edger <- function(dge_list_filtered, group)
{
    dge_list_norm <- calcNormFactors(dge_list_filtered)
    design_matrix <- model.matrix(~group)
    dge_list_disp <- estimateDisp(dge_list_norm,design_matrix)
    fit <- glmQLFit(dge_list_disp,design_matrix)
    lrt <- glmQLFTest(fit,coef=2)
    top_tags <- topTags(lrt, n = nrow( lrt$table ), sort.by = "PValue", p.value = 0.05)$table
    result <- top_tags[which(top_tags$FDR < 0.05),]
    write.csv(result, "edgeR_result.csv")
}

voom_limma <- function(dge_list_filtered, group)
{
    dge_list_norm <- calcNormFactors(dge_list_filtered)
    design_matrix <- model.matrix(~group)
    y <- voom(dge_list_norm, design_matrix)
    fit <- lmFit(y, design_matrix)
    fit <- eBayes(fit)
    top_table <- topTable(fit, sort.by = "P", n = Inf)
    result <-top_table[which(top_table$adj.P.Val < 0.05),]
    write.csv(result, "voom_result.csv")

}

deseq2 <- function(dge_list_filtered, group)
{
    condition <- factor(group)
    colData <- data.frame(row.names=colnames(dge_list_filtered), condition)
    deseqdata <- DESeqDataSetFromMatrix(countData=dge_list_filtered$counts, colData=colData, design=~condition)
    dds <-DESeq(deseqdata)
    top_table <- lfcShrink(dds, coef=paste0("condition_",levels(condition)[2],"_vs_",levels(condition)[1]), type="apeglm")
    top_table <- top_table[order(top_table$padj),]
    result <- top_table[which(top_table$padj < 0.05),]
    write.csv(result, "deseq2_result.csv")
}

option_list = list(
    make_option(c("-w", "--working_dir"), type="character", default=NULL,
        help="Working directory", metavar="character"),
    make_option(c("-d", "--design"), type="character", default=NULL,
        help="Design dataframe: bam_filename, group_id columns", metavar="character")
    )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$design)){
    print_help(opt_parser)
    stop("Design matrix should be supplied", call.=FALSE)
}

if (is.null(opt$working_dir)){
    print_help(opt_parser)
    stop("Working directory should be supplied", call.=FALSE)
}

# Set working directory
setwd(opt$working_dir)

counts <- read.table("counts_table.csv", header=TRUE, sep=' ')
design <- read.table(opt$design, header=TRUE, sep=' ')
group <- design$group_id

# filter
dge_list <- DGEList(counts=counts,group=group)
keep <- filterByExpr(dge_list)
dge_list_filtered <- dge_list[keep,,keep.lib.sizes=FALSE]

# DGE
edger(dge_list_filtered, group)
voom_limma(dge_list_filtered, group)
deseq2(dge_list_filtered, group)
