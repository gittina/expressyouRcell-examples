setwd("hereyouhavetowriteyourfolderpath")

library(data.table)
library(ggplot2)
library(expressyouRcell)

load("ensg_dt.RData")

# create association genes and localizations (just once) ----

gene_loc_table <- map_gene_localization(gene_set = "gencode.vM22.primary_assembly.annotation.gtf")
#save(gene_loc_table, file = "gene_loc_table.RData")

# figure 1) trap seq data from neurons ----

trap_data <- as.data.table(readxl::read_excel("data/dataset_trap_neurons_figure2.xlsx", skip = 3))
colnames(trap_data) <- c("gene_symbol","description","refSeq","entrez","log2FC","FDR",
                         "basal average normalized counts","LTP average normalized counts","Time point")


time_points <- unique(trap_data$`Time point`)

trap_l <- list()
for (tp in time_points){
  trap_l[[tp]] <- trap_data[`Time point` == tp]
}

trap_fc <- color_cell(timepoint_list = trap_l,
                      pictogram = "neuron",
                      gene_loc_table = gene_loc_table,
                      group_by="class",
                      thr = 0.4,
                      pval_col="FDR",
                      pval_thr=0.1,
                      coloring_mode="mean",
                      col_name="log2FC",
                      colors = list("+" = c("#eaf3ea", "#307e2d"),
                                    "-" = c("#f3eaea", "#7e302d")))

# generate animations for up and down genes with fold change option
animate(data = trap_fc,
        timepoints=names(trap_fc$final_dt)[names(trap_fc$final_dt) %like% '\\-'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 6, width = 12,
        filename = "movie_S1_trap_fc_down",
        names = c("30'","60'","120'"),
        format = "video")

animate(data = trap_fc,
        timepoints=names(trap_fc$final_dt)[names(trap_fc$final_dt) %like% '\\+'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 6, width = 12,
        filename = "movie_S2_trap_fc_up",
        names = c("30'","60'","120'"),
        format = "video")

trap_enr <- color_cell(timepoint_list = trap_l,
                       pictogram = "neuron",
                       gene_loc_table = gene_loc_table,
                       group_by="class",
                       thr = 0.4,
                       pval_col="FDR",
                       pval_thr=0.1,
                       coloring_mode="enrichment",
                       col_name="log2FC")

# generate animations for up and down genes with enrichment option

animate(data = trap_enr,
        timepoints=names(trap_enr$final_dt)[names(trap_enr$final_dt) %like% '\\-'],
        seconds = 3, fps = 15,
        input_dir = "movies", height = 6, width = 12,
        filename = "movie_S3_trap_enr_down",
        names = c("30'","60'","120'"),
        format = "video")

animate(data = trap_enr,
        timepoints=names(trap_enr$final_dt)[names(trap_enr$final_dt) %like% '\\+'],
        seconds = 3, fps = 15,
        input_dir = "movies", height = 6, width = 12,
        filename = "movie_S4_trap_enr_up",
        names = c("30'","60'","120'"),
        format = "video")



# figure s1) mus musculus nine timepoints ----

# from https://www.nature.com/articles/s41467-018-04559-0
# TPM measured by salmon tool
# mus musculus ensembl annotatation version 97

# all counts table
mm_9timepoints <- as.data.table(readxl::read_excel("data/dataset1_musmusculus_corticalneurons_9development_figureS1.xlsx"))

# filtered by the authors and considered for downstream analysis
subset <- as.data.table(readxl::read_excel("data/dataset1_musmusculus_corticalneurons_9development_supplementary_figureS1.xlsx"))
subset <- subset[!is.na(subset$`Gene Symbol`)]

cols <- split(colnames(mm_9timepoints[, -"Gene"]), tstrsplit(colnames(mm_9timepoints[, -"Gene"]), "a|b"))

# compute means of replicas
for (col in names(cols)){
  mm_9timepoints[, ((col)) := rowMeans(.SD), .SDcols = unlist(cols[col])]
}

mean9 <- mm_9timepoints[, c("Gene", names(cols)), with=FALSE]
setcolorder(mean9, neworder = c("Gene", "E14.5", "E16.5", "P0", "P4", "P7", "P15", "P30", "P110", "21M"))

# assign gene symbol
mean9 <- merge.data.table(mean9, ensg_dt[, c("external_gene_name", "ensembl_gene_id")],
                          by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)

mean9 <- mean9[external_gene_name %in% subset$`Gene Symbol`]

mean9_list <- list()
for (tp in colnames(mean9[, -c("Gene", "external_gene_name")])){
  mean9_list[[tp]] <- mean9[, c("external_gene_name", tp), with=FALSE]
  mean9_list[[tp]] <- setnames(mean9_list[[tp]], c("gene_symbol", "tpm"))
}

# create pictograms 
mm_9_static <- color_cell(timepoint_list = mean9_list,
                          pictogram = "neuron",
                          gene_loc_table = gene_loc_table,
                          coloring_mode = "mean",
                          col_name = "tpm",
                          legend=TRUE)

animate(data = mm_9_static,
        timepoints=names(mm_9_static$final_dt),
        seconds = 2, fps = 15,
        input_dir = "movies", height = 5, width = 10,
        filename = "movie_S5_neuron_development",
        names = c("E14.5","E16.5","P0","P4","P7", "P15","P30","P110","21M"),
        format = "video")

# figure s2) rna seq microglia ----
m2 <- readxl::read_excel("data/dataset3_microglia_figureS2.xlsx ",
                         sheet = "Month 2")
m4 <- readxl::read_excel("data/dataset3_microglia_figureS2.xlsx",
                         sheet = "Month 4")
m6 <- readxl::read_excel("data/dataset3_microglia_figureS2.xlsx",
                         sheet = "Month 6")
m8 <- readxl::read_excel("data/dataset3_microglia_figureS2.xlsx",
                         sheet = "Month 8")

microglia_data <- list("m2"=m2,
                       "m4"=m4,
                       "m6"=m6,
                       "m8"=m8)

microglia_data <- lapply(microglia_data, function(x) setnames(as.data.table(x), old = "GeneSymbol", new = "gene_symbol"))

microglia_data_meanfc <- color_cell(timepoint_list = microglia_data,
                                    pictogram = "microglia",
                                    gene_loc_table = gene_loc_table,
                                    col_name="foldchange",
                                    group_by="class",
                                    thr = 1.5,
                                    pval_col="pvalue",
                                    pval_thr=0.05,
                                    coloring_mode="mean",
                                    colors = list("+" = c("#eaf3ea", "#307e2d"),
                                                  "-" = c("#f3eaea", "#7e302d")))

animate(data = microglia_data_meanfc,
        timepoints=names(microglia_data_meanfc$final_dt)[names(microglia_data_meanfc$final_dt) %like% '\\-'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 7, width = 10,
        filename = "movie_S6_microglia_fc_down",
        names = c("M2","M4","M6","M8"),
        format = "video")

animate(data = microglia_data_meanfc,
        timepoints=names(microglia_data_meanfc$final_dt)[names(microglia_data_meanfc$final_dt) %like% '\\+'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 7, width = 10,
        filename = "movie_S7_microglia_fc_up",
        names = c("M2","M4","M6","M8"),
        format = "video")
 

microglia_data_enr <- color_cell(timepoint_list = microglia_data,
                                 pictogram = "microglia",
                                 gene_loc_table = gene_loc_table,
                                 group_by="class",
                                 thr = 1.5,
                                 pval_col="pvalue",
                                 pval_thr=0.05,
                                 coloring_mode="enrichment",
                                 col_name="foldchange")

animate(data = microglia_data_enr,
        timepoints=names(microglia_data_enr$final_dt)[names(microglia_data_enr$final_dt) %like% '\\-'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 7, width = 10,
        filename = "movie_S8_microglia_enr_down",
        names = c("M2","M4","M6","M8"),
        format = "video")

animate(data = microglia_data_enr,
        timepoints=names(microglia_data_enr$final_dt)[names(microglia_data_enr$final_dt) %like% '\\+'],
        seconds = 3, fps = 10,
        input_dir = "movies", height = 7, width = 10,
        filename = "movie_S9_microglia_enr_up",
        names = c("M2","M4","M6","M8"),
        format = "video")

# figure s3) clements data  ----
wound_time_df <- readxl::read_excel("data/dataset4_SC_wound_figureS3.xlsx", skip = 7)

wound_time_df <- setDT(wound_time_df)[,c("Gene ID","ratio_D1-I","ratio_D2-I",	"ratio_D4-I","ratio_D6-I"	,"ratio_D8-I","ratio_D10-I",
                                         "pval_D1-I","pval_D2-I","pval_D4-I","pval_D6-I","pval_D8-I","pval_D10-I")]

setnames(wound_time_df, new = c("gene_symbol", 
                                paste0(c("D1", "D2", "D4", "D6", "D8", "D10"), "_log2FC"),
                                paste0(c("D1", "D2", "D4", "D6", "D8", "D10"), "_pval")))

d_times <- c("D1", "D2", "D4", "D6", "D8", "D10")
wound_list <- list()

for (t in d_times){
  wound_list[[t]] <- wound_time_df[, c("gene_symbol", colnames(wound_time_df)[colnames(wound_time_df) %like% t]), 
                                   with=FALSE]
  wound_list[[t]][,  paste0(t, "_log2FC") := as.numeric(get(paste0(t, "_log2FC")))]
  wound_list[[t]][,  paste0(t, "_pval") := as.numeric(get(paste0(t, "_pval")))]
}

wound_list[["D1"]]$D10_log2FC <- NULL
wound_list[["D1"]]$D10_pval <- NULL

lapply(wound_list,
       function(x) setnames(x, old = colnames(x)[colnames(x) %like% "log2FC"], new = "log2FC"))
lapply(wound_list,
       function(x) setnames(x, old = colnames(x)[colnames(x) %like% "pval"], new = "pval"))

wound_enr <- color_cell(timepoint_list = wound_list,
                        pictogram = "cell",
                        gene_loc_table = gene_loc_table,
                        group_by="class",
                        coloring_mode="enrichment",
                        col_name = "log2FC",
                        thr = 2.5,
                        pval_col="pval",
                        pval_thr=0.05,
                        colors = list("+" = c("#eaf3ea", "#557C55"),
                                      "-" = c("#f3eaea", "#7e302d")),
                        grouping_vars = list("class"=c("+","-")))


animate(data = wound_enr,
        timepoints=names(wound_enr$final_dt)[names(wound_enr$final_dt) %like% '\\-'],
        seconds = 4, fps = 10,
        input_dir = "movies", height = 5, width = 8,
        filename = "movie_S10_wound_down_enr",
        names = d_times,
        format = "video")

animate(data = wound_enr,
        timepoints=names(wound_enr$final_dt)[names(wound_enr$final_dt) %like% '\\+'],
        seconds = 4, fps = 10,
        input_dir = "movies", height = 5, width = 8,
        filename = "movie_S11_wound_up_enr",
        names = d_times,
        format = "video")



# figure s4) doktor samples ----
doktor_scp1 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Spinal_cord_PND1_gene_expressio"))
doktor_scp5 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Spinal_cord_PND5_gene_expressio"))

doktor_bp1 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Brain_PND1_gene_expression"))
doktor_bp5 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Brain_PND5_gene_expression"))

doktor_lp1 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Liver_PND1_gene_expression"))
doktor_lp5 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                                sheet = "Liver_PND5_gene_expression"))

doktor_mp1 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                               sheet = "Muscle_PND1_gene_expression"))
doktor_mp5 <- as.data.table(readxl::read_excel("data/dataset5_SMA_rnaseq_figureS4.xlsx",
                                               sheet = "Muscle_PND5_gene_expression"))


doktor_neuro <- list("scp1"=doktor_scp1, "scp5"=doktor_scp5,
                     "bp1"=doktor_bp1, "bp5"=doktor_bp5)

doktor_peripheric <- list("lp1"=doktor_lp1, "lp5"=doktor_lp5,
                          "mp1"=doktor_mp1, "mp5"=doktor_mp5)

lapply(doktor_neuro, function(x) setnames(x, old="Gene symbol", new="gene_symbol"))
lapply(doktor_peripheric, function(x) setnames(x, old="Gene symbol", new="gene_symbol"))

doktor_all_enr_n <- color_cell(timepoint_list = doktor_neuro,
                               pictogram = "neuron",
                               gene_loc_table = gene_loc_table,
                               coloring_mode = "enrichment",
                               group_by = "class",
                               col_name = "log2FoldChange",
                               thr = 0.12,
                               pval_col="padj",
                               pval_thr=0.1)


#spinal cord down
animate(data = doktor_all_enr_n,
        timepoints=c("scp1-", "scp5-"),
        seconds = 4, fps = 10,
        input_dir = "movies", height = 5, width = 10,
        filename = "movie_S12_sc_down_enr",
        names = c("PND1","PND5"),
        format = "video")

# spinal cord up
animate(data = doktor_all_enr_n,
        timepoints=c("scp1+", "scp5+"),
        seconds = 4, fps = 10,
        input_dir = "movies", height = 5, width = 10,
        filename = "movie_S13_sc_up_enr",
        names = c("PND1","PND5"),
        format = "video")



doktor_all_enr_p <- color_cell(timepoint_list = doktor_peripheric,
                               pictogram = "cell",
                               gene_loc_table = gene_loc_table,
                               coloring_mode = "enrichment",
                               group_by = "class",
                               col_name = "log2FoldChange",
                               thr = 0.12,
                               pval_col="padj",
                               pval_thr=0.1)

# muscle down
animate(data = doktor_all_enr_p,
        timepoints=c("mp1-", "mp5-"),
        seconds = 4, fps = 10,
        input_dir = "movies", height = 6, width = 10,
        filename = "movie_S14_m_down_enr",
        names = c("PND1","PND5"),
        format = "video")

# muscle up
animate(data = doktor_all_enr_p,
        timepoints=c("mp1+", "mp5+"),
        seconds = 4, fps = 10,
        input_dir = "movies", height = 6, width = 10,
        filename = "movie_S15_m_up_enr",
        names = c("PND1","PND5"),
        format = "video")







