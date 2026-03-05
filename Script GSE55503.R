#Dataset: GSE55503 (Breasts Cancer in siControl vs Breast Cancer siPRKCD)
#Platform: Microarray (Illumina GPL10558)
#Problem Statement: Mengidentifikasi gen yang mengalami perubahan ekspresi 
# antara breast cancer siControl dan siPRKCD (Protein kinase C delta)

#Persiapan Lingkungan (Memanggil Library/Unduh Dahulu)
#1. Instal GEOuery
if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

#2. Instal GEOquerry LIMMA
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

#3. Install paket CRAN untuk visualisasi dan manipulasi data 
# a. phetmap: heatmap ekspresi gen 
# b. ggplot2: grafik (volcano plot)
# c. dplyr: manipulasi tabel data 
install.packages(c("pheatmap", "ggplot2", "dplyr"))

#4. Instal UMAP
#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#5. Instal annotation package sesuai platform
# a. Cek Platform GSE55503
annotation(gset)
# b. Diperoleh: GPL10558 [Illumina HumanHT-12 V4.0 expression beadchip]
BiocManager::install("illuminaHumanv4.db")

#6. #Install package enrichment
BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"),
                     ask = FALSE, update = FALSE)
n

#7. Memanggil library
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(illuminaHumanv4.db)
library(AnnotationDbi)
library(umap)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)


#Pengambilan Data dari GEO
options(download.file.method.GEOquery = "libcurl")
options(RCurlOptions = list(httpversion = 1))
options(timeout = 600)
gset <- getGEO("GSE55503", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

#Metadata
head(pData(gset))

#Buat Group
group <- ifelse(grepl("siControl", gset$title, ignore.case = TRUE),
                "Control",
                "siPRKCD")

gset$group <- factor(group)

table(gset$group)
#Hasilnya Harus 6 Control VS 6 siPRKCD

#Matriks Ekspresi/Data Ekspresi
ex <- exprs(gset)

#Log2FC Transform
#Log2 digunakan untuk:
#1. Menstabilkan variasi
#2. Mendekati asumsi model linear
#3. Memudahkan interpretasi log fold change
#quantile(): menghitung nilai kuantil (persentil)
#as.numeric(): mengubah hasil quantile (yang berupa named vector)
#menjadi vektor numerik biasa agar mudah dibandingkan

qx <- as.numeric(quantile(ex, c(0,0.25,0.5,0.75,0.99,1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6]-qx[1] > 50 & qx[2] > 0)

#Jika True:
# Nilai <= 0 tidak boleh di-log, maka diubah menjadi NA
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

#Pembagian Grup antara siControl VS siPRCKD
group <- ifelse(grepl("siControl", gset$title),
                "Control",
                "siPRKCD")
gset$group <- factor(group)

#Desain Matriks/Kerangka Statistik
gset$group <- relevel(gset$group, ref = "Control")
design <- model.matrix(~ gset$group)
colnames(design)

#Analisis LIMMA
fit <- lmFit(ex, design)
fit <- eBayes(fit)

topTableResults <- topTable(fit, coef = 2, number = Inf)
head(topTableResults)

summary(topTableResults$logFC)
summary(topTableResults$adj.P.Val)

deg <- topTableResults[
  topTableResults$adj.P.Val < 0.05 &
    abs(topTableResults$logFC) > 1,
]

nrow(deg)
#Hasilnya 52 yang artinya jumlah DEG teratasnya adalah 52#

#Cek Arah Perubahan
table(deg$logFC > 0)
  #TRUE → upregulated di siPRKCD
  #FALSE → downregulated di siPRKCD

#Jika TRUE/UPREGULATED:
#Setelah PRKCD di-knockdown, ekspresi gen meningkat.
#Hal ini mengindikasikan bahwa PRKCD kemungkinan berperan sebagai 
# repressor terhadap gen tersebut.
#Jika FALSE/DOWNREGULATED:
#Setelah PRKCD di-knockdown, ekspresi gen menurun.
#Hal ini mengindikasikan bahwa PRKCD kemungkinan berperan sebagai 
# aktivator terhadap gen tersebut.


#FILTER DEG
deg_filtered <- topTableResults %>%
  filter(adj.P.Val < 0.05 & abs(logFC) > 1)

#VOLCANO PLOT
volcano_data <- topTableResults
volcano_data$status <- "NO"

volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"

ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN"="blue","NO"="grey","UP"="red")) +
  theme_minimal() +
  ggtitle("Volcano Plot GSE55503")

#Karena lebih banyak gen yang mengalami penurunan ekspresi setelah 
# knockdown PRKCD, hal ini mengindikasikan bahwa PRKCD kemungkinan 
# berperan sebagai regulator positif (aktivator) terhadap banyak gen 
# dalam sistem ini.

#HEATMAP TOP 50
top50 <- head(topTableResults[order(topTableResults$adj.P.Val), ], 50)

mat_heatmap <- ex[rownames(top50), ]

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(mat_heatmap,
         scale = "row",
         annotation_col = annotation_col,
         show_colnames = FALSE,
         main = "Top 50 DEG")

#Di bagian atas heatmap terdapat Control biru dan siPRKCD merah yang 
# artinya gen tersebut naik di siPRKCD (upregulated).
# Di bagian bawah terdapat Control (merah) dan siPRKCD (biru) yang
# artinya gen tersebut turun di siPRKCD (downregulated).
  
#Jika dikaitkan dengan hasil Volcano Plot sebelumnya, terdapat 12 gen naik
# dan 40 gen turun yang dapat terlihat bahwa mayoritas memang turun
# dan terlihat di heatmap.


#ANOTASI NAMA GEN 
#Mengambil ID probe dari hasil DEG
probe_ids <- rownames(topTableResults)

#Mapping probe -> gene symbol & gene name
gene_annotation <- AnnotationDbi::select(
  illuminaHumanv4.db,
  keys = probe_ids,
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

#Gabungkan dengan hasil limma
topTableResults$PROBEID <- rownames(topTableResults)

topTableResults <- merge(
  topTableResults,
  gene_annotation,
  by = "PROBEID",
  all.x = TRUE
)

#Cek hasil
head(topTableResults[, c("PROBEID", "SYMBOL", "GENENAME")])
#Top Table Results
#Merupakan daftar gen/probe yang paling signifikan secara statistik 
# setelah analisis differential expression.
# Beberapa top table results yang diperoleh dapat dilihat pada tabel "topTableResults"


#BOXPLOT DISTRIBUSI NILAI EKSPRESI 
#Boxplot digunakan untuk:
#- Mengecek distribusi nilai ekspresi antar sampel
#- Melihat apakah ada batch effect
#- Mengevaluasi apakah normalisasi/log-transform sudah wajar

#Set warna berdasarkan grup
group_colors <- as.numeric(gset$group)

boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)

legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

#Penjelasan Singkat dari Hasil Box Plot
# Distribusi terlihat merata antara kotak hitam dan pink yang dapat terlihat 
# semua kotaknya berada pada ketinggian yang sama/sejajar


#DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT) 
#Density plot menunjukkan sebaran global nilai ekspresi gen
#Digunakan untuk:
#- Mengecek efek log-transform
#- Membandingkan distribusi antar grup

#Gabungkan ekspresi & grup ke data frame
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)

ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen",
    x = "Expression Value (log2)",
    y = "Density"
  )
#Penjelasan SIngkat dari Hasil
#1. Distribusi ekspresi gen pada Control dan siPRKCD relatif mirip.
#2. Bentuk kurvanya ke arah kanan dan ada yang meninggi di satu titik
#   yang dapat diartikan bahwa banyak gen memiliki ekspresi rendah–sedang
#   dengan sedikit gen yang sangat tinggi ekspresinya.


#UMAP (VISUALISASI DIMENSI RENDAH)
#UMAP digunakan untuk:
#- Mereduksi ribuan gen menjadi 2 dimensi
#- Melihat pemisahan sampel secara global
#- Alternatif PCA (lebih sensitif ke struktur lokal)

umap_input <- t(ex)

#Jalankan UMAP
umap_result <- umap(umap_input, n_neighbors = 5)
#Note: n_neighbour dibuat 5 (atau rentang 5 ≤ x < 12) karena data UMAP 
#      hanya 12, 'umap_input' digunakan untuk data yang >15

#Simpan hasil ke data frame
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)

#Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot Sampel Berdasarkan Ekspresi Gen",
    x = "UMAP 1",
    y = "UMAP 2"
  )

#Penjelasan Singkat:
#1. UMAP menunjukkan adanya pemisahan sampel yang kemungkinan besar 
#   dipengaruhi oleh perbedaan cell line. Namun demikian, analisis 
#   diferensial tetap mengidentifikasi gen-gen yang berubah signifikan antara kondisi siControl dan siPRKCD.
#2. Perbedaan jenis cell line berdasarkan database
#   berupa MDA-MB-468 dan BT-549 yang masing-masing terdapat 3 siControl
#   dan 3 siPRKCD, sehingga total ada 12 data gen.


#ENRICHMENT (GO & KEGG Pathway)
# Menggunakan P.Value dan logFC > 0.5 
# untuk menangkap lebih banyak gen karena sinyal siPRKCD lemah.
deg <- topTableResults[topTableResults$P.Value < 0.05 & abs(topTableResults$logFC) > 0.5, ]

# Cek jumlah gen yang lolos (> 10-20 gen)
message(paste("Jumlah gen untuk enrichment:", nrow(deg)))
# Hasilnya ada 1001

# Ambil gene symbol dan bersihkan
gene_symbols <- unique(deg$SYMBOL)
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

# Convert SYMBOL -> ENTREZID
gene_entrez <- bitr(gene_symbols,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

entrez_list <- gene_entrez$ENTREZID

# GO Enrichment
# GO (Gene Ontology) adalah sistem klasifikasi yang digunakan untuk 
# mengelompokkan gen berdasarkan fungsi biologisnya.

ego <- enrichGO(gene          = entrez_list,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP", #Biological Process: proses biologis yang dilakukan gen
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.1,
                qvalueCutoff  = 0.2)


# KEGG Enrichment
# KEGG (Kyoto Encyclopedia of Genes and Genomes) adalah database yang 
# memetakan gen ke dalam jalur biologis (pathway).Hubungan antar gen 
# dalam satu jalur reaksi atau sinyal biologis yang lengkap.

ekegg <- enrichKEGG(gene         = entrez_list,
                    organism     = "hsa",
                    pvalueCutoff = 0.1)

#VISUALISASI
# Fungsi dotplot akan error jika objek kosong, jadi diberi kondisi 'IF'
if(!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
  dotplot(ego, showCategory = 15) + ggtitle("GO Enrichment (Biological Process)")
} else {
  message("GO Enrichment tetap kosong: Gen signifikan terlalu sedikit untuk membentuk jalur.")
}

if(!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  dotplot(ekegg, showCategory = 15) + ggtitle("KEGG Pathway Enrichment")
} else {
  message("KEGG Enrichment tetap kosong.")
}