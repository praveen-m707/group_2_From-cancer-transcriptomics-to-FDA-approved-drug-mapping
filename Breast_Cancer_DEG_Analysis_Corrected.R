# Load breast cancer expression data
expression_data <- read.csv(
  "Breast_Cancer_Expression.csv",
  check.names = FALSE
)
# Load sample information
metadata <- read.csv("Breast_Cancer_Metadata.csv")
dim(expression_data)
dim(metadata)
table(metadata$Group)
all(colnames(expression_data)[-1] == metadata$Sample_ID)
# Check the first few probe IDs
head(expression_data$ID_REF)

# Check the expression values
summary(as.numeric(unlist(expression_data[1, -1])))

# Check whether expression columns are numeric
str(expression_data[, 2:6])

# Check duplicate probe IDs
sum(duplicated(expression_data$ID_REF))
head(expression_data$ID_REF)
summary(as.numeric(unlist(expression_data[1, -1])))
summary(as.numeric(unlist(expression_data[, -1])))
# Convert expression measurements into a matrix
expr_matrix <- as.matrix(expression_data[, -1])

# Find the largest expression value
max(expr_matrix, na.rm = TRUE)

# Find where the largest value occurs
which(expr_matrix == max(expr_matrix, na.rm = TRUE),
      arr.ind = TRUE)
# Check expression ranges for each sample
summary(colMeans(expr_matrix, na.rm = TRUE))

# Check the first few rows and columns
expr_matrix[1:5, 1:5]
# Locate the highest expression value
max_position <- which(
  expr_matrix == max(expr_matrix, na.rm = TRUE),
  arr.ind = TRUE
)

# Display its row and column
max_position

# Identify the probe
expression_data$ID_REF[max_position[1, "row"]]

# Identify the sample
colnames(expr_matrix)[max_position[1, "col"]]

# Display the actual value
expr_matrix[max_position]
expression_data$ID_REF[32934]
colnames(expr_matrix)[96]
expr_matrix[32934, 96]
summary(as.numeric(expr_matrix[32934, ]))
summary(expr_matrix[, 96])
sum(expr_matrix > 100, na.rm = TRUE)
# Check whether this sample is normal or tumour
metadata[metadata$Sample_ID == "GSM564015", ]

# Check the expression values in this sample
summary(expr_matrix[, 96])

# Count how many values are above 100 in this sample
sum(expr_matrix[, 96] > 100, na.rm = TRUE)

# Check the dimensions of our expression matrix
dim(expr_matrix)

# Check the number of samples in metadata
dim(metadata)
# Find the 10 unusually high expression values
high_values <- which(expr_matrix[, 96] > 100)

# Display their probe IDs and expression values
data.frame(
  Probe_ID = expression_data$ID_REF[high_values],
  Expression = expr_matrix[high_values, 96]
)
# Identify whether sample 96 is normal or tumour
metadata[metadata$Sample_ID == "GSM564015", ]
# Check this probe across all 186 samples
probe_values <- expr_matrix[
  expression_data$ID_REF == "A_24_P937751", 
]

# Summarise its expression values
summary(probe_values)

# Find samples where this probe has expression above 100
which(probe_values > 100)

# Display the sample IDs and expression values
data.frame(
  Sample_ID = names(probe_values)[probe_values > 100],
  Expression = probe_values[probe_values > 100]
)


# Display the 10 highest values for this probe
sort(probe_values, decreasing = TRUE)[1:10]

# Calculate the median expression
median(probe_values, na.rm = TRUE)

# Count values above 10
sum(probe_values > 10, na.rm = TRUE)

# Count missing values
sum(is.na(probe_values))
# Get all expression values from sample GSM564015
sample_96 <- expr_matrix[, 96]

# Display the 20 highest expression values
sort(sample_96, decreasing = TRUE)[1:20]

# Count values above 10
sum(sample_96 > 10, na.rm = TRUE)

# Count values above 1000
sum(sample_96 > 1000, na.rm = TRUE)

# Calculate the median and mean
median(sample_96, na.rm = TRUE)
mean(sample_96, na.rm = TRUE)

# Calculate median expression for every sample
sample_medians <- apply(expr_matrix, 2, median, na.rm = TRUE)

# Calculate mean expression for every sample
sample_means <- colMeans(expr_matrix, na.rm = TRUE)

# Display results
sample_QC <- data.frame(
  Sample_ID = colnames(expr_matrix),
  Median = sample_medians,
  Mean = sample_means
)

# Sort samples by median expression
sample_QC[order(sample_QC$Median), ]

# Compare the median of all samples
summary(sample_medians)

# Check the median of GSM564015
sample_medians["GSM564015"]
boxplot(
  expr_matrix,
  outline = FALSE,
  las = 2,
  main = "Expression Distribution Across All Samples",
  xlab = "Samples",
  ylab = "Expression Values"
)
boxplot(
  sample_medians,
  main = "Distribution of Sample Median Expression",
  ylab = "Median Expression"
)
# Check the minimum and maximum values
range(expr_matrix, na.rm = TRUE)

# Check whether any values are zero or negative
sum(expr_matrix <= 0, na.rm = TRUE)

# Examine the expression value distribution
summary(as.vector(expr_matrix))
log_expr_matrix <- log2(expr_matrix)

# Summary of log2-transformed expression values
summary(as.vector(log_expr_matrix))

# Check the range after transformation
range(log_expr_matrix, na.rm = TRUE)

# Plot transformed expression distributions
boxplot(
  log_expr_matrix,
  outline = FALSE,
  main = "Log2 Expression Distribution Across Samples",
  xlab = "Samples",
  ylab = "Log2 Expression Values"
)

log_expr_matrix[32934, 96]

save(
  expr_matrix,
  log_expr_matrix,
  metadata,
  sample_QC,
  file = "Breast_Cancer_QC_Progress.RData"
)
sample_cor <- cor(
  log_expr_matrix,
  method = "pearson",
  use = "pairwise.complete.obs"
)

# Examine correlation values
summary(as.vector(sample_cor[upper.tri(sample_cor)]))

hist(
  sample_cor[upper.tri(sample_cor)],
  breaks = 30,
  main = "Sample-to-Sample Correlation",
  xlab = "Pearson Correlation"
)

# Identify sample groups
tumor_samples <- metadata$Group == "primary breast tumor"
normal_samples <- metadata$Group == "normal breast tissue"

# Correlation among tumour samples
tumor_cor <- cor(
  log_expr_matrix[, tumor_samples],
  method = "pearson"
)

# Correlation among normal samples
normal_cor <- cor(
  log_expr_matrix[, normal_samples],
  method = "pearson"
)

# Correlation between tumour and normal samples
cross_cor <- cor(
  log_expr_matrix[, tumor_samples],
  log_expr_matrix[, normal_samples],
  method = "pearson"
)

# Summarize each correlation group
summary(tumor_cor[upper.tri(tumor_cor)])

summary(normal_cor[upper.tri(normal_cor)])

summary(as.vector(cross_cor))
# Average correlation of each tumour sample
tumor_avg_cor <- rowMeans(tumor_cor, na.rm = TRUE)

# Remove self-correlation (diagonal = 1)
diag(tumor_cor) <- NA
tumor_avg_cor <- rowMeans(tumor_cor, na.rm = TRUE)

# Display 10 least-correlated tumour samples
head(sort(tumor_avg_cor), 10)
# Remove self-correlation
diag(normal_cor) <- NA

# Average correlation of each normal sample
normal_avg_cor <- rowMeans(normal_cor, na.rm = TRUE)

# Display 5 least-correlated normal samples
head(sort(normal_avg_cor), 5)
# List the 10 unusual tumour samples
low_cor_samples <- names(head(sort(tumor_avg_cor), 10))

# Compare their median and mean expression
low_cor_QC <- sample_QC[
  sample_QC$Sample_ID %in% low_cor_samples,
]

# Add their average correlation
low_cor_QC$Average_Correlation <-
  tumor_avg_cor[low_cor_QC$Sample_ID]

# Display results
low_cor_QC[
  order(low_cor_QC$Average_Correlation),
]
# Check metadata of the 10 least-correlated tumour samples
metadata[
  match(low_cor_samples, metadata$Sample_ID),
]

  heatmap(
    sample_cor,
    symm = TRUE,
    scale = "none",
    margins = c(8, 8),
    main = "Hierarchical Clustering of Breast Cancer Samples"
  )
  
  
  # Calculate variance for each probe
  probe_variance <- apply(log_expr_matrix, 1, var, na.rm = TRUE)
  
  # Select the 500 most variable probes
  top_500 <- order(probe_variance, decreasing = TRUE)[1:500]
  
  # Create a matrix containing these probes
  pca_data <- t(log_expr_matrix[top_500, ])
  
  # Perform PCA
  pca_result <- prcomp(pca_data, center = TRUE, scale. = FALSE)
  
  # Assign colours to sample groups
  group_colors <- ifelse(
    metadata$Group == "primary breast tumor",
    "red",
    "blue"
  )
  
  # Plot PCA
  plot(
    pca_result$x[, 1],
    pca_result$x[, 2],
    col = group_colors,
    pch = 19,
    xlab = "Principal Component 1",
    ylab = "Principal Component 2",
    main = "PCA of Breast Cancer Samples"
  )
  
  legend(
    "topright",
    legend = c("Tumour", "Normal"),
    col = c("red", "blue"),
    pch = 19
  ) 
  # Calculate percentage variance explained
  pca_variance <- pca_result$sdev^2 /
    sum(pca_result$sdev^2) * 100
  
  # Display first 10 principal components
  round(pca_variance[1:10], 2)  
  
  # Percentage explained by first two components
  sum(pca_variance[1:2])
  # Load limma for differential expression analysis
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  if (!requireNamespace("limma", quietly = TRUE)) {
    BiocManager::install("limma")
  }
  
  library(limma)
  
  # Prepare log2 expression matrix
  expr_for_limma <- as.matrix(log_expr_matrix)
  
  # Check dimensions and sample groups
  dim(expr_for_limma)
  table(metadata$Group)
  
  # Confirm sample order
  all(colnames(expr_for_limma) == metadata$Sample_ID)
  
  # Define sample groups
  group <- factor(
    metadata$Group,
    levels = c("normal breast tissue",
               "primary breast tumor")
  )
  
  # Create design matrix
  design <- model.matrix(~ 0 + group)
  colnames(design) <- c("Normal", "Tumor")
  
  # View design matrix
  head(design)
  dim(design)

  # Fit linear model
  fit <- lmFit(expr_for_limma, design)
  
  # Compare Tumour against Normal
  contrast_matrix <- makeContrasts(
    Tumor - Normal,
    levels = design
  )
  
  # Apply contrast
  fit2 <- contrasts.fit(fit, contrast_matrix)
  
  # Empirical Bayes moderation
  fit2 <- eBayes(fit2)
  
  # Extract all differential expression results
  deg_results <- topTable(
    fit2,
    coef = 1,
    number = Inf,
    adjust.method = "BH"
  )
  
  # View top 10 results
  head(deg_results, 10)
  
  # Check total results
  nrow(deg_results)

  # Significant probes
  sum(deg_results$adj.P.Val < 0.05, na.rm = TRUE)
  
  # Upregulated probes
  sum(deg_results$adj.P.Val < 0.05 &
        deg_results$logFC > 0, na.rm = TRUE)
  
  # Downregulated probes
  sum(deg_results$adj.P.Val < 0.05 &
        deg_results$logFC < 0, na.rm = TRUE)


  png_plot <- "GSE22820_Volcano_Plot.png"
  
  # Create and save the volcano plot
  png(
    filename = png_plot,
    width = 1200,
    height = 900,
    res = 150
  )
  
  plot(
    deg_results$logFC,
    -log10(pmax(deg_results$adj.P.Val, .Machine$double.xmin)),
    pch = 20,
    cex = 0.5,
    col = ifelse(
      deg_results$adj.P.Val < 0.05 &
        deg_results$logFC > 0, "red",
      ifelse(
        deg_results$adj.P.Val < 0.05 &
          deg_results$logFC < 0, "blue",
        "grey"
      )
    ),
    xlab = "Log2 Fold Change (Tumour vs Normal)",
    ylab = "-Log10 Adjusted P-value",
    main = "Volcano Plot of Breast Cancer Differential Expression"
  )
  
  abline(v = c(-1, 1), lty = 2)
  abline(h = -log10(0.05), lty = 2)
  
  dev.off()
# Map differential expression results to GPL6480 annotation
  
  deg_results$PROBEID <- expression_data$ID_REF[
    as.integer(rownames(deg_results))
  ]
  
  deg_results$SYMBOL <- annotation_data$GENE_SYMBOL[
    match(deg_results$PROBEID, annotation_data$ID)
  ]
  
  deg_results$GENENAME <- annotation_data$GENE_NAME[
    match(deg_results$PROBEID, annotation_data$ID)
  ]
  
  deg_results$ENSEMBL_ID <- annotation_data$ENSEMBL_ID[
    match(deg_results$PROBEID, annotation_data$ID)
  ]
  
  # Save all differential expression results
  write.csv(
    deg_results,
    "GSE22820_All_DEG_Results.csv",
    row.names = FALSE
  )
  
  # Filter significant probes
  significant_DEGs <- subset(
    deg_results,
    adj.P.Val < 0.05
  )
  
  write.csv(
    significant_DEGs,
    "GSE22820_Significant_DEGs.csv",
    row.names = FALSE
  )
  
  # Separate upregulated and downregulated probes
  upregulated <- subset(
    significant_DEGs,
    logFC > 0
  )
  
  downregulated <- subset(
    significant_DEGs,
    logFC < 0
  )
  
  write.csv(upregulated, "GSE22820_Upregulated.csv", row.names = FALSE)
  write.csv(downregulated, "GSE22820_Downregulated.csv", row.names = FALSE)
  # Add probe IDs to differential expression results
  deg_results$PROBEID <- expression_data$ID_REF[
    as.integer(rownames(deg_results))
  ]
  
  # Add GPL6480 gene annotation
  deg_results$GENE_SYMBOL <- annotation_data$GENE_SYMBOL[
    match(deg_results$PROBEID, annotation_data$ID)
  ]
  
  deg_results$GENE_NAME <- annotation_data$GENE_NAME[
    match(deg_results$PROBEID, annotation_data$ID)
  ]
  
  # Keep significant differentially expressed probes
  significant_DEGs <- subset(deg_results, adj.P.Val < 0.05)
  
  # Keep probes with an annotated gene symbol
  mapped_DEGs <- subset(
    significant_DEGs,
    !is.na(GENE_SYMBOL) & GENE_SYMBOL != ""
  )
  
  # Rank genes by adjusted p-value and retain one row per gene
  mapped_DEGs <- mapped_DEGs[
    order(mapped_DEGs$adj.P.Val),
  ]
  
  unique_ranked_genes <- mapped_DEGs[
    !duplicated(mapped_DEGs$GENE_SYMBOL),
  ]
  
  # Save unique ranked genes
  write.csv(
    unique_ranked_genes,
    "GSE22820_Unique_Ranked_Genes.csv",
    row.names = FALSE
  )
  
  # Save top 50 candidate genes
  top50_genes <- head(unique_ranked_genes, 50)
  
  write.csv(
    top50_genes,
    "GSE22820_Top50_Candidate_Genes.csv",
    row.names = FALSE
  )
  # ==========================================
  # DRUG-GENE INTERACTION ANALYSIS USING DGIdb
  # ==========================================
  
  # Load DGIdb interaction results
  drug_data <- read.delim(
    "C:/Users/Jyotsna/Downloads/gene_interaction_results-20_11_2026.tsv",
    check.names = FALSE
  )
  
  # Match top 50 candidate genes with drug interactions
  matched_genes <- merge(
    top50_genes,
    drug_data,
    by.x = "GENE_SYMBOL",
    by.y = "gene"
  )
  
  # Save all matched drug-gene interactions
  write.csv(
    matched_genes,
    "GSE22820_Matched_Drug_Gene_Interactions.csv",
    row.names = FALSE
  )
  
  # Extract interactions labelled as approved
  approved_drugs <- subset(
    matched_genes,
    `regulatory approval` == "Approved"
  )
  
  write.csv(
    approved_drugs,
    "GSE22820_Approved_Drug_Interactions.csv",
    row.names = FALSE
  )
  # ==========================================
  # IDENTIFY POTENTIAL ONCOLOGY DRUG MATCHES
  # ==========================================
  
  # Extract unique approved drug-gene interactions
  approved_summary <- unique(
    approved_drugs[, c(
      "GENE_SYMBOL",
      "drug",
      "indication",
      "interaction score"
    )]
  )
  
  # Identify oncology-related indications
  oncology_matches <- approved_summary[
    grepl(
      "antineoplastic|cancer|neoplasm|oncolog|tumor|tumour",
      approved_summary$indication,
      ignore.case = TRUE
    ),
  ]
  
  # Save potential oncology drug matches
  write.csv(
    oncology_matches,
    "GSE22820_Potential_Oncology_Drug_Matches.csv",
    row.names = FALSE
  )
  
  # Display results
  print(oncology_matches)
    






















































































































