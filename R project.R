# Load the ggplot2 library for data visualization 
library(ggplot2)

# Load the reshape2 library for data manipulation
library(reshape2)

# Read the gNSAF2.csv file into a data frame called 'data'
# The file likely contains protein data with Accession names and values for two groups:
# - Control cancer group
# - Treatment group
data <- read.csv("gNSAF2.csv")

# Remove the first row of the data frame
data <- data.frame(data[-1, ])

# Convert all columns except the first one (Accession names) to numeric
# sapply applies as.numeric(as.character(x)) to each column to ensure values are numeric for analysis
# as.character() prevents issues with factors, and as.numeric() converts strings to numbers
data[, -1] <- sapply(data[, -1], function(x) as.numeric(x))

################ Boxplot for the Original Data #######################
# Purpose: Create boxplots to visualize the distribution of protein expression values
# for each protein across control ('C') and treatment ('T') groups

# Transpose the data frame so rows become columns and vice versa
# This makes each row represent a sample and each column a protein (Accession)
data_transposed <- data.frame(t(data))

# Set column names of the transposed data to the Accession names from the first row
# Assumes the first row of data_transposed contains protein IDs
names(data_transposed) <- data_transposed[1, ]

# Remove the first row (Accession names) since it’s now used as column headers
data_transposed <- data_transposed[-1, ]

# Convert all columns except the first one to numeric for plotting
# as.character() handles factors, and as.numeric() ensures values are suitable for analysis
data_transposed[, -1] <- sapply(data_transposed[, -1], function(x) as.numeric(as.character(x)))

# Add a 'group' column to label samples as control ('C') or treatment ('T')
# The order assumes 8 samples: 2 control, 4 treatment, 2 control (based on provided vector)
data_transposed$group <- c("C", "C", "T", "T", "T", "T", "C", "C")

# Store protein Accession names from column names for use in plotting
proteins <- colnames(data_transposed)

# Example boxplot for one protein (column 4) to compare groups
# aes() maps 'group' to x-axis, protein values to y-axis, and fills boxes by group

ggplot(data_transposed, aes(x = group, y = data_transposed[, 4], fill = group)) +
  geom_boxplot() +  
  labs(x = "Groups", y = "Value") +  # Labels axes for clarity
  ggtitle(proteins[4])  # Titles plot with the protein’s Accession name

# Loop to generate boxplots for the first 5 proteins
# Each plot compares control vs. treatment group distributions
for (i in 1:5) {
  plt <- ggplot(data_transposed, aes(x = group, y = data_transposed[, i], fill = group)) +
    geom_boxplot() +  # Boxplot for the i-th protein
    labs(x = "Groups", y = "Value") +  # Consistent axis labels
    ggtitle(proteins[i])  # Title with the protein’s Accession name
  print(plt)  # Display the plot
}
##########################################################


############# preprocessing #####################

# PQN normalizarion
# Perform PQN to correct for technical variation
# in protein expression data , ensuring samples are comparable
# while preserving biological differences between control ('C') and treatment ('T') groups

data_norm <- data  # Create a copy of the original data frame to store normalized values, preserving 'data' for reference
ref_median <- median(data[, 2], na.rm = TRUE)  
# Calculate the reference median from the second column
# na.rm = TRUE handles missing values (NA)
# Note: Choosing column 2 is arbitrary; it should represent a typical sample (e.g., control)

# Loop through sample columns (2 to ncol(data)), skipping column 1 (protein IDs)
for (i in 2:ncol(data)) {  
  data_norm[, i] <- data[, i] / median(data[, i] / ref_median, na.rm = TRUE)  # Normalize each sample
  # Steps: (1) Compute ratios of sample values to ref_median
  #        (2) Find the median of these ratios (scaling factor)
  #        (3) Divide sample values by this median to align with the reference
  # na.rm = TRUE ensures missing values don’t break the calculation
}

############# Filtration process  #####################


# Divide the normalized dataset into two groups: Control and Treatment
# df1: Control group samples (X10C, X18C, X5C, X8C)
# df2: Treatment group samples (X18T, X19T, X3T, X4T)
df1 <- data.frame(data_norm$Accessions, data_norm$X10C, data_norm$X18C, data_norm$X5C, data_norm$X8C)
df2 <- data.frame(data_norm$Accessions, data_norm$X18T, data_norm$X19T, data_norm$X3T, data_norm$X4T)

# Remove the prefix "data_norm." from the column names for easier readability
names(df1) <- gsub("data_norm.", "", names(df1), fixed = TRUE)
names(df2) <- gsub("data_norm.", "", names(df2), fixed = TRUE)

# Filter out rows with 50% or more missing values (i.e., 3 or more NAs)
# Step 1: Count NAs for each row in the Control group
na_count <- rowSums(is.na(df1))
# Step 2: Keep rows with fewer than 3 NAs
keep_rows <- na_count < 3
df1 <- df1[keep_rows, ]

# Repeat the same NA filtering for the Treatment group
na_count <- rowSums(is.na(df2))
keep_rows <- na_count < 3
df2 <- df2[keep_rows, ]

############# Imputation of missing values  #####################


# Impute missing values in the Control group (df1) by the median of each row
# Loop through each row of df1 by row names
for(i in rownames(df1)) {
  
  # Calculate the median of the current row (excluding the first column: Accessions)
  row_median <- median(as.numeric(df1[i, -1]), na.rm = TRUE)
  
  # Loop through each column in the row
  for(j in colnames(df1)) {
    
    # If the cell is NA, replace it with the row median
    if(is.na(df1[i, j])) {
      df1[i, j] <- row_median
    }
    
  }
}

# Repeat the imputation process for the Treatment group (df2)
for(i in rownames(df2)) {
  
  row_median <- median(as.numeric(df2[i, -1]), na.rm = TRUE)
  
  for(j in colnames(df2)) {
    
    if(is.na(df2[i, j])) {
      df2[i, j] <- row_median
    }
    
  }
}


#########################################################################
#using joins
library(dplyr)
# full join to show unique protiens that appear in one of the groups
full_join_df <- full_join(df1, df2, by = "Accessions") 

#inner join to take common protiens 
inner_join_df <-inner_join(df1, df2, by = "Accessions")

####################Statistical Analysis###############################################

rows_mean <- rowMeans(inner_join_df[,-1])
shapiro.test(rows_mean) # p < 0.05   so it is non parametric 

#scale data (optinal)
#inner_join_df[,-1] = scale(inner_join_df[,-1])


# --- Non-paired Wilcoxon test for differential expression ---

# 'ptable' will store the results: Accessions, raw p-values, adjusted p-values, and fold changes
ptable <- data.frame(inner_join_df[,1])  # Start by copying the Accessions column
colnames(ptable)[1] = "Accessions"       # Rename the first column for clarity

# Calculate mean expression for each protein across Control and Treatment samples
group_c <- apply(inner_join_df[,2:5], 1, mean)  # Columns 2 to 5 are Control samples
group_T <- apply(inner_join_df[,6:9], 1, mean)  # Columns 6 to 9 are Treatment samples

# Compute fold change as Treatment mean / Control mean for each protein
foldchange <- group_T / group_c

# Loop through each row (protein) to perform statistical testing and build the ptable
for (i in rownames(inner_join_df)) {
  
  # Perform a non-paired Wilcoxon rank-sum test between Control and Treatment groups
  ptable[i, 2] <- wilcox.test(
    as.numeric(inner_join_df[i, 2:5]),
    as.numeric(inner_join_df[i, 6:9]),
    paired = FALSE,
    exact = FALSE
  )$p.value
  
  # Adjust the raw p-value using Benjamini-Hochberg (FDR control)
  ptable[i, 3] <- p.adjust(ptable[i, 2], method = "BH")
  
  # Add fold change value for this protein
  ptable[i, 4] <- foldchange[as.numeric(i)]
}

# Rename columns for better readability
colnames(ptable)[2] = "pvalue"
colnames(ptable)[3] = "p-adj"
colnames(ptable)[4] = "FC"  # Fold Change

###################################################################################

############# Volcano Plot to Visualize Significant Proteins ###################

# Create basic volcano plot: log2(FC) vs. -log10(p-adjusted)
p <- ggplot(ptable, aes(x = log2(FC), y = -log10(`p-adj`))) +
  geom_point() +  # plot points for each protein
  labs(x = "Fold change (log2)", y = "-log10(p-adj)") +
  theme_minimal()  # clean background

# Add red threshold lines:
# Vertical lines at log2(FC) = ±1 (2-fold change up/down)
# Horizontal line at -log10(0.05) = significance threshold
p2 <- p + 
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red")

# Add a new column 'diffexpressed' to classify significance
ptable$diffexpressed <- "NO"  # default: not significantly differentially expressed

# Mark as "UP" if FC > 1  and p-adj < 0.05
ptable$diffexpressed[log2(ptable$FC) > 1 & ptable$`p-adj` < 0.05] <- "UP"

# Mark as "DOWN" if FC < -1  and p-adj < 0.05
ptable$diffexpressed[log2(ptable$FC) < -1 & ptable$`p-adj` < 0.05] <- "DOWN"

# Plot with colored points according to significance category
p <- ggplot(data = ptable, aes(x = log2(FC), y = -log10(`p-adj`), col = diffexpressed)) +
  geom_point() + 
  theme_minimal()

# Add threshold lines again to this new colored plot
p2 <- p + 
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red")

# Assign manual colors: blue for DOWN, red for UP, black for NO change
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)

# Optional: Add protein labels only for significant ones
ptable$ptablelabel <- NA
ptable$ptablelabel[ptable$diffexpressed != "NO"] <- ptable$Accessions[ptable$diffexpressed != "NO"]

# Plot with labels (without overlap protection)
ggplot(data = ptable, aes(x = log2(FC), y = -log10(`p-adj`), col = diffexpressed, label = ptablelabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()  # can be messy if many points are labeled

# Use ggrepel for better label placement (prevents overlapping)
library(ggrepel)
ggplot(data = ptable, aes(x = log2(FC), y = -log10(`p-adj`), col = diffexpressed, label = ptablelabel)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = mycolors) +  # use manual color scheme
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  labs(x = "Fold change (log2)", y = "-log10(p-adj)") +
  geom_text_repel()  # smart text label placement

#######################################################################################

######### Heatmap for Significant Proteins per Group ####################

# Filter only the significant proteins (adjusted p-value < 0.05)
significant_proteins <- ptable[ptable$`p-adj` < 0.05,]

# Extract corresponding expression data from the main data frame
sig_data <- inner_join_df[inner_join_df$Accessions %in% significant_proteins$Accessions,]

# Load the pheatmap library for plotting
library(pheatmap)

# Transpose the data frame so that:
# - Rows become samples
# - Columns become proteins (required format for pheatmap)
data_trans <- data.frame(t(sig_data))

# Rename columns using the first row (protein Accessions)
names(data_trans) <- data_trans[1,]  # Set column names from the first row
data_trans <- data_trans[-1,]        # Remove the first row (was used for naming)

# Convert all data to numeric (as pheatmap expects numeric matrix)
data_trans[] <- lapply(data_trans, function(x) as.numeric(as.character(x)))

# Optional: Scale the data (z-score normalization) across rows (i.e., across samples)
data_trans = scale(data_trans)

# Plot the heatmap
pheatmap(data_trans, 
         main = "Significant Proteins Heatmap")  # Title for the heatmap

############################################################

############### Correlation Between Control and Treatment Groups ###############

# Calculate the column-wise means (i.e., average expression per sample) for the control group (columns 2-5)
group_C_colmeans <- colMeans(inner_join_df[,2:5])
# Calculate the column-wise means for the treatment group (columns 6-9)
group_T_colmeans <- colMeans(inner_join_df[,6:9])

# Perform Pearson correlation test between the control and treatment group column means
# This evaluates the similarity of average expression across the groups at the sample level
res <- cor.test(group_C_colmeans, group_T_colmeans, method = "pearson")
# View the correlation test result (output includes correlation coefficient, p-value, etc.)
res

# Calculate the row-wise means (i.e., average expression per protein) for the control group (columns 2-5)
group_C_rowmeans <- rowMeans(inner_join_df[,2:5])
# Calculate the row-wise means for the treatment group (columns 6-9)
group_T_rowmeans <- rowMeans(inner_join_df[,6:9])
# Perform Pearson correlation test between the control and treatment group row means
# This evaluates how similar the average expression of each protein is across the two groups
cor.test(group_C_rowmeans, group_T_rowmeans, method = "pearson")

################### biological analysis #############################

# Load the dataset containing information on biological pathways and proteins
# The encoding is set to "UTF-8" to avoid issues with special characters in the file.
data2 <- read.csv("gProfiler_hsapiens.csv", fileEncoding = "UTF-8")

# Filter significant entries based on the adjusted p-value < 0.05
# This will select only the pathways that have significant associations (after correction for multiple testing)
significant <- data2[data2$adjusted_p_value < 0.05,]

# Extract only the Gene Ontology (GO) pathways from the dataset
# The 'source' column contains the type of pathway, and we are interested in those that start with "GO"
GO_pathes <- data2[grepl('^GO', data2$source), ]

# Split the proteins listed in the 'intersections' column for each GO pathway into individual proteins
# The 'intersections' column contains a comma-separated list of proteins involved in each pathway
proteins_list <- strsplit(GO_pathes$intersections, ",")

# Create a data frame that associates each protein with its corresponding GO pathway (term_id)
# The 'term_id' is replicated for each protein listed under the GO pathway in the 'intersections' column
pathway_proteins <- data.frame(
  term_id = rep(GO_pathes$term_id, lengths(proteins_list)), # Repeat term_id for each protein in the pathway
  proteins = unlist(proteins_list) # Unlist the proteins from the split strings into a single vector
)

# Diagram 1: Visualization of proteins in each GO pathway

# ggplot function to create a scatter plot
ggplot(pathway_proteins, aes(x = term_id, y = proteins)) + 
  # Use geom_point() to create a scatter plot showing proteins on the y-axis and pathways (term_id) on the x-axis
  geom_point() + 
  
  # Set the x-axis label to "Pathways"
  xlab("Pathways") + 
  
  # Set the y-axis label to "Proteins"
  ylab("Proteins") + 
  
  # Add a title to the plot
  ggtitle("Proteins in Each Pathway") + 
  
  # Adjust the x-axis labels to make them readable by rotating the text by 90 degrees
  # This is useful when pathway names (term_id) are long and may overlap
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



# Diagram 2: Sankey diagram to visualize the relationship between proteins and GO pathways

# Load the necessary package for Sankey diagram
library(ggsankey)

# Prepare the data for the Sankey diagram by converting it into long format
# The 'make_long' function creates two columns: one for proteins and one for pathways (term_id)
data_pathway <- pathway_proteins %>%
  make_long(proteins, 
            term_id)

# Create the Sankey diagram with ggplot
ggplot(tail(data_pathway, n = 60), 
       aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
  
  # Use geom_sankey to create the flow visualization with some transparency (flow.alpha)
  geom_sankey(flow.alpha = 0.75, node.color = 1, space = 70) +
  
  # Add labels to the nodes using geom_sankey_label and adjust label size and color
  geom_sankey_label(size = 6, color = "white", space = 70) +
  
  # Use a color scale for the nodes with the viridis color palette
  scale_fill_viridis_d(option = "H", alpha = 0.95) +
  
  # Apply the theme for the Sankey diagram to adjust base text size
  theme_sankey(base_size = 20) +
  
  # Remove the x-axis label and adjust the title
  labs(x = NULL) +
  
  # Customize the plot appearance by removing the legend and centering the title
  theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
  
  # Add a title to the plot
  ggtitle("Pathways of Proteins")

# Diagram 3: Sankey diagram using networkD3 to visualize the relationship between proteins and GO pathways

# Load the networkD3 package, which is used for creating interactive Sankey diagrams
library(networkD3)

# Create a 'links' data frame to represent the connections between proteins and pathways
# Each row in 'links' represents a relationship between a protein and a GO pathway
# The 'source' column corresponds to proteins, and the 'target' column corresponds to GO pathways (term_id)
links <- data.frame(source = paste0(pathway_proteins$proteins),
                    target   = paste0(pathway_proteins$term_id))

# Convert the 'source' and 'target' columns into character data type for consistency
# This is required for proper matching of node names when building the Sankey diagram
links$source <- as.character(links$source)
links$target <- as.character(links$target)

# Create a 'nodes' data frame, which contains all unique proteins and pathways (source and target combined)
# The 'name' column contains the node identifiers (protein names and pathway term IDs)
nodes <- data.frame(name = unique(c(links$source, links$target)))

# Update the 'source' and 'target' columns in the 'links' data frame to match the indices in the 'nodes' data frame
# This step ensures that the 'source' and 'target' columns contain node indices starting from 0
links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1

# Add a 'value' column to the 'links' data frame, where each connection between a protein and a pathway has a value of 1
# This indicates that there is a direct connection, though the value could be adjusted based on other criteria (e.g., counts or weights)
links$value <- 1 

# Use the sankeyNetwork function from the networkD3 package to generate the Sankey diagram
# 'Links' refers to the data frame containing the connections between nodes (proteins and pathways),
# 'Nodes' is the data frame of unique node identifiers, and the other parameters specify which columns to use for source, target, and values
# Additional styling parameters like 'fontSize' and 'nodeWidth' adjust the appearance of the diagram
graph <- sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                       Target = 'target', Value = 'value', NodeID = 'name', fontSize = 12, nodeWidth = 30)

# Display the Sankey diagram in the output
graph

