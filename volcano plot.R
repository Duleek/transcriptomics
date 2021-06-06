readdata <- read.csv("genes.csv", head = TRUE)
class(readdata)
data <- read.csv("genes.csv", head = TRUE)
de <- data[complete.cases(data), ]

#create a new column of data to track which genes are upregulated and downregulated
#set thresholds for what is considered up/down regulated based on log fold change in expression and p-value
#in this case, genes are significant if they have double or half the control's expression (base 2 log of 1 and -1 respectively)
#p must be < 0.05 to be significant
#if adjusted p values are available, use those instead (much more accurate)

de$diffexpressed <- "NO"
de$diffexpressed[de$logFC > 1 & de$P.Value <0.05] <- "UP"
de$diffexpressed[de$logFC < -1  & de$P.Value <0.05] <- "DOWN"


#plot all points (each point is a gene)
#log 2 fold change in the x axis, -log10(p value) on the y axis

p <- ggplot(data=de, aes(x=de$logFC, y=-log10(de$P.Value), col=diffexpressed)) + geom_point() + theme_minimal()
p

#add vertical line for log 2 fold change threshold
#add horizontal line for p threshold: p = 0.05 so -log10(p) = 1.3
p2 <- p + geom_vline(xintercept=c(-1, 1), col="black") +
    geom_hline(yintercept=-log10(0.05), col="black")
p2

#change colour based on category
mycolors <- c("red", "green", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

#add labels to significant genes
de$delabel <- NA
de$gene <- as.character(de$gene)
de$delabel[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]

ggplot(data=de, aes(x=de$logFC, y=-log10(de$P.Value), col=diffexpressed, label=delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text()

#clean up organization of gene labels and make a final plot

library(ggrepel)

p4 <- ggplot(data=de, aes(x=de$logFC, y=-log10(de$P.Value), col=diffexpressed, label=delabel)) +
        geom_point() + 
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("red", "black", "green")) +
        geom_vline(xintercept=c(-1, 1), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black")

p4 + labs (x = "Log 2-fold change", y = " -log10(p-value)")
