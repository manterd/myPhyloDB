require(phyloseq)
require(ape)
require(ggplot2)
require(picante)
require(FD)
require(reshape2)

#######################################################


# subset data from phyloseq file using assigned metadata, then run anova as normal
# use Anova data to create standard R plots, as well as a dataframe to hand back to python
# back in python, use dataframe to create highcharts graphs as per normal

print(getwd())

biomLoc <- paste(pathToBiom,'phyloseq.biom', sep='/')
print('Hello, world!')

P <- import_biom(biomLoc)
print('Imported!')

names(tax_table(P))
print('Taxa!')

#colnames(tax_table(P)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species","OTU99")
print('Columns!')

dat <- sample_data(P)
print('Samples!')

#fullTree <- read.tree("gg_13_5_99.tree")
print('Need!')

taxa <- taxa_names(P)
print('Better!')

#tree <- drop.tip(fullTree,fullTree$tip.label[-match(taxa,fullTree$tip.label)])
print('Debug!')

#P <- merge_phyloseq(P,tree)
print('Tools!')
