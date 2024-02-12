##Reported Odds ratio calculation
#########
#function
#######
# d is a data frame with 4 columns
# d$x gives variable names
# d$y gives center point
# d$ylo gives lower limits
# d$yhi gives upper limits
forestplot <- function(d, xlab="Odds Ratio", ylab="Study"){
  require(ggplot2)
  p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi, fill=x)) + 
    geom_pointrange(width=0.1, color="blue", fill="blue") + 
    coord_flip() +
    geom_hline(aes(x=1, yintercept= 1), lty=2) +
    theme_bw()+
    ylab(xlab) +
    xlab(ylab)#switch because of the coord_flip() above
  return(p)
}

file= read_xlsx("~/Desktop/Liver_sc_Atlas/Figures/Drug_effect_liver_injury.xlsx", sheet = 2)
file= data.frame(file)

colnames(file) = c( "x", "ylo", "yhi","y") ### x= name of the feature, ylo=min CI value, yhi= max CI value, and y=ROR value
#png("~/Desktop/ECMO_Vs_MODS_data/OR_for_key_genes.png")
forestplot(file)+theme(axis.text = element_text(size=14, face="bold"), panel.border = element_rect(fill=NA, colour = "black", size=1))
ggsave("~/Desktop/Ke_metastatic_work/Figures/Drug_effect_liver_injury.pdf", width = 7, height = 4)