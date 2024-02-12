#read all the libraries
library(readxl)
library(stringr)
library(dplyr)
library(tableone)

#read the file containing information regrading different drug class
drug_list <- read_excel("list.xlsx",sheet="drug_list")

#Store the exposed and control drug list in variable drug_exp and drug_cntrl
# For Sitagliptin vs (Metformin,glimepiride,pioglitazone,glipizide) which are store in Exposed5 and Control3,Control4,Control5,Control6
drug_exp5 <- na.omit(tolower(drug_list$Exposed_5))
drug_cntrl3 <- na.omit(tolower(drug_list$control_3))
drug_cntrl4 <- na.omit(tolower(drug_list$control_4))
drug_cntrl5 <- na.omit(tolower(drug_list$control_5))
drug_cntrl6 <- na.omit(tolower(drug_list$control_6))

# removes the attributes from the data 
attributes(drug_exp5) <- NULL
attributes(drug_cntrl3) <- NULL
attributes(drug_cntrl4) <- NULL
attributes(drug_cntrl5) <- NULL
attributes(drug_cntrl6) <- NULL

#Read file containing information regarding primary adverse event
adverse_rec <- read_excel("list.xlsx",sheet="adverse_event")
adverse_rec1 <- na.omit(tolower(adverse_rec$Reaction_1)) #keywords to identify liver metastasis
adverse_rec2 <- na.omit(tolower(adverse_rec$Reaction_2)) #keywords to identify pancreatic cancer
adverse_rec3 <- na.omit(tolower(adverse_rec$Reaction_3)) #keywords to identify colorectal cancer
adverse_rec4 <- na.omit(tolower(adverse_rec$Reaction_4)) #keywords to identify gastric cancer

# removes the attributes from the data 
attributes(adverse_rec1) <- NULL
attributes(adverse_rec2) <- NULL
attributes(adverse_rec3) <- NULL
attributes(adverse_rec4) <- NULL

# adverse event with both liver metastasis and pancreatic/colorectal/gastric cancer
panc_rec <- c(adverse_rec1,adverse_rec2)
col_rec <- c(adverse_rec1,adverse_rec3)
gas_rec <- c(adverse_rec1,adverse_rec4)

#read the Exposed data downloaded from FAERs public dashboard  for - year 2012-2022, age group 18-85

#Sitagliptin Users
exposed5 <- read_excel("exposed5.xlsx")
#duplicate  records delete
stglptn <- exposed5[!duplicated(exposed5),]

#read the Control data downloaded from FAERs public dashboard - year 2012-2022, age group 18-85

#Metformin
control3 <- read_excel("control3.xlsx")
metformin <- control3[!duplicated(control3),]      #duplicate  records delete

#Glimepiride
control4 <- read_excel("control4.xlsx")
glimepiride <- control4[!duplicated(control4),]      #duplicate  records delete

#Pioglitazone
control5 <- read_excel("control5.xlsx")
pioglitazone <- control5[!duplicated(control5),]      #duplicate  records delete

#glipizide
control6 <- read_excel("control6.xlsx")
glipizide <- control6[!duplicated(control6),]      #duplicate  records delete

#Function to count the number of cases for both exposed and control group

count_exp <- function(dataset_list,concmed)
{
  # Initialize an empty list to store the results
  result <- list()
  
  # Iterate through the datasets 
  for (i in seq_along(dataset_list))
  {
    dataset <- dataset_list[[i]]
    
    #Remove all the cases with compare drug as among suspected product
    for(j in seq_along(concmed))
    {
      drug<-concmed[[j]]
      #Remove cases having suspect product(concomitant medication) from compare drug
      dataset$`Suspect Product Active Ingredients` <- tolower(dataset$`Suspect Product Active Ingredients`)
      dataset$presence <- ifelse(rowSums(sapply(drug, function(word) grepl(word,dataset$`Suspect Product Active Ingredients`))) > 0, 1, 0)
      dataset <- dataset[dataset$presence==0,-c(25)]
      
      #"Reason for use" is a column which store the names of diseases for which primary suspect drug has been used 
      #Remove all cases which have liver metastasis or pancreatic/colorectal/gastric cancer disease mention in "Reason for use"
      dataset$`Reason for Use` <- tolower(dataset$`Reason for Use`)
      
      dataset$panc <- ifelse(rowSums(sapply(panc_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for pancreatic cancer
      dataset_panc <- dataset[dataset$panc==0,-c(25)]
      
      dataset$col <- ifelse(rowSums(sapply(col_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for colorectal cancer
      dataset_col <- dataset[dataset$col==0,-c(25)]
      
      dataset$gas <- ifelse(rowSums(sapply(gas_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for gastric cancer
      dataset_gas <- dataset[dataset$gas==0,-c(25)]
      
      #consider only records with either both liver metastasis pancreatic/colorectal/gastric cancer adverse event or no reaction
      dataset_panc$Reactions <- tolower(dataset_panc$Reactions)
      dataset_panc$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_panc$Reactions))) > 0, 1, 0)
      dataset_panc$presence2 <- ifelse(rowSums(sapply(adverse_rec2, function(word) grepl(word,dataset_panc$Reactions))) > 0, 1, 0)
      
      dataset_col$Reactions <- tolower(dataset_col$Reactions)
      dataset_col$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_col$Reactions))) > 0, 1, 0)
      dataset_col$presence2 <- ifelse(rowSums(sapply(adverse_rec3, function(word) grepl(word,dataset_col$Reactions))) > 0, 1, 0)
      
      dataset_gas$Reactions <- tolower(dataset_gas$Reactions)
      dataset_gas$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_gas$Reactions))) > 0, 1, 0)
      dataset_gas$presence2 <- ifelse(rowSums(sapply(adverse_rec4, function(word) grepl(word,dataset_gas$Reactions))) > 0, 1, 0)
      
      #Count the number of cases with or without adverse event for each of the categories
      n1 <- nrow(dataset_panc[dataset_panc$presence1==1 & dataset_panc$presence2==1,])
      n2 <- nrow(dataset_panc[!(dataset_panc$presence1==1 & dataset_panc$presence2==1),])
      
      n3 <- nrow(dataset_col[dataset_col$presence1==1 & dataset_col$presence2==1,])
      n4 <- nrow(dataset_col[!(dataset_col$presence1==1 & dataset_col$presence2==1),])
      
      n5 <- nrow(dataset_gas[dataset_gas$presence1==1 & dataset_gas$presence2==1,])
      n6 <- nrow(dataset_gas[!(dataset_gas$presence1==1 & dataset_gas$presence2==1),])
      
      comb_nam <- paste0(name1[i],"_",name2[j])
      result[[comb_nam]] <- list(exp_pan=n1,exp_no_pan=n2,exp_col=n3,exp_no_col=n4,exp_gas=n5,exp_no_gas=n6)
    }
  }
  
  # Return the results as a named list
  return(result)
}

#***************Pass arguments****************

# For exposed dataset
dataset <- list(stglptn)    #Pass all the exposed dataset names
concmed <- list(drug_cntrl3,drug_cntrl4,drug_cntrl5,drug_cntrl6)    #Pass all the concomitant  medication list
#Pass all the names of Exposed and Control drug
name1 <- c("sitagliptin")
name2 <- c("metformin","glimepiride","pioglitazone","glipizide")
result1 <- count_exp(dataset,concmed)   #Call the fucntion
#Get a table from the list
max_length <- max(sapply(result1, length))
count_ae1 <- do.call(rbind, lapply(result1, function(x) c(x, rep(NA, max_length - length(x)))))
count_ae1 <- as.data.frame(count_ae1)

#Count number of cases for control groups   
count_cntrl <- function(dataset_list,concmed)
{
  # Initialize an empty list to store the results
  result <- list()
  
  # Iterate through the dataset 
  for (i in seq_along(dataset_list))
  {
    dataset <- dataset_list[[i]]
    
    #Remove all the cases with compare drug as among suspected product
    for(j in seq_along(concmed))
    {
      drug<-concmed[[j]]
      #Remove cases having suspect product(concomitant medication) from compare drug
      dataset$`Suspect Product Active Ingredients` <- tolower(dataset$`Suspect Product Active Ingredients`)
      dataset$presence <- ifelse(rowSums(sapply(drug, function(word) grepl(word,dataset$`Suspect Product Active Ingredients`))) > 0, 1, 0)
      dataset <- dataset[dataset$presence==0,-c(25)]
      
      #"Reason for use" is a column which store the names of diseases for which primary suspect drug has been used 
      #Remove all cases which have liver metastasis or pancreatic/colorectal/gastric cancer disease mention in "Reason for use"
      dataset$`Reason for Use` <- tolower(dataset$`Reason for Use`)
      
      dataset$panc <- ifelse(rowSums(sapply(panc_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for pancreatic cancer
      dataset_panc <- dataset[dataset$panc==0,-c(25)]
      
      dataset$col <- ifelse(rowSums(sapply(col_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for colorectal cancer
      dataset_col <- dataset[dataset$col==0,-c(25)]
      
      dataset$gas <- ifelse(rowSums(sapply(gas_rec, function(word) grepl(word,dataset$`Reason for Use`))) > 0, 1, 0) # for gastric cancer
      dataset_gas <- dataset[dataset$gas==0,-c(25)]
      
      #consider only records with either both liver metastasis pancreatic/colorectal/gastric cancer adverse event or no reaction
      dataset_panc$Reactions <- tolower(dataset_panc$Reactions)
      dataset_panc$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_panc$Reactions))) > 0, 1, 0)
      dataset_panc$presence2 <- ifelse(rowSums(sapply(adverse_rec2, function(word) grepl(word,dataset_panc$Reactions))) > 0, 1, 0)
      
      dataset_col$Reactions <- tolower(dataset_col$Reactions)
      dataset_col$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_col$Reactions))) > 0, 1, 0)
      dataset_col$presence2 <- ifelse(rowSums(sapply(adverse_rec3, function(word) grepl(word,dataset_col$Reactions))) > 0, 1, 0)
      
      dataset_gas$Reactions <- tolower(dataset_gas$Reactions)
      dataset_gas$presence1 <- ifelse(rowSums(sapply(adverse_rec1, function(word) grepl(word,dataset_gas$Reactions))) > 0, 1, 0)
      dataset_gas$presence2 <- ifelse(rowSums(sapply(adverse_rec4, function(word) grepl(word,dataset_gas$Reactions))) > 0, 1, 0)
      
      #Count the number of cases with or without adverse event for each of the categories
      n1 <- nrow(dataset_panc[dataset_panc$presence1==1 & dataset_panc$presence2==1,])
      n2 <- nrow(dataset_panc[!(dataset_panc$presence1==1 & dataset_panc$presence2==1),])
      
      n3 <- nrow(dataset_col[dataset_col$presence1==1 & dataset_col$presence2==1,])
      n4 <- nrow(dataset_col[!(dataset_col$presence1==1 & dataset_col$presence2==1),])
      
      n5 <- nrow(dataset_gas[dataset_gas$presence1==1 & dataset_gas$presence2==1,])
      n6 <- nrow(dataset_gas[!(dataset_gas$presence1==1 & dataset_gas$presence2==1),])
      
      comb_nam <- paste0(name2[j],"_",name1[i])
      result[[comb_nam]] <- list(cntr_pan=n1,cntr_no_pan=n2,cntr_col=n3,cntr_no_col=n4,cntr_gas=n5,cntr_no_gas=n6)
    }
  }
  
  # Return the results as a named list
  return(result)
}

# For control dataset
dataset <- list(metformin,glimepiride,pioglitazone,glipizide) #Pass all the exposed dataset names
concmed <- list(drug_exp5) #Concomitant  medication list
#Names of Exposed and Control drug
name1 <- c("metformin","glimepiride","pioglitazone","glipizide")
name2 <- c("sitagliptin")
result2 <- count_cntrl(dataset,concmed)   #Call the fucntion
#Get a table from the list
max_length <- max(sapply(result2, length))
count_ae2 <- do.call(rbind, lapply(result2, function(x) c(x, rep(NA, max_length - length(x)))))
count_ae2 <- as.data.frame(count_ae2)

#Combine both the Exposed and control dataframe into one 
count_ae_pan <-  cbind(count_ae1[,c(1,2)], count_ae2[,c(1,2)])
count_ae_col <-  cbind(count_ae1[,c(3,4)], count_ae2[,c(3,4)])
count_ae_gas <-  cbind(count_ae1[,c(5,6)], count_ae2[,c(5,6)])

#Derive ROR and CI values
#install.packages("epitools")
library(epitools)

drug <- c("Interested","Control")
outcome <- c("Liver_mets","No_Liver_mets")
p_val <- NULL
name <- NULL
ROR <- NULL
CI_low <- NULL
CI_high <- NULL

info_tab <- function(dataset)
{
  for(i in 1:nrow(dataset))
  {
    if (dataset[i, 1] == 0 || dataset[i, 3] == 0) {
      next  # Skip the current iteration and move to the next one
    }
    data <- matrix(c(dataset[[i,1]],dataset[[i,2]],dataset[[i,3]],dataset[[i,4]]),nrow=2, ncol=2, byrow=TRUE)
    dimnames(data) <- list('Outcome'=outcome,'Drug'=drug)
    info <- oddsratio(data)
    name[i] <- drug_comb[i]
    p_val[i] <- round(ifelse(dataset[[i,1]] < 5 || dataset[[i,3]] < 5, info$p.value[2,2], info$p.value[2,3]),5)
    ROR[i] <- round(info$measure[2,1],3)
    CI_low[i] <- round(info$measure[2,2],3)
    CI_high[i] <- round(info$measure[2,3],3)
  }
  df <- as.data.frame(cbind(name,ROR,p_val,CI_low,CI_high))
  return(df)
}
drug_comb <- c("sitagliptin_metformin","sitagliptin_glimepiride","sitagliptin_pioglitazone","sitagliptin_glipizide")
info_tab_gas <- info_tab(count_ae_gas)    #pass arguments one-by-one