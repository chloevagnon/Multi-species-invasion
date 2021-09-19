#####################################################################################################################
#                                          R code example for :                                                     #
#  Inferring the trophic attributes of new co-occurring invaders in the largest natural French lake (Lake Bourget)  #
#                                                                                                                   #
# Author: Chloe Vagnon                                                                                              #
# Date: September 2021                                                                                              #
#                                                                                                                   #
# This scipt provides an example for using the different functions of the aNM applied for inferring the trophic     #
# attributes of consummers.                                                                                         #
# This script is seperated into three parts divided in different steps:                                             #
#   1) Inferring and refining all the trophic links                                                                 #
#      A) The prey body sizes range are reconstructed for consumers based on their category and their body size     #
#      B) All the possible trophic links are inferred using QRs calibrated for the aNM (Vagnon et al. 2021)         #
#      C) The links are refined based on the consumer diet type (e.g., carnivorous)                                 #
#      D) The resulting links are refined based on the consumer habitat (e.g., pelagic)                             #
#      E) The final links are weigthed                                                                              #
#                                                                                                                   #
#   2) Example of data analyses considering one species                                                             # 
#      A) Generate binary matrices with the Bernouilli trials from the weighted links                               #
#      B) Consider the 4th species as a predator in the food web                                                    #
#      C) Consider the 4th species as a prey in the food web                                                        #
#                                                                                                                   #
#   3) Example of data analyses considering a gradient of body size for a small invasive fish                       #
#      A) Infer the diet range for a body size gradient of invaders                                                 #
#      B) Get information on prey                                                                                   #
#      C) Other methods to analyze information on prey                                                              #
#                                                                                                                   #
# Literature cited:                                                                                                 #
# Vagnon C, Cattanéo F, Goulon C, Grimardias D, Guillard J, Frossard V (2021) An allometric niche model for species # 
#       interactions in temperate freshwater ecosystems. Ecosphere 12:e03420                                        #
#                                                                                                                   #
#####################################################################################################################


###### 1) Inferrence and refinement of the trophic links among species in an inventory #####

#Loading library
library(stringr)
library(ggplot2)
library(Rlab)
library(dplyr)
library(data.table)
library(vegan)
library(FactoMineR)
library(factoextra)

# Functions used to infer trophic links
source("cv_Functions_aNM_Multispecies.R")

#Parameters to calculate the consumers diet range (regressions quantiles)
load("Param_reginvert.Rdata") # For invertebrates
load("Param_regvert.Rdata")   # For vertebrates

#Loading of the species inventory created for the example
# /!\ Species have to be ranged according by decreasing body size
load("Example_SpInventory.Rdata")


#### A) The resource body sizes range are reconstructed for consumers based on : 
#1. species category (i.e., primary producers, zooplankton,invertebrate, vertebrate)
#2. body size 

#Note that rows corresponding to primary producers are automatically filled with 0.
bs_r<-get_niche_attributes(species_name=DATA$Species,
                           body_size = DATA$Log10bs,
                           species_category = DATA$Category)


#### B) Inferrences of all trophic links (binary matrix)
Bmat<-L_fn2(name=bs_r$name,n=bs_r$n,c=bs_r$c,low=bs_r$low,high=bs_r$high,table="NO")


#### C) Links refinement based on consumers" diet
Bmat_Diet<-Ref_L_Diet(Bmat=Bmat, diet=DATA$Diet,Table="NO")


#### D) Links refinement based on consumers" habitat
Bmat_Hab<-Ref_L_Hab(Bmat=Bmat_Diet,habitat =DATA$Habitat,Table="NO")
# Note that supplementary refinement can be implemented depending on the ecosystems or the
# species studied.In this example fish and predatory invertebrates are not allowed to consume 
# primary producers.
# Final binary matrix Mb is given by: 
Mb<-Bmat_Hab


#### E) Weighting of refined trophic links
Mb_W<-Weighting(Niche_attributes = bs_r,Bmat=Mb)



##### 2) Example of data analyses considering one species in analyses #####

####A) Generate the 100 matrices from the Bernouilli trials 
BernMat<-make_bern(n=100,Wmat=Mb_W)


##### B) Consider the 4th species as a predator in the food web (one species)
BernSP4_Pred<-data.frame(Simu=NA,Prey=NA)
for (i in 1:100){
  vecPrey<-which(BernMat[i,,which(colnames(Mb_W)=="SP4")]!=0)
  if (length(vecPrey)!=0){
    Preyname<-data.frame(Simu=i,Prey=rownames(Mb_W)[vecPrey])
    
  }else{Preyname<-data.frame(Simu=i,Prey=NA)}
  
  BernSP4_Pred<-rbind(BernSP4_Pred,Preyname)  
  
}
BernSP4_Pred<-na.omit(BernSP4_Pred)
BernSP4_Pred$Habitat<-DATA$Habitat[match(BernSP4_Pred$Prey,DATA$Species)]
BernSP4_Pred$Category<-DATA$Category[match(BernSP4_Pred$Prey,DATA$Species)]
BernSP4_Pred$PreySize<-DATA$bs[match(BernSP4_Pred$Prey,DATA$Species)]
BernSP4_Pred$Log10PreySize<-log10(DATA$bs[match(BernSP4_Pred$Prey,DATA$Species)]) 

# Plot of the prey body size distribution accouting for all the simulations
ggplot(data=BernSP4_Pred,aes(x=Log10PreySize),breaks=seq(0, 4, by=0.2))+
  geom_histogram()+
  scale_x_continuous("Log10(Prey body size, µm)",breaks = seq(0, 4, by=0.2))+
  ylab("Counts")+
  theme(axis.text = element_text(size=16,colour = "black"),
        axis.title = element_text(size=18,colour = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_rect(colour = "black", fill=NA,size=1),
        panel.grid.major = element_line(size = 0.5, linetype = "solid",colour = "grey86"),
        legend.background = element_rect(fill = "transparent", colour = NA))


# Prey summarized by habitat
Prey_Habitat<-BernSP4_Pred%>%group_by(Simu)%>%summarise(data.frame(table(Habitat)))

#Mean number of prey in each habitat and plot
Mean_Hab<-data.frame(Habitat=unique(Prey_Habitat$Habitat),
                     Mean=tapply(Prey_Habitat$Freq, Prey_Habitat$Habitat, mean),
                     Sd=tapply(Prey_Habitat$Freq, Prey_Habitat$Habitat, sd))
ggplot(data = Mean_Hab)+
  geom_bar(aes(x=Habitat,y=Mean),stat="identity")+
  theme(axis.text = element_text(size=16,colour = "black"),
        axis.title = element_text(size=18,colour = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_rect(colour = "black", fill=NA,size=1),
        panel.grid.major = element_line(size = 0.5, linetype = "solid",colour = "grey86"),
        legend.background = element_rect(fill = "transparent", colour = NA))

# Prey summarized by category
Prey_Category<-BernSP4_Pred%>%group_by(Simu)%>%summarise(data.frame(table(Category)))

#Mean number of prey in each category and plot
Mean_Cat<-data.frame(Category=unique(Prey_Category$Category),
                     Mean=tapply(Prey_Category$Freq, Prey_Category$Category, mean),
                     Sd=tapply(Prey_Category$Freq, Prey_Category$Category, sd))
ggplot(data = Mean_Cat)+
  geom_bar(aes(x=Category,y=Mean),stat="identity")+
  theme(axis.text = element_text(size=16,colour = "black"),
        axis.title = element_text(size=18,colour = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_rect(colour = "black", fill=NA,size=1),
        panel.grid.major = element_line(size = 0.5, linetype = "solid",colour = "grey86"),
        legend.background = element_rect(fill = "transparent", colour = NA))


##### C) Consider the 4th species as a prey in the food web (one species)
# Do the same analysis as previously presented 
BernSP4_Prey<-data.frame(Simu=NA,Pred=NA)
for (i in 1:100){
  vecPred<-which(BernMat[i,which(rownames(Mb_W)=="SP4"),]!=0)
  if (length(vecPred)!=0){
    Predname<-data.frame(Simu=i,Pred=colnames(Mb_W)[vecPred])
    
  }else{Predname<-data.frame(Simu=i,Pred="No predator")}
  
  BernSP4_Prey<-rbind(BernSP4_Prey,Predname)  
  
}
BernSP4_Prey<-na.omit(BernSP4_Prey)


##### 3) Example of data analyses considering a body size gradient of a small invasive fish #####
# A) Infer the diet range for a body size gradient of invaders
# Create a data frame of the new species characteristics
Invader<-data.frame(Species=paste(rep("Invader",6),seq(1,6,1),sep = "_"),
                    Category=rep("Vertebrate",6),
                    bs=seq(20000,70000,10000),
                    Diet= rep("omnivorous",6), 
                    Habitat = rep("pel/litto",6),
                    Carnivorous = rep("0",6),
                    Log10bs=log10(seq(20000,70000,10000)) )

# Combine both data frames
DATA2<-rbind(DATA,Invader)

# Reorder the new data frame including the species gradient body size
DATA2<- DATA2[order(DATA2$bs, decreasing = TRUE),]


#Do the same steps as with one species 
bs_r2<-get_niche_attributes(species_name=DATA2$Species,
                           body_size = DATA2$Log10bs,
                           species_category = DATA2$Category)

# Inferrences of all trophic links (binary matrix)
Bmat2<-L_fn2(name=bs_r2$name,n=bs_r2$n,c=bs_r2$c,low=bs_r2$low,high=bs_r2$high,table="NO")

#Links refinement based on consumers" diet
Bmat_Diet2<-Ref_L_Diet(Bmat=Bmat2, diet=DATA2$Diet,Table="NO")

#Links refinement based on consumers" habitat
Bmat2_Hab<-Ref_L_Hab(Bmat=Bmat_Diet2,habitat =DATA2$Habitat,Table="NO")
Mb2<-Bmat2_Hab

#Weighting of refined trophic links
Mb_W2<-Weighting(Niche_attributes = bs_r2,Bmat=Mb2)

#Generate the 100 matrices from the Bernouilli trials 
BernMat2<-make_bern(n=100,Wmat=Mb_W2)


# B) Get information on prey of the invasive fish for the gradient of body sizes
INVnames<-DATA2$Species[which(str_detect(DATA2$Species, pattern="Inv", negate = FALSE))]
INV_Pred<-list()
for(j in INVnames){
  BernINV_Pred<-data.frame(Simu=NA,Prey=NA,Habitat=NA,Category=NA,PreySize=NA,Log10PreySize=NA,PredSize=NA,Log10PredSize=NA)
  PCA<-data.frame(Simu=NA,SizeINV=NA,Prey_richness=NA,InvPerc=NA,FishPerc=NA,Canibalism=NA,ZooPerc=NA,
                  LittoPerc=NA,PelPerc=NA,PelLittoPerc=NA)
  
  for (i in 1:100){
    vecPrey<-which(BernMat[i,,which(colnames(Mb_W2)==j)]!=0)
    
    if(length(vecPrey)!=0){
      Preyname<-data.frame(Simu=i,
                         Prey=rownames(Mb_W2)[vecPrey],
                         Habitat= DATA2$Habitat[match(rownames(Mb_W2)[vecPrey],DATA2$Species)],
                         Category=DATA2$Category[match(rownames(Mb_W2)[vecPrey],DATA2$Species)],
                         PreySize = DATA2$bs[match(rownames(Mb_W2)[vecPrey],DATA2$Species)], 
                         Log10PreySize=DATA2$Log10bs[match(rownames(Mb_W2)[vecPrey],DATA2$Species)],
                         PredSize=DATA2$bs[match(j,DATA2$Species)],
                         Log10PredSize=DATA2$Log10bs[match(j,DATA2$Species)])
    BernINV_Pred<-rbind(BernINV_Pred,Preyname);BernINV_Pred<-na.omit(BernINV_Pred)
    
    PCAbis<-data.frame(Simu=i,Prey_richness=nrow(Preyname),
                       SizeINV=Preyname$PredSize[1],
                       InvPerc=(100*length(which(Preyname$Category=="Invertebrate")))/nrow(Preyname),
                       FishPerc=(100*length(which(Preyname$Category=="Vertebrate"& !is.element(Preyname$Prey,Preyname$Prey[which(str_detect(Preyname$Prey, pattern="Inv", negate = FALSE))]))))/nrow(Preyname),
                       Canibalism=(100*length(which(str_detect(Preyname$Prey, pattern="Inv", negate = FALSE))))/nrow(Preyname),
                       ZooPerc=(100*length(which(Preyname$Category=="Zooplancton"))/nrow(Preyname)),
                       LittoPerc=(100*length(which(Preyname$Habitat=="litto")))/nrow(Preyname),
                       PelPerc=(100*length(which(Preyname$Habitat=="pel")))/nrow(Preyname),
                       PelLittoPerc=(100*length(which(Preyname$Habitat=="pel/litto")))/nrow(Preyname))
    }else{PCAbis<-NA}
    
    PCA<-rbind(PCA,PCAbis);PCA<-na.omit(PCA)
    
  }
  INV_Pred[[j]][["Prey"]]<-BernINV_Pred # save the list of prey for each simulation and each invasive fish body size
  INV_Pred[[j]][["PCA"]]<-PCA # save the list of prey metrics for each simulation and each invasive fish body size
}



# C) Other methods to analyze information on prey 
DAT_Prey<- rbindlist(lapply(INV_Pred, '[[', 2),use.names=TRUE) # Extract from the list

# Find the number of prey for each body size of the invasive fish 
TG_INV<-data.frame(SizeINV=unique(DAT_Prey$SizeINV),Mean=NA,SD=NA)
for (i in TG_INV$SizeINV){
  TG_INV$Mean[TG_INV$SizeINV==i]<-mean(DAT_Prey$Prey_richness[DAT_Prey$SizeINV==i])
  TG_INV$SD[TG_INV$SizeINV==i]<-sd(DAT_Prey$Prey_richness[DAT_Prey$SizeINV==i])
}

#Plot it
ggplot()+
  geom_point(data=TG_INV,aes(x=SizeINV/10000,y=Mean),size=3)+
  geom_errorbar(data=TG_INV,aes(x=SizeINV/10000,y=Mean,ymin=Mean-SD,ymax=Mean+SD))+
  xlab("Body size (cm)")+ylab("Mean number of prey")+
  theme(axis.text.x = element_text(size=16,colour = "black"),axis.text.y = element_text(size=16,colour = "black"),
        axis.title = element_text(size=18,colour = "black"),
        panel.background = element_rect(fill = "transparent", colour = NA),
        panel.border=element_rect(colour = "black", fill=NA,size=1),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',colour = "grey86"),
        legend.background = element_rect(fill = "transparent", colour = NA))


# Perform a PCA
DAT_PCA<- DAT_Prey
ColKeep<-which(colSums(DAT_PCA)!=0 & !colnames(DAT_PCA)%in%c("Simu","SizeINV"))
DAT_PCA2<-DAT_PCA[,..ColKeep]

# Standardize data an compute the PCA
resh <- decostand(DAT_PCA2, "standardize")
res.pca <- PCA(resh, graph = FALSE)

# variance explained (%) by axes
get_eigenvalue(res.pca)
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 65))

# PCA Biplot
fviz_pca_biplot(res.pca, col.var = "black", col.ind = DAT_PCA$SizeINV/10000, label="var",pointsize=3,
                addEllipses =F,ellipse.level = 0.95,repel = T)+ 
  scale_colour_gradient(low = "#FFD140", high = "#330086")+
  theme(axis.text= element_text(face="bold",color="black",size=16))


