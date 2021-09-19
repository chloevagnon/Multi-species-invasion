
##########################################################################################################################################
# Infer the niche attributes of one or more consumers 
# Necessitates to load the packages "stringr". 
# Capital or lowercase letters work for the species category. 
#
# Input:
# species_name = the name of the species to use
# body_size = log10(vector of consumer body sizes in µm)
# species_category = "vertebrate" or "vertebre" or "invertebrate" or "invertebre" or other
#    /!\ If "other" is mentionned, the species will automatically be considered as producer and not as consumer
#
# Returns : 
# a data frame with the Niche attributes of a species

get_niche_attributes<-function(species_name,body_size, species_category) {
  
  Niche<-data.frame(name=NA,n=NA,c=NA,low=NA,high=NA)
  
  for(i in 1:length(body_size)){
    
    Niche2<-data.frame(name=NA,n=NA,c=NA,low=NA,high=NA)
    
    if(!str_detect(species_category[i], regex("vertebrate", ignore_case = TRUE)) | 
       !str_detect(species_category[i], regex("invertebrate", ignore_case = TRUE)) |
       !str_detect(species_category[i], regex("vertebre", ignore_case = TRUE)) | 
       !str_detect(species_category[i], regex("invertebre", ignore_case = TRUE))|
       !str_detect(species_category[i], regex("zoop", ignore_case = TRUE))){
      # Estimate the parameters for the niche model
      Niche2$name=species_name[i]
      Niche2$n = body_size[i]						# The niche n
      Niche2$low = 0	# The lower limit of the range
      Niche2$high = 0	# The higher limit of the range
      Niche2$c = 0			# The centroid c
    }
    
    if(str_detect(species_category[i], regex("vertebrate", ignore_case = TRUE))&
       !str_detect(species_category[i], regex("inv", ignore_case = TRUE))|
       str_detect(species_category[i], regex("vertebre", ignore_case = TRUE))&
       !str_detect(species_category[i], regex("inv", ignore_case = TRUE))){
      # Unwrap the input parameters
      qrsup = Param_regvert[[2]]
      qrinf = Param_regvert[[3]]
      # Estimate the parameters for the niche model
      Niche2$name=species_name[i]
      Niche2$n = body_size[i]						# The niche n
      Niche2$low = qrinf[1] + qrinf[2]*body_size[i]	# The lower limit of the range
      Niche2$high = qrsup[1] + qrsup[2]*body_size[i]	# The higher limit of the range
      Niche2$c = Niche2$low+(Niche2$high-Niche2$low)/2			# The centroid c
    }
    
    if(str_detect(species_category[i], regex("invertebrate", ignore_case = TRUE))|
       str_detect(species_category[i], regex("invertebre", ignore_case = TRUE))|
       str_detect(species_category[i], regex("zoop", ignore_case = TRUE))){
      # Unwrap the input parameters
      qrsup = Param_reginvert[[2]]
      qrinf = Param_reginvert[[3]]
      # Estimate the parameters for the niche model
      Niche2$name=species_name[i]
      Niche2$n = body_size[i]						# The niche n
      Niche2$low = qrinf[1] + qrinf[2]*body_size[i]	# The lower limit of the range
      Niche2$high = qrsup[1] + qrsup[2]*body_size[i]	# The higher limit of the range
      Niche2$c = Niche2$low+(Niche2$high-Niche2$low)/2			# The centroid c
      
    }
    Niche<-rbind(Niche,Niche2)
  }
  
  
  return(na.omit(Niche))	
}


##########################################################################################################################################
# Transform the parameters into an interaction matrix (the metaweb) from Gravel et al. 2013
# Modified to also obtain the data frame corresponding to the binary matrix
# Input:
# name=the name of the species to use
# n = vector of size S with the parameter n for each of the S species
# c = vector of size S with the parameter c for each of the S species
# low = vector of size S with the parameter low for each of the S species
# high = vector of size S with the parameter high for each of the S species
# table="YES" or "NO" if the matrix is also wanted in data.frame counting one line per interaction
#
#
# Returns : 
# a SxS matrix with 0 indicating absence of a link and 1 indicating the presence of a link
# Predators are on columns, preys are on rows
# or a list with the SxS matrix + A table with the name of the prey, the predator and the body size of the prey and the predator

L_fn2 = function(name,n,c,low,high,table) {
  
  S = length(n)   	
  L = matrix(0,nr=S,nc=S)
  
  for(i in 1:S)
    for(j in 1:S)
      if(n[j]>low[i] && n[j]<high[i]) L[j,i] = 1
  colnames(L)<-name
  rownames(L)<-name
  
  Table<-data.frame(Prey=NA,Pred=NA,Log_Size_Prey=NA,Log_Size_Pred=NA)
  for(i in 1:S){
    if(length(which(L[,i]==1))!=0){
      Table2<-data.frame(Prey=names(which(L[,i]==1)),
                         Pred=rep(colnames(L)[i],length(which(L[,i]==1))),
                         Log_Size_Prey=n[which(L[,i]==1)],
                         Log_Size_Pred=n[i])}else{Table2<-data.frame(Prey=NA,Pred=NA,Log_Size_Prey=NA,Log_Size_Pred=NA)}
    Table<-rbind(Table,Table2)}
  if(table=="NO"){
    return(L)
  }
  if(table=="YES"){ 
    return(list(Bmat=L,Table=na.omit(Table)))
  }
}


##########################################################################################################################################
# Refine links for impossible results according to species diet trait
# Fish do not eat primary producers
# Carnivorous macroinvertebrates do not eat primary producers
# Other diet refinement must be done after applying the function as they are considered as site-specific
# Inputs:
# Bmat = binary matrix of trophic links inferred from L_fn2 function with :
#        colnames and rownames = names of species in the inventory
# diet = "invP" or "inv" for Invertebrates | "p"=piscivorous or "o"=omnivorous for fish | "prod" for primary producers 
# Table = "YES" or "NO" to obtain or note the matrix with links refined
# LinksTab = if Table = "YES" provide the table of links obtained with the function L_fn2
#
# Returns : 
# if table="NO", returns th Binary matrix after link refinement
# if table="YES", returns th Binary matrix after link refinement + the table of links after links refinement

Ref_L_Diet = function(Bmat, diet,Table,LinksTab) {
  if(Table=="NO"){
    for(i in 1:ncol(Bmat)){
      for (j in 1:nrow(Bmat)){
        if(length(which(Bmat[,i]==1))!=0){
          
          if(diet[i]=="omnivorous" & diet[j]=="prod" | 
             diet[i]=="piscivorous" & diet[j]=="prod"){
            Bmat[j,i]<-0
          }else{Bmat[j,i]<-Bmat[j,i]}
          
          if(diet[i]=="invP" & 
             diet[j]=="prod"){
            Bmat[j,i]<-0
          }else{Bmat[j,i]<-Bmat[j,i]}
          
          # if(diet[i]=="inv" & 
          #    diet[j]!="prod"){
          #   Bmat[j,i]<-0
          # }else{Bmat[j,i]<-Bmat[j,i]}
          
        }else{Bmat[,i]<-0}
      }
    }
    
    return(Bmat_ref=Bmat)	
  }
  
  if (Table == "YES") { 
    for(i in 1:ncol(Bmat)){
      for (j in 1:nrow(Bmat)){
        if(length(which(Bmat[,i]==1))!=0){
          
          if(diet[i]=="omnivorous" & diet[j]=="prod" | 
             diet[i]=="piscivorous" & diet[j]=="prod"){
            Bmat[j,i]<-0
          }else{Bmat[j,i]<-Bmat[j,i]}
          
          if(diet[i]=="invP" & 
             diet[j]=="prod"){
            Bmat[j,i]<-0
          }else{Bmat[j,i]<-Bmat[j,i]}
          
        }else{Bmat[,i]<-0}
      }
    }    
    for (i in 1:ncol(Bmat)) {
      for (j in 1:nrow(Bmat)) {
        if(Bmat[j,i]==0){
          LinksTab$Prey[LinksTab$Prey==rownames(Bmat)[j]& LinksTab$Pred==colnames(Bmat)[i]]<-NA
        }else{
          
        }
        
      }
    }
    return(list(Bmat_ref = Bmat,Table_ref=na.omit(LinksTab)))
  }
  
}


##########################################################################################################################################
# Refine links for impossible results according to the habitat trait 
# Inputs:
# Bmat = binary matrix of trophic links inferred from L_fn2 function with :
#        colnames and rownames = names of species in the inventory
# habitat = "pel" or "ben" or "pel/ben" or "pel" or "litto" or "pel/litto" 
# Table = "YES" or "NO" to obtain or note the matrix with links refined
#
# Returns : 
# if Table="NO", returns th Binary matrix after link refinement
# if Table="YES", returns th Binary matrix after link refinement + the table of links after links refinement

Ref_L_Hab = function(Bmat, habitat, Table,LinksTab) {
  if(Table=="NO"){
    for(i in 1:ncol(Bmat)){
      for (j in 1:nrow(Bmat)){
        if(length(which(Bmat[,i]==1))!=0){
          
          if(habitat[i]=="pel" & 
             habitat[j]=="ben" |
             habitat[i]=="ben" &
             habitat[j]=="pel"|
             
             habitat[i]=="pel" & 
             habitat[j]=="litto" |
             habitat[i]=="litto" &
             habitat[j]=="pel"){
            Bmat[j,i]<-0
          }else{Bmat[j,i]<-Bmat[j,i]}
          
        }else{Bmat[,i]<-0}
      }
    }
    
    return(Bmat_ref=Bmat)	
  }
  
  if (Table == "YES") { 
    for (i in 1:ncol(Bmat)) {
      for (j in 1:nrow(Bmat)) {
        if (length(which(Bmat[, i] == 1)) != 0) {
          if (habitat[i] == "pel" &
              habitat[j] == "ben" |
              habitat[i] == "ben" &
              habitat[j] == "pel"|
              habitat[i]=="pel" & 
              habitat[j]=="litto" |
              habitat[i]=="litto" &
              habitat[j]=="pel"
              
              ) {
            Bmat[j, i] <- 0
          } else{
            Bmat[j, i] <- Bmat[j, i]
          }
          
        } else{
          Bmat[, i] <- 0
        }
      }
    }
    
    for (i in 1:ncol(Bmat)) {
      for (j in 1:nrow(Bmat)) {
        if(Bmat[j,i]==0){
          LinksTab$Prey[LinksTab$Prey==rownames(Bmat)[j]& LinksTab$Pred==colnames(Bmat)[i]]<-NA
        }else{
          
        }
        
      }
    }
    return(list(Bmat_ref = Bmat,Table_ref=na.omit(LinksTab)))
  }
  
}


###########################################################################################################
# Weighting procedure: 
#
# Inputs: 
# Niche_attributes = data frame resulting from the function "get_niche_attribute" with: 
#     names = species name also used as colnames/rownames for the binary matrix
#     n = log10(species body size (µm))
#     c = optimal center of the niche (log10(µm))
#     low = lower bound of the niche range (log10(µm))
#     high =higher bound of the niche range (log10(µm))

# Bmat = initial binary interaction matrix

# Returns : the weighted interaction matrix according to the method choosen

Weighting<-function(Niche_attributes,Bmat){
  
    for(i in 1:ncol(Bmat)){
      j=0
      if(length(which(Bmat[,i]==1))!=0){
        prob_interact<-vector()
        for(j in 1:length(which(Bmat[,i]==1))){
        prey_range_vector<-seq(from=Niche_attributes$low[i],to=Niche_attributes$high[i],length.out=10000)
        prob<-scales::rescale(dnorm(seq(from=5,to=15,length.out=10000),mean=10,sd=2),to=c(0.1,0.95)) # Or 0.01 and 0.99 (Pomeranz et al. 2020)
        prob<-prob[which.min(abs(Niche_attributes$n[which(Bmat[,i]==1)[j]]-prey_range_vector))]
        prob_interact[j]<-prob
        }
        Bmat[which(Bmat[,i]==1),i]<-prob_interact
        
      }else{Bmat[which(Bmat[,i]!=1),i]<-0}
    }

  return(Bmat)
}



###########################################################################################################
# Create multiple binary matrices of interactions from weighted links
# This function necessitates to load the package Rlab
# Inputs : 
# n = number of matrix to simulate
# Wmat = initial squared matrix containing weighted links generated from the function Weighting
#
#Returns : 
# a large array of dimension n,nrow(Wmat),ncol(Wmat) where each n = a simulated interaction matrix
#

library(Rlab)

# make multiple binary matrix from a weighted matrix.
make_bern<-function(n,Wmat){
  mat_ber<-array(NA,c(n,ncol(Wmat),nrow(Wmat)))
  for(k in 1:n){
    mat_inter<-Wmat
    
    for (j in 1:ncol(Wmat)){
      for (i in 1:nrow(Wmat)){
        
        mat_inter[i,j]<- rbern(1, Wmat[i,j])
      }
    }
    mat_ber[k,,]<-  mat_inter}
  return(mat_ber)
}


