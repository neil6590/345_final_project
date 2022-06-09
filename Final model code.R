### CELL LIST FUNCTIONS ###

free_cells_fun <- function(dfname) {
  dfname <- as.data.frame(which(dfname == 0, arr.ind = TRUE)) 
  dfname <- dfname[order(dfname[,1],dfname[,2], decreasing=FALSE),]
  return(dfname)
}

occu_cells_fun <- function(dfname) {
  dfname <- as.data.frame(which(dfname != 0, arr.ind = TRUE)) 
  dfname <- dfname[order(dfname[,1],dfname[,2], decreasing=FALSE),]
  return(dfname)
}


### DEATHSTEP FUNCTION FOR INDIVIDUAL DEATH ###


Polyp_DeathStep <- function(community){

  DeathRow <- sample(nrow(community),1)
  DeathCol <- sample(ncol(community),1) # so sampling a random point in the community, could be filled or could be empty
  
  if(community[DeathRow, DeathCol] !=0){     # action only takes place if it is filled. If it is empty it returns the matrix as is
    community[DeathRow, DeathCol] <- 0

  }
  
  return(community)
}


### DEATHSTEP FUNCTION FOR META-INDIVIDUAL DEATH ###

ClearClump <- function(x,y,Communitymatrix,speciesID) { #this function will be used inside the death step
  
  
  if(y >= 1 && y <= ncol(Communitymatrix) && x >= 1 && x <= nrow(Communitymatrix)) { # so if the point is in the matrix then move forward
    
    if(Communitymatrix[x,y] == speciesID){  # double checking that it is the same species ID
    
      Communitymatrix[x,y] <- 0  #kill point in community

      Communitymatrix <- ClearClump((x-1),y, Communitymatrix,speciesID)  # self-referring function will carry the same actions out in each compass direction 
      Communitymatrix <- ClearClump((x+1),y, Communitymatrix, speciesID) # and for each neighbouring cell it checks whether or not they meet the conditions
      Communitymatrix <- ClearClump(x,(y-1), Communitymatrix, speciesID)
      Communitymatrix <- ClearClump(x,(y+1), Communitymatrix, speciesID)
      
    }
  }
  return(Communitymatrix)
  
}



Colony_DeathStep <- function(community,polyp_prob=0.95){  # so polyp probability is probability of death carried out being an individual death
  DeathRow <- sample(nrow(community),1)
  DeathCol <- sample(ncol(community),1)
  
  
  if(community[DeathRow, DeathCol]!=0)
  {
    x <- runif(1)               # running here the probability argument to decide which death step it is going to undertake
    
    if(x < polyp_prob){
      community[DeathRow, DeathCol] <- 0
      
    }else{
      speciesid <- community[DeathRow, DeathCol]        # assigning species ID here as it needs to be passed into argument of clear clump function
      community <- ClearClump(DeathRow,DeathCol, community, speciesid)
      
    }
    
    print(c(DeathRow, DeathCol))
    
  }
  return(community)
}


### BIRTHSTEP FUNCTION FOR EITHER MODEL + PREREQUISITES ###


BirthStep <- function(community, speciation =0.005, immigration = 0.01)
{
  free_cell_list <- free_cells_fun(community)         #Always picks a free cell ensuring 1 birth in each timestep
  birth_success <- 0                                  # this will be the condition to break the while loop below
  while(birth_success == 0){
    
    birthp <- sample(1:nrow(free_cell_list),1)        # so for each running of the while loop, will choose a new birth cell to act on
    BirthRow <- free_cell_list[birthp,1]
    BirthCol <- free_cell_list[birthp,2]
    
    #while loop here for birth success, and until it is achieved, run the steps below sequentially
    #at the bottom of all, have attached an if statement saying birth success <- 1
    
    
    # while loop also makes biological sense as the more surrounded by 0s, the lower the chance it is populated
    
    h <- runif(1)       # this value will decide which if statement the function goes on to follow/how it will fill up the empty cell
    
    if(h <= speciation){ # relates to prob of speciation
      
      community[BirthRow,BirthCol] <-  runif(1)
      birth_success <- 1
      
      if(birth_success ==1){
        break
      }
      
    } else if(h <= immigration){
      
      occu_cell_list <- occu_cells_fun(community)
      motherp <- sample(1:nrow(occu_cell_list),1) #Choosing random mother cell from local species pool to replace it with
      MotherRow <- occu_cell_list[motherp,1]
      MotherCol <- occu_cell_list[motherp,2]
      community[BirthRow, BirthCol] <- community[MotherRow, MotherCol] 
      
      birth_success <- 1                                                                                                 
      if(birth_success ==1){
        break
      }
      
    } else { #colonisation of empty space by a neighbouring species 
      
      w <- runif(1)
      
      if(w <=0.25){
        
        if(BirthRow+1 >= 1 && (BirthRow+1) <= nrow(community) && (BirthCol >= 1) && BirthCol <= ncol(community) && (community[BirthRow+1, BirthCol] != 0)){ 
          community[BirthRow, BirthCol] <- community[BirthRow+1, BirthCol]   }  # above are the conditions ensuring the neighbouring cell is within the 
                                                                                # bounds of the matrix (imp for birth cells on the edges of the matrix)
      }                                                                         # and also checks that the neighbouring cell is not equal to 0
      
      else if(w <= 0.5){  # each neighbour has an equal chance of being chosen
        
        if(BirthRow-1 >= 1 && BirthRow-1 <= nrow(community) && BirthCol >= 1 && BirthCol <= ncol(community) && community[BirthRow-1, BirthCol] !=0){
          community[BirthRow, BirthCol] <- community[BirthRow-1, BirthCol]     }
      }
      
      else if(w <= 0.75){
        
        if(BirthRow >= 1 && BirthRow <= nrow(community) && BirthCol+1 >= 1 && BirthCol+1 <= ncol(community) && community[BirthRow, BirthCol+1] != 0){
          community[BirthRow,BirthCol] <- community[BirthRow,BirthCol+1]      }
      }
      
      else{
        
        if(BirthRow >= 1 && BirthRow <= nrow(community) && BirthCol-1 >= 1 && BirthCol-1 <= ncol(community) && community[BirthRow, BirthCol-1] !=0){
          community[BirthRow, BirthCol] <- community[BirthRow, BirthCol-1]      }
        
      }
      
      birth_success <- 1
      if(birth_success ==1){
        break
        
        
      }
    }
    
  }
  return(community)
  
}



### RUNSIM FUNCTION ###

 # this death step has only individual death
RunSimulation_polyp <- function(community,speciation=0.001, immigration = 0.01,NumSteps=100,000)
{

  for (i in 1:(NumSteps)){
    
    community <- Polyp_DeathStep(community)
    
    
    community <- BirthStep(community, speciation, immigration)
    
    
  }
  return(community) 
}


  #so this death step has chance of both individual and meta-individual death
RunSimulation_colony <- function(community,speciation=0.001, immigration = 0.01, polyp_prob = 0.9, NumSteps=100)
{
  for (i in 1:(NumSteps)){ 
    community <- Colony_DeathStep(community, polyp_prob)
    community <- BirthStep(community, speciation, immigration)
  }
  return(community) 
}






### CREATING AND INITIALISING LANDSCAPES ###

SetupCommunityMatrix <- function(width,height,NumberSpecies=2)
{
  CommunityMatrix <- matrix(sample(x = 0:NumberSpecies, size = (width*height), replace = TRUE),nrow=height,ncol=width) #So matrix with sampling between 1-no of species)
  return(CommunityMatrix)
}

CommunityMatrix1 <- SetupCommunityMatrix(250,250,30)   # Three 250X250 matrices created
 
CommunityMatrix2 <- SetupCommunityMatrix(250,250,30)

CommunityMatrix3 <- SetupCommunityMatrix(250,250,30)

## BURN-IN FOR LANDSCAPES ##

initial_comm_1 <- RunSimulation_polyp(CommunityMatrix1, NumSteps=100000)

initial_comm_2 <- RunSimulation_polyp(CommunityMatrix2, NumSteps=100000)

initial_comm_3 <- RunSimulation_polyp(CommunityMatrix3, NumSteps=100000)




### RUNNING THE MODELS ON THE LANDSCAPES ###



comm_im1_pp95_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.1)
comm_im1_pp95_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.1)
comm_im1_pp95_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.1)


comm_im1_pp99_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.1)
comm_im1_pp99_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.1)
comm_im1_pp99_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.1)


comm_im1_pp1_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 1, immigration = 0.1)
comm_im1_pp1_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 1, immigration = 0.1)
comm_im1_pp1_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 1, immigration = 0.1)





comm_im01_pp95_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.01)
comm_im01_pp95_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.01)
comm_im01_pp95_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 0.95, immigration = 0.01)


comm_im01_pp99_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.01)
comm_im01_pp99_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.01)
comm_im01_pp99_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 0.99, immigration = 0.01)


comm_im01_pp1_1 <- RunSimulation_colony(initial_comm_1, NumSteps = 200000, polyp_prob = 1, immigration = 0.01)
comm_im01_pp1_2 <- RunSimulation_colony(initial_comm_2, NumSteps = 200000, polyp_prob = 1, immigration = 0.01)
comm_im01_pp1_3 <- RunSimulation_colony(initial_comm_3, NumSteps = 200000, polyp_prob = 1, immigration = 0.01)



####### CALCULATING INDICES AND CREATING VISUAL PLOTS #######

## Species richness ##

species_richness <- function(community){
  mylist <- c()
  for(r in 1:nrow(community)){
    for(c in 1:ncol(community)){
      mylist <- rbind(mylist, community[r,c])
    }
  }
  return(length(unique(mylist))-1) # -1 here to not include 0 values in the empty cells
  
}

##CALCULATING SPECIES RICHNESS FOR INDIVIDUALS##

comm_im1_pp95_1.SR <- species_richness(comm_im1_pp95_1)
comm_im1_pp95_2.SR <- species_richness(comm_im1_pp95_2)
comm_im1_pp95_3.SR <- species_richness(comm_im1_pp95_3)
comm_im1_pp95_avg.err.SR <- std.error(c(comm_im1_pp95_1.SR,comm_im1_pp95_2.SR,comm_im1_pp95_3.SR)) # calculating std.error

comm_im1_pp95_avg.SR <- (comm_im1_pp95_1.SR + comm_im1_pp95_2.SR + comm_im1_pp95_3.SR)/3  # Calculating value for each replicate individually, then dividing it by 3


comm_im1_pp99_1.SR <- species_richness(comm_im1_pp99_1)
comm_im1_pp99_2.SR <- species_richness(comm_im1_pp99_2)
comm_im1_pp99_3.SR <- species_richness(comm_im1_pp99_2)
comm_im1_pp99_avg.err.SR <- std.error(c(comm_im1_pp99_1.SR,comm_im1_pp99_2.SR,comm_im1_pp99_3.SR))


comm_im1_pp99_avg.SR <- (comm_im1_pp99_1.SR+comm_im1_pp99_2.SR+comm_im1_pp99_3.SR)/3


comm_im1_pp1_1.SR <- species_richness(comm_im1_pp1_1)
comm_im1_pp1_2.SR <- species_richness(comm_im1_pp1_2)
comm_im1_pp1_3.SR <- species_richness(comm_im1_pp1_3)
comm_im1_pp1_avg.err.SR <- std.error(c(comm_im1_pp1_1.SR,comm_im1_pp1_2.SR,comm_im1_pp1_3.SR))

comm_im1_pp1_avg.SR <- (comm_im1_pp100_1.SR+comm_im1_pp100_2.SR+comm_im1_pp100_3.SR)/3



comm_im01_pp95_1.SR <- species_richness(comm_im01_pp95_1)
comm_im01_pp95_2.SR <- species_richness(comm_im01_pp95_2)
comm_im01_pp95_3.SR <- species_richness(comm_im01_pp95_3)
comm_im01_pp95_avg.err.SR <- std.error(c(comm_im01_pp95_1.SR,comm_im01_pp95_2.SR,comm_im01_pp95_3.SR))

comm_im01_pp95_avg.SR <- (comm_im01_pp95_1.SR + comm_im01_pp95_2.SR + comm_im01_pp95_3.SR)/3


comm_im01_pp99_1.SR <- species_richness(comm_im01_pp99_1)
comm_im01_pp99_2.SR <- species_richness(comm_im01_pp99_2)
comm_im01_pp99_3.SR <- species_richness(comm_im01_pp99_2)
comm_im01_pp99_avg.err.SR <- std.error(c(comm_im01_pp99_1.SR,comm_im01_pp99_2.SR,comm_im01_pp99_3.SR))

comm_im01_pp99_avg.SR <- (comm_im01_pp99_1.SR+comm_im01_pp99_2.SR+comm_im01_pp99_3.SR)/3


comm_im01_pp1_1.SR <- species_richness(comm_im01_pp1_1)
comm_im01_pp1_2.SR <- species_richness(comm_im01_pp1_2)
comm_im01_pp1_3.SR <- species_richness(comm_im01_pp1_3)
comm_im01_pp1_avg.err.SR <- std.error(c(comm_im01_pp1_1.SR,comm_im01_pp1_2.SR,comm_im01_pp1_3.SR))

comm_im01_pp1_avg.SR <- (comm_im01_pp1_1.SR+comm_im01_pp1_2.SR+comm_im01_pp1_3.SR)/3

SR.no<- c(comm_im1_pp95_avg.SR,comm_im1_pp99_avg.SR,comm_im1_pp1_avg.SR,comm_im01_pp95_avg.SR,comm_im01_pp99_avg.SR,comm_im01_pp1_avg.SR)
SR.names <- c("comm_im1_pp95_avg.SR", "comm_im1_pp99_avg.SR", "comm_im1_pp1_avgSR", "comm_im01_pp95_avg.SR", "comm_im01_pp99_avg.SR", "comm_im01_pp1_avg.SR" )
SR.error <- c(comm_im1_pp95_avg.err.SR, comm_im1_pp99_avg.err.SR, comm_im1_pp1_avg.err.SR, comm_im01_pp95_avg.err.SR, comm_im01_pp99_avg.err.SR, comm_im01_pp1_avg.err.SR)
SR <- rbind(SR.names, SR.no, SR.error)  #compiling data into a table






## shannon diversity + Effective no of species ## 

eff_no_species <- function(df){
  df.N <- sum(df$Freq)
  df.prop <- df$Freq/df.N
  df.mult <- df.prop*log(df.prop)
  df.shannon <- 0 - sum(df.mult)
  df.effective <- exp(df.shannon)
  return(df.effective)
  
}



## Calculating values for each sim ## 

comm_im1_pp95_1.EF <- eff_no_species(comm_im1_pp95_1.spp)
comm_im1_pp95_2.EF <- eff_no_species(comm_im1_pp95_2.spp)
comm_im1_pp95_3.EF <- eff_no_species(comm_im1_pp95_3.spp)
comm_im1_pp95_avg.err.EF <- std.error(c(comm_im1_pp95_1.EF,comm_im1_pp95_2.EF,comm_im1_pp95_3.EF)) # carrying out same actions as above, just for EF values now

comm_im1_pp95_avg.EF <- (comm_im1_pp95_1.EF+comm_im1_pp95_2.EF+comm_im1_pp95_3.EF)/3

comm_im1_pp99_1.EF <- eff_no_species(comm_im1_pp99_1.spp)
comm_im1_pp99_2.EF <- eff_no_species(comm_im1_pp99_2.spp)
comm_im1_pp99_3.EF <- eff_no_species(comm_im1_pp99_3.spp)
comm_im1_pp99_avg.err.EF <- std.error(c(comm_im1_pp99_1.EF,comm_im1_pp99_2.EF,comm_im1_pp99_3.EF))

comm_im1_pp99_avg.EF <- (comm_im1_pp99_1.EF+comm_im1_pp99_2.EF+comm_im1_pp99_3.EF)/3

comm_im1_pp1_1.EF <- eff_no_species(comm_im1_pp1_1.spp)
comm_im1_pp1_2.EF <- eff_no_species(comm_im1_pp1_2.spp)
comm_im1_pp1_3.EF <- eff_no_species(comm_im1_pp1_3.spp)
comm_im1_pp1_avg.err.EF <- std.error(c(comm_im1_pp1_1.EF,comm_im1_pp1_2.EF,comm_im1_pp1_3.EF))

comm_im1_pp1_avg.EF <- (comm_im1_pp1_1.EF+comm_im1_pp1_2.EF+comm_im1_pp1_3.EF)/3



comm_im01_pp95_1.EF <- eff_no_species(comm_im01_pp95_1.spp)
comm_im01_pp95_2.EF <- eff_no_species(comm_im01_pp95_2.spp)
comm_im01_pp95_3.EF <- eff_no_species(comm_im01_pp95_3.spp)
comm_im01_pp95_avg.err.EF <- std.error(c(comm_im01_pp95_1.EF,comm_im01_pp95_2.EF,comm_im01_pp95_3.EF))

comm_im01_pp95_avg.EF <- (comm_im01_pp95_1.EF+comm_im01_pp95_2.EF+comm_im01_pp95_3.EF)/3

comm_im01_pp99_1.EF <- eff_no_species(comm_im01_pp99_1.spp)
comm_im01_pp99_2.EF <- eff_no_species(comm_im01_pp99_2.spp)
comm_im01_pp99_3.EF <- eff_no_species(comm_im01_pp99_3.spp)
comm_im01_pp99_avg.err.EF <- std.error(c(comm_im01_pp99_1.EF,comm_im01_pp99_2.EF,comm_im01_pp99_3.EF))

comm_im01_pp99_avg.EF <- (comm_im01_pp99_1.EF+comm_im01_pp99_2.EF+comm_im01_pp99_3.EF)/3

comm_im01_pp1_1.EF <- eff_no_species(comm_im01_pp1_1.spp)
comm_im01_pp1_2.EF <- eff_no_species(comm_im01_pp1_2.spp)
comm_im01_pp1_3.EF <- eff_no_species(comm_im01_pp1_3.spp)
comm_im01_pp1_avg.err.EF <- std.error(c(comm_im01_pp1_1.EF,comm_im01_pp1_2.EF,comm_im01_pp1_3.EF))

comm_im01_pp1_avg.EF <- (comm_im01_pp1_1.EF+comm_im01_pp1_2.EF+comm_im01_pp1_3.EF)/3



EF.no<- c(comm_im1_pp95_avg.EF,comm_im1_pp99_avg.EF,comm_im1_pp1_avg.EF,comm_im01_pp95_avg.EF,comm_im01_pp99_avg.EF,comm_im01_pp1_avg.EF)
EF.names <- c("comm_im1_pp95_avg.EF", "comm_im1_pp99_avg.EF", "comm_im1_pp1_avgEF", "comm_im01_pp95_avg.EF", "comm_im01_pp99_avg.EF", "comm_im01_pp1_avg.EF" )
EF.error <- c(comm_im1_pp95_avg.err.EF, comm_im1_pp99_avg.err.EF, comm_im1_pp1_avg.err.EF, comm_im01_pp95_avg.err.EF, comm_im01_pp99_avg.err.EF, comm_im01_pp1_avg.err.EF)
EF <- rbind(EF.names, EF.no, EF.error)
View(EF)
###


## Berger-Parker Index of Dominance ## 

BP_IoD <- function(df){
  df.BP <- max(df$Freq)/sum(df$Freq)
  return(df.BP)
}


}

## For Real Sims   # carrying out same actions as above, just for BP values now


comm_im1_pp95_1.BP <- BP_IoD(comm_im1_pp95_1.spp)
comm_im1_pp95_2.BP <- BP_IoD(comm_im1_pp95_2.spp)
comm_im1_pp95_3.BP <- BP_IoD(comm_im1_pp95_3.spp)
comm_im1_pp95_avg.err.BP <- std.error(c(comm_im1_pp95_1.BP,comm_im1_pp95_2.BP,comm_im1_pp95_3.BP))

comm_im1_pp95_avg.BP <- (comm_im1_pp95_1.BP+comm_im1_pp95_2.BP+comm_im1_pp95_3.BP)/3

comm_im1_pp99_1.BP <- BP_IoD(comm_im1_pp99_1.spp)
comm_im1_pp99_2.BP <- BP_IoD(comm_im1_pp99_2.spp)
comm_im1_pp99_3.BP <- BP_IoD(comm_im1_pp99_3.spp)
comm_im1_pp99_avg.err.BP <- std.error(c(comm_im1_pp99_1.BP,comm_im1_pp99_2.BP,comm_im1_pp99_3.BP))

comm_im1_pp99_avg.BP <- (comm_im1_pp99_1.BP+comm_im1_pp99_2.BP+comm_im1_pp99_3.BP)/3

comm_im1_pp1_1.BP <- BP_IoD(comm_im1_pp1_1.spp)
comm_im1_pp1_2.BP <- BP_IoD(comm_im1_pp1_2.spp)
comm_im1_pp1_3.BP <- BP_IoD(comm_im1_pp1_3.spp)
comm_im1_pp1_avg.err.BP <- std.error(c(comm_im1_pp1_1.BP,comm_im1_pp1_2.BP,comm_im1_pp1_3.BP))

comm_im1_pp1_avg.BP <- (comm_im1_pp1_1.BP+comm_im1_pp1_2.BP+comm_im1_pp1_3.BP)/3



comm_im01_pp95_1.BP <- BP_IoD(comm_im01_pp95_1.spp)
comm_im01_pp95_2.BP <- BP_IoD(comm_im01_pp95_2.spp)
comm_im01_pp95_3.BP <- BP_IoD(comm_im01_pp95_3.spp)
comm_im01_pp95_avg.err.BP <- std.error(c(comm_im01_pp95_1.BP,comm_im01_pp95_2.BP,comm_im01_pp95_3.BP))

comm_im01_pp95_avg.BP <- (comm_im01_pp95_1.BP+comm_im01_pp95_2.BP+comm_im01_pp95_3.BP)/3

comm_im01_pp99_1.BP <- BP_IoD(comm_im01_pp99_1.spp)
comm_im01_pp99_2.BP <- BP_IoD(comm_im01_pp99_2.spp)
comm_im01_pp99_3.BP <- BP_IoD(comm_im01_pp99_3.spp)
comm_im01_pp99_avg.err.BP <- std.error(c(comm_im01_pp99_1.BP,comm_im01_pp99_2.BP,comm_im01_pp99_3.BP))

comm_im01_pp99_avg.BP <- (comm_im01_pp99_1.BP+comm_im01_pp99_2.BP+comm_im01_pp99_3.BP)/3

comm_im01_pp1_1.BP <- BP_IoD(comm_im01_pp1_1.spp)
comm_im01_pp1_2.BP <- BP_IoD(comm_im01_pp1_2.spp)
comm_im01_pp1_3.BP <- BP_IoD(comm_im01_pp1_3.spp)
comm_im01_pp1_avg.err.BP <- std.error(c(comm_im01_pp1_1.BP,comm_im01_pp1_2.BP,comm_im01_pp1_3.BP))

comm_im01_pp1_avg.BP <- (comm_im01_pp1_1.BP+comm_im01_pp1_2.BP+comm_im01_pp1_3.BP)/3



BP.no<- c(comm_im1_pp95_avg.BP,comm_im1_pp99_avg.BP,comm_im1_pp1_avg.BP,comm_im01_pp95_avg.BP,comm_im01_pp99_avg.BP,comm_im01_pp1_avg.BP)
BP.names <- c("comm_im1_pp95_avg.BP", "comm_im1_pp99_avg.BP", "comm_im1_pp1_avgBP", "comm_im01_pp95_avg.BP", "comm_im01_pp99_avg.BP", "comm_im01_pp1_avg.BP" )
BP.error <- c(comm_im1_pp95_avg.err.BP, comm_im1_pp99_avg.err.BP, comm_im1_pp1_avg.err.BP, comm_im01_pp95_avg.err.BP, comm_im01_pp99_avg.err.BP, comm_im01_pp1_avg.err.BP)
BP <- rbind(BP.names, BP.no, BP.error)
View(BP)





### REQUIRED TRANSFORMATIONS TO BE ABLE TO PLOT RESULTS AS SPECIES ABUNDANCE DISTRIBUTIONS (SADs) ###


comm_im1_pp95_1.freq <- as.data.frame(table(comm_im1_pp95_1))  #converting table of frequencies for each species ID to a data frame
comm_im1_pp95_1.sort <- comm_im1_pp95_1.freq[order(-comm_im1_pp95_1.freq$Freq),]  #re-ordering the dataframe in terms of decreasing rank abundance
comm_im1_pp95_1.spp <- subset(comm_im1_pp95_1.sort, comm_im1_pp95_1.sort$comm_im1_pp95_1!=0)   #removing frequencies of 0s in the dataframe
comm_im1_pp95_1.spp$rank <- c(1:nrow(comm_im1_pp95_1.spp))      # assigning rank values according to re-ordering that happened for "comm_im1_pp95_1.sort"


comm_im1_pp95_2.freq <- as.data.frame(table(comm_im1_pp95_2))    #doing this for each replicate of every combination
comm_im1_pp95_2.sort <- comm_im1_pp95_2.freq[order(-comm_im1_pp95_2.freq$Freq),]
comm_im1_pp95_2.spp <- subset(comm_im1_pp95_2.sort, comm_im1_pp95_2.sort$comm_im1_pp95_2!=0)
comm_im1_pp95_2.spp$rank <- c(1:nrow(comm_im1_pp95_2.spp))


comm_im1_pp95_3.freq <- as.data.frame(table(comm_im1_pp95_3))
comm_im1_pp95_3.sort <- comm_im1_pp95_3.freq[order(-comm_im1_pp95_3.freq$Freq),]
comm_im1_pp95_3.spp <- subset(comm_im1_pp95_3.sort, comm_im1_pp95_3.sort$comm_im1_pp95_3!=0)
comm_im1_pp95_3.spp$rank <- c(1:nrow(comm_im1_pp95_3.spp))



comm_im1_pp99_1.freq <- as.data.frame(table(comm_im1_pp99_1))
comm_im1_pp99_1.sort <- comm_im1_pp99_1.freq[order(-comm_im1_pp99_1.freq$Freq),]
comm_im1_pp99_1.spp <- subset(comm_im1_pp99_1.sort, comm_im1_pp99_1.sort$comm_im1_pp99_1!=0)
comm_im1_pp99_1.spp$rank <- c(1:nrow(comm_im1_pp99_1.spp))


comm_im1_pp99_2.freq <- as.data.frame(table(comm_im1_pp99_2))
comm_im1_pp99_2.sort <- comm_im1_pp99_2.freq[order(-comm_im1_pp99_2.freq$Freq),]
comm_im1_pp99_2.spp <- subset(comm_im1_pp99_2.sort, comm_im1_pp99_2.sort$comm_im1_pp99_2!=0)
comm_im1_pp99_2.spp$rank <- c(1:nrow(comm_im1_pp99_2.spp))


comm_im1_pp99_3.freq <- as.data.frame(table(comm_im1_pp99_3))
comm_im1_pp99_3.sort <- comm_im1_pp99_3.freq[order(-comm_im1_pp99_3.freq$Freq),]
comm_im1_pp99_3.spp <- subset(comm_im1_pp99_3.sort, comm_im1_pp99_3.sort$comm_im1_pp99_3!=0)
comm_im1_pp99_3.spp$rank <- c(1:nrow(comm_im1_pp99_3.spp))



comm_im1_pp1_1.freq <- as.data.frame(table(comm_im1_pp1_1))
comm_im1_pp1_1.sort <- comm_im1_pp1_1.freq[order(-comm_im1_pp1_1.freq$Freq),]
comm_im1_pp1_1.spp <- subset(comm_im1_pp1_1.sort, comm_im1_pp1_1.sort$comm_im1_pp1_1!=0)
comm_im1_pp1_1.spp$rank <- c(1:nrow(comm_im1_pp1_1.spp))


comm_im1_pp1_2.freq <- as.data.frame(table(comm_im1_pp1_2))
comm_im1_pp1_2.sort <- comm_im1_pp1_2.freq[order(-comm_im1_pp1_2.freq$Freq),]
comm_im1_pp1_2.spp <- subset(comm_im1_pp1_2.sort, comm_im1_pp1_2.sort$comm_im1_pp1_2!=0)
comm_im1_pp1_2.spp$rank <- c(1:nrow(comm_im1_pp1_2.spp))


comm_im1_pp1_3.freq <- as.data.frame(table(comm_im1_pp1_3))
comm_im1_pp1_3.sort <- comm_im1_pp1_3.freq[order(-comm_im1_pp1_3.freq$Freq),]
comm_im1_pp1_3.spp <- subset(comm_im1_pp1_3.sort, comm_im1_pp1_3.sort$comm_im1_pp1_3!=0)
comm_im1_pp1_3.spp$rank <- c(1:nrow(comm_im1_pp1_3.spp))





comm_im01_pp95_1.freq <- as.data.frame(table(comm_im01_pp95_1))
comm_im01_pp95_1.sort <- comm_im01_pp95_1.freq[order(-comm_im01_pp95_1.freq$Freq),]
comm_im01_pp95_1.spp <- subset(comm_im01_pp95_1.sort, comm_im01_pp95_1.sort$comm_im01_pp95_1!=0)
comm_im01_pp95_1.spp$rank <- c(1:nrow(comm_im01_pp95_1.spp))


comm_im01_pp95_2.freq <- as.data.frame(table(comm_im01_pp95_2))
comm_im01_pp95_2.sort <- comm_im01_pp95_2.freq[order(-comm_im01_pp95_2.freq$Freq),]
comm_im01_pp95_2.spp <- subset(comm_im01_pp95_2.sort, comm_im01_pp95_2.sort$comm_im01_pp95_2!=0)
comm_im01_pp95_2.spp$rank <- c(1:nrow(comm_im01_pp95_2.spp))


comm_im01_pp95_3.freq <- as.data.frame(table(comm_im01_pp95_3))
comm_im01_pp95_3.sort <- comm_im01_pp95_3.freq[order(-comm_im01_pp95_3.freq$Freq),]
comm_im01_pp95_3.spp <- subset(comm_im01_pp95_3.sort, comm_im01_pp95_3.sort$comm_im01_pp95_3!=0)
comm_im01_pp95_3.spp$rank <- c(1:nrow(comm_im01_pp95_3.spp))



comm_im01_pp99_1.freq <- as.data.frame(table(comm_im01_pp99_1))
comm_im01_pp99_1.sort <- comm_im01_pp99_1.freq[order(-comm_im01_pp99_1.freq$Freq),]
comm_im01_pp99_1.spp <- subset(comm_im01_pp99_1.sort, comm_im01_pp99_1.sort$comm_im01_pp99_1!=0)
comm_im01_pp99_1.spp$rank <- c(1:nrow(comm_im01_pp99_1.spp))


comm_im01_pp99_2.freq <- as.data.frame(table(comm_im01_pp99_2))
comm_im01_pp99_2.sort <- comm_im01_pp99_2.freq[order(-comm_im01_pp99_2.freq$Freq),]
comm_im01_pp99_2.spp <- subset(comm_im01_pp99_2.sort, comm_im01_pp99_2.sort$comm_im01_pp99_2!=0)
comm_im01_pp99_2.spp$rank <- c(1:nrow(comm_im01_pp99_2.spp))


comm_im01_pp99_3.freq <- as.data.frame(table(comm_im01_pp99_3))
comm_im01_pp99_3.sort <- comm_im01_pp99_3.freq[order(-comm_im01_pp99_3.freq$Freq),]
comm_im01_pp99_3.spp <- subset(comm_im01_pp99_3.sort, comm_im01_pp99_3.sort$comm_im01_pp99_3!=0)
comm_im01_pp99_3.spp$rank <- c(1:nrow(comm_im01_pp99_3.spp))



comm_im01_pp1_1.freq <- as.data.frame(table(comm_im01_pp1_1))
comm_im01_pp1_1.sort <- comm_im01_pp1_1.freq[order(-comm_im01_pp1_1.freq$Freq),]
comm_im01_pp1_1.spp <- subset(comm_im01_pp1_1.sort, comm_im01_pp1_1.sort$comm_im01_pp1_1!=0)
comm_im01_pp1_1.spp$rank <- c(1:nrow(comm_im01_pp1_1.spp))


comm_im01_pp1_2.freq <- as.data.frame(table(comm_im01_pp1_2))
comm_im01_pp1_2.sort <- comm_im01_pp1_2.freq[order(-comm_im01_pp1_2.freq$Freq),]
comm_im01_pp1_2.spp <- subset(comm_im01_pp1_2.sort, comm_im01_pp1_2.sort$comm_im01_pp1_2!=0)
comm_im01_pp1_2.spp$rank <- c(1:nrow(comm_im01_pp1_2.spp))


comm_im01_pp1_3.freq <- as.data.frame(table(comm_im01_pp1_3))
comm_im01_pp1_3.sort <- comm_im01_pp1_3.freq[order(-comm_im01_pp1_3.freq$Freq),]
comm_im01_pp1_3.spp <- subset(comm_im01_pp1_3.sort, comm_im01_pp1_3.sort$comm_im01_pp1_3!=0)
comm_im01_pp1_3.spp$rank <- c(1:nrow(comm_im01_pp1_3.spp))



comm_im01_pp1_1.spp.graph <- comm_im01_pp1_1.spp[1:131,1:3]  #selecting only for relevant columns from the dataframe
comm_im01_pp1_2.spp.graph <- comm_im01_pp1_2.spp[1:131,1:3]
comm_im01_pp1_3.spp.graph <- comm_im01_pp1_3.spp[1:131,1:3]

comm_im01_pp1_avg.spp.graph <- (comm_im01_pp1_1.spp.graph+comm_im01_pp1_2.spp.graph+comm_im01_pp1_3.spp.graph)/3   # averaging rank abundances across the replicates


comm_im01_pp95_1.spp.graph <- comm_im01_pp95_1.spp[1:121,1:3]
comm_im01_pp95_2.spp.graph <- comm_im01_pp95_2.spp[1:121,1:3]
comm_im01_pp95_3.spp.graph <- comm_im01_pp95_3.spp[1:121,1:3]

comm_im01_pp95_avg.spp.graph <- (comm_im01_pp95_1.spp.graph+comm_im01_pp95_2.spp.graph+comm_im01_pp95_3.spp.graph)/3



comm_im01_pp99_1.spp.graph <- comm_im01_pp99_1.spp[1:137,1:3]
comm_im01_pp99_2.spp.graph <- comm_im01_pp99_2.spp[1:137,1:3]
comm_im01_pp99_3.spp.graph <- comm_im01_pp99_3.spp[1:137,1:3]

comm_im01_pp99_avg.spp.graph <- (comm_im01_pp99_1.spp.graph+comm_im01_pp99_2.spp.graph+comm_im01_pp99_3.spp.graph)/3




comm_im1_pp1_1.spp.graph <- comm_im1_pp1_1.spp[1:139,1:3]
comm_im1_pp1_2.spp.graph <- comm_im1_pp1_2.spp[1:139,1:3]
comm_im1_pp1_3.spp.graph <- comm_im1_pp1_3.spp[1:139,1:3]

comm_im1_pp1_avg.spp.graph <- (comm_im1_pp1_1.spp.graph+comm_im1_pp1_2.spp.graph+comm_im1_pp1_3.spp.graph)/3



comm_im1_pp95_1.spp.graph <- comm_im1_pp95_1.spp[1:133,1:3]
comm_im1_pp95_2.spp.graph <- comm_im1_pp95_2.spp[1:133,1:3]
comm_im1_pp95_3.spp.graph <- comm_im1_pp95_3.spp[1:133,1:3]

comm_im1_pp95_avg.spp.graph <- (comm_im1_pp95_1.spp.graph+comm_im1_pp95_2.spp.graph+comm_im1_pp95_3.spp.graph)/3



comm_im1_pp99_1.spp.graph <- comm_im1_pp99_1.spp[1:141,1:3]
comm_im1_pp99_2.spp.graph <- comm_im1_pp99_2.spp[1:141,1:3]
comm_im1_pp99_3.spp.graph <- comm_im1_pp99_3.spp[1:141,1:3]

comm_im1_pp99_avg.spp.graph <- (comm_im1_pp99_1.spp.graph+comm_im1_pp99_2.spp.graph+comm_im1_pp99_3.spp.graph)/3



library(ggplot2)  # loading ggplot package


  
comm_im01_avg.spp.graph <- cbind(comm_im01_pp95_avg.spp.graph,comm_im01_pp99_avg.spp.graph[1:121,2],comm_im01_pp1_avg.spp.graph[1:121,2])
#combining all the 3 averaged datasets into one 


comm_im_01_graph <- ggplot(comm_im01_avg.spp.graph, aes(x = reorder(rank, -Freq), group=1)) +  #ensuring it is in decreasing order of frequency on y axis
  
  geom_line(aes(y = Freq), size=1.5, colour ="black") +
  geom_line(aes(y = comm_im01_pp99_avg.spp.graph[1:121, 2]), size=1.5, colour ="blue") +
  geom_line(aes(y = comm_im01_pp1_avg.spp.graph[1:121, 2]), size=1.5, colour ="green") +   #specific colours chosen not to disadvantage those who may have certain types of colourblindness
  
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 15))+ #adjusting axis text
  
  xlab("Rank") + ylab("Frequency of Individual Occurence in Community")+  #axis labels
  
  scale_x_discrete(breaks=seq(0,131,10))  #defining intervals at which axis text will be visualised

comm_im_01_graph  # graph


comm_im1_avg.spp.graph <- cbind(comm_im1_pp95_avg.spp.graph,comm_im1_pp99_avg.spp.graph[1:133,2],comm_im1_pp1_avg.spp.graph[1:133,2])  #same here as above

comm_im_1_graph <- ggplot(comm_im1_avg.spp.graph, aes(x = reorder(rank, -Freq), group=1)) +
  
  geom_line(aes(y = Freq), size=1.5, colour ="black") +
  geom_line(aes(y = comm_im1_pp99_avg.spp.graph[1:133, 2]), size=1.5, colour ="blue") +
  geom_line(aes(y = comm_im1_pp1_avg.spp.graph[1:133, 2]), size=1.5, colour ="green") +
  
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 15))+
  
  xlab("Rank") + ylab("Frequency of Inividual Occurence in Community")+
  
  scale_x_discrete(breaks=seq(0,121,10))

comm_im_1_graph




####NOW TO DO EVERYTHING ABOVE BUT CONSIDERING A CLUMP AS AN INDIVIDUAL INSTEAD ####

clump_counter <- function(df){   #created clump counter function to take a community as an input and output a dataframe of all the species ID along with their clump frequencies
  
  freq_df <- as.data.frame(table(df)) # taking frequency of individuals as a dataframe from the table function
  freq_df[,1] <- as.numeric(as.character(freq_df[,1]))    #making sure the column with the species ID is a numeric (comes up as a factor initially so need to use as.character() to convert)
  freq_df$clump <- 0  # adding a new column for clump frequencies to the dataframe
  for(x in 1:nrow(df)){
    for(y in 1:ncol(df)){  #so checking each cell in the inputted matrix
      if(df[x,y]!=0){      #making sure the cell is not already 0
        for(i in 1:nrow(freq_df)){    #so for each cell in the matrix, going through the species IDs of the frequency data column
          if(freq_df[i,1]==df[x,y]){   # if it finds a species id that matches, 
            freq_df[i,3] <- freq_df[i,3]+1  # will add 1 to its clump frequency 
            spp <- df[x,y]   #species ID for input into clear clump
            df <-  ClearClump(x,y,df,spp)  # running clear clump floodfill function to then kill all connected cells with the same species ID as the clump that has just been counted. Ensures no recounting of clumps
          }
        }
        
      }
    }
  }
  return(freq_df)
}


Comm_1_clumps <- clump_counter(CommunityMatrix1)



### Now calculating it all over again for clumps!### The only thing different is the column and variable names, for notes on the actions taken below please consult the corresponding section for the individual calculations above.


comm_im1_pp95_1.freq.clump <- clump_counter(comm_im1_pp95_1)
comm_im1_pp95_1.sort.clump <- comm_im1_pp95_1.freq.clump[order(-comm_im1_pp95_1.freq.clump$clump),]
comm_im1_pp95_1.spp.clump <- subset(comm_im1_pp95_1.sort.clump, comm_im1_pp95_1.sort.clump$df!=0)
comm_im1_pp95_1.spp.clump$rank <- c(1:nrow(comm_im1_pp95_1.spp.clump))

comm_im1_pp95_2.freq.clump <- clump_counter(comm_im1_pp95_2)
comm_im1_pp95_2.sort.clump <- comm_im1_pp95_2.freq.clump[order(-comm_im1_pp95_2.freq.clump$clump),]
comm_im1_pp95_2.spp.clump <- subset(comm_im1_pp95_2.sort.clump, comm_im1_pp95_2.sort.clump$df!=0)
comm_im1_pp95_2.spp.clump$rank <- c(1:nrow(comm_im1_pp95_2.spp.clump))

comm_im1_pp95_3.freq.clump <- clump_counter(comm_im1_pp95_3)
comm_im1_pp95_3.sort.clump <- comm_im1_pp95_3.freq.clump[order(-comm_im1_pp95_3.freq.clump$clump),]
comm_im1_pp95_3.spp.clump <- subset(comm_im1_pp95_3.sort.clump, comm_im1_pp95_3.sort.clump$df!=0)
comm_im1_pp95_3.spp.clump$rank <- c(1:nrow(comm_im1_pp95_3.spp.clump))



comm_im1_pp99_1.freq.clump <- clump_counter(comm_im1_pp99_1)
comm_im1_pp99_1.sort.clump <- comm_im1_pp99_1.freq.clump[order(-comm_im1_pp99_1.freq.clump$clump),]
comm_im1_pp99_1.spp.clump <- subset(comm_im1_pp99_1.sort.clump, comm_im1_pp99_1.sort.clump$df!=0)
comm_im1_pp99_1.spp.clump$rank <- c(1:nrow(comm_im1_pp99_1.spp.clump))

comm_im1_pp99_2.freq.clump <- clump_counter(comm_im1_pp99_2)
comm_im1_pp99_2.sort.clump <- comm_im1_pp99_2.freq.clump[order(-comm_im1_pp99_2.freq.clump$clump),]
comm_im1_pp99_2.spp.clump <- subset(comm_im1_pp99_2.sort.clump, comm_im1_pp99_2.sort.clump$df!=0)
comm_im1_pp99_2.spp.clump$rank <- c(1:nrow(comm_im1_pp99_2.spp.clump))

comm_im1_pp99_3.freq.clump <- clump_counter(comm_im1_pp99_3)
comm_im1_pp99_3.sort.clump <- comm_im1_pp99_3.freq.clump[order(-comm_im1_pp99_3.freq.clump$clump),]
comm_im1_pp99_3.spp.clump <- subset(comm_im1_pp99_3.sort.clump, comm_im1_pp99_3.sort.clump$df!=0)
comm_im1_pp99_3.spp.clump$rank <- c(1:nrow(comm_im1_pp99_3.spp.clump))



comm_im1_pp1_1.freq.clump <- clump_counter(comm_im1_pp1_1)
comm_im1_pp1_1.sort.clump <- comm_im1_pp1_1.freq.clump[order(-comm_im1_pp1_1.freq.clump$clump),]
comm_im1_pp1_1.spp.clump <- subset(comm_im1_pp1_1.sort.clump, comm_im1_pp1_1.sort.clump$df!=0)
comm_im1_pp1_1.spp.clump$rank <- c(1:nrow(comm_im1_pp1_1.spp.clump))

comm_im1_pp1_2.freq.clump <- clump_counter(comm_im1_pp1_2)
comm_im1_pp1_2.sort.clump <- comm_im1_pp1_2.freq.clump[order(-comm_im1_pp1_2.freq.clump$clump),]
comm_im1_pp1_2.spp.clump <- subset(comm_im1_pp1_2.sort.clump, comm_im1_pp1_2.sort.clump$df!=0)
comm_im1_pp1_2.spp.clump$rank <- c(1:nrow(comm_im1_pp1_2.spp.clump))

comm_im1_pp1_3.freq.clump <- clump_counter(comm_im1_pp1_3)
comm_im1_pp1_3.sort.clump <- comm_im1_pp1_3.freq.clump[order(-comm_im1_pp1_3.freq.clump$clump),]
comm_im1_pp1_3.spp.clump <- subset(comm_im1_pp1_3.sort.clump, comm_im1_pp1_3.sort.clump$df!=0)
comm_im1_pp1_3.spp.clump$rank <- c(1:nrow(comm_im1_pp1_3.spp.clump))






comm_im01_pp95_1.freq.clump <- clump_counter(comm_im01_pp95_1)
comm_im01_pp95_1.sort.clump <- comm_im01_pp95_1.freq.clump[order(-comm_im01_pp95_1.freq.clump$clump),]
comm_im01_pp95_1.spp.clump <- subset(comm_im01_pp95_1.sort.clump, comm_im01_pp95_1.sort.clump$df!=0)
comm_im01_pp95_1.spp.clump$rank <- c(1:nrow(comm_im01_pp95_1.spp.clump))

comm_im01_pp95_2.freq.clump <- clump_counter(comm_im01_pp95_2)
comm_im01_pp95_2.sort.clump <- comm_im01_pp95_2.freq.clump[order(-comm_im01_pp95_2.freq.clump$clump),]
comm_im01_pp95_2.spp.clump <- subset(comm_im01_pp95_2.sort.clump, comm_im01_pp95_2.sort.clump$df!=0)
comm_im01_pp95_2.spp.clump$rank <- c(1:nrow(comm_im01_pp95_2.spp.clump))

comm_im01_pp95_3.freq.clump <- clump_counter(comm_im01_pp95_3)
comm_im01_pp95_3.sort.clump <- comm_im01_pp95_3.freq.clump[order(-comm_im01_pp95_3.freq.clump$clump),]
comm_im01_pp95_3.spp.clump <- subset(comm_im01_pp95_3.sort.clump, comm_im01_pp95_3.sort.clump$df!=0)
comm_im01_pp95_3.spp.clump$rank <- c(1:nrow(comm_im01_pp95_3.spp.clump))



comm_im01_pp99_1.freq.clump <- clump_counter(comm_im01_pp99_1)
comm_im01_pp99_1.sort.clump <- comm_im01_pp99_1.freq.clump[order(-comm_im01_pp99_1.freq.clump$clump),]
comm_im01_pp99_1.spp.clump <- subset(comm_im01_pp99_1.sort.clump, comm_im01_pp99_1.sort.clump$df!=0)
comm_im01_pp99_1.spp.clump$rank <- c(1:nrow(comm_im01_pp99_1.spp.clump))

comm_im01_pp99_2.freq.clump <- clump_counter(comm_im01_pp99_2)
comm_im01_pp99_2.sort.clump <- comm_im01_pp99_2.freq.clump[order(-comm_im01_pp99_2.freq.clump$clump),]
comm_im01_pp99_2.spp.clump <- subset(comm_im01_pp99_2.sort.clump, comm_im01_pp99_2.sort.clump$df!=0)
comm_im01_pp99_2.spp.clump$rank <- c(1:nrow(comm_im01_pp99_2.spp.clump))

comm_im01_pp99_3.freq.clump <- clump_counter(comm_im01_pp99_3)
comm_im01_pp99_3.sort.clump <- comm_im01_pp99_3.freq.clump[order(-comm_im01_pp99_3.freq.clump$clump),]
comm_im01_pp99_3.spp.clump <- subset(comm_im01_pp99_3.sort.clump, comm_im01_pp99_3.sort.clump$df!=0)
comm_im01_pp99_3.spp.clump$rank <- c(1:nrow(comm_im01_pp99_3.spp.clump))



comm_im01_pp1_1.freq.clump <- clump_counter(comm_im01_pp1_1)
comm_im01_pp1_1.sort.clump <- comm_im01_pp1_1.freq.clump[order(-comm_im01_pp1_1.freq.clump$clump),]
comm_im01_pp1_1.spp.clump <- subset(comm_im01_pp1_1.sort.clump, comm_im01_pp1_1.sort.clump$df!=0)
comm_im01_pp1_1.spp.clump$rank <- c(1:nrow(comm_im01_pp1_1.spp.clump))

comm_im01_pp1_2.freq.clump <- clump_counter(comm_im01_pp1_2)
comm_im01_pp1_2.sort.clump <- comm_im01_pp1_2.freq.clump[order(-comm_im01_pp1_2.freq.clump$clump),]
comm_im01_pp1_2.spp.clump <- subset(comm_im01_pp1_2.sort.clump, comm_im01_pp1_2.sort.clump$df!=0)
comm_im01_pp1_2.spp.clump$rank <- c(1:nrow(comm_im01_pp1_2.spp.clump))

comm_im01_pp1_3.freq.clump <- clump_counter(comm_im01_pp1_3)
comm_im01_pp1_3.sort.clump <- comm_im01_pp1_3.freq.clump[order(-comm_im01_pp1_3.freq.clump$clump),]
comm_im01_pp1_3.spp.clump <- subset(comm_im01_pp1_3.sort.clump, comm_im01_pp1_3.sort.clump$df!=0)
comm_im01_pp1_3.spp.clump$rank <- c(1:nrow(comm_im01_pp1_3.spp.clump))




### CLUMP EF VALUES###

eff_no_species.clump <- function(df){
  df.N <- sum(df$clump)
  df.prop <- df$clump/df.N
  df.mult <- df.prop*log(df.prop)
  df.shannon <- 0 - sum(df.mult)
  df.effective <- exp(df.shannon)
  return(df.effective)
  
}


comm_im1_pp95_1.EF.clump <- eff_no_species.clump(comm_im1_pp95_1.spp.clump)
comm_im1_pp95_2.EF.clump <- eff_no_species.clump(comm_im1_pp95_2.spp.clump)
comm_im1_pp95_3.EF.clump <- eff_no_species.clump(comm_im1_pp95_3.spp.clump)
comm_im1_pp95_avg.EF.clump <- (comm_im1_pp95_1.EF.clump+comm_im1_pp95_2.EF.clump+comm_im1_pp95_3.EF.clump)/3
comm_im1_pp95_avg.EF.clump.err <- std.error(c(comm_im1_pp95_1.EF.clump,comm_im1_pp95_2.EF.clump,comm_im1_pp95_3.EF.clump))


comm_im1_pp99_1.EF.clump <- eff_no_species.clump(comm_im1_pp99_1.spp.clump)
comm_im1_pp99_2.EF.clump <- eff_no_species.clump(comm_im1_pp99_2.spp.clump)
comm_im1_pp99_3.EF.clump <- eff_no_species.clump(comm_im1_pp99_3.spp.clump)
comm_im1_pp99_avg.EF.clump <- (comm_im1_pp99_1.EF.clump+comm_im1_pp99_2.EF.clump+comm_im1_pp99_3.EF.clump)/3
comm_im1_pp99_avg.EF.clump.err <- std.error(c(comm_im1_pp99_1.EF.clump,comm_im1_pp99_2.EF.clump,comm_im1_pp99_3.EF.clump))

comm_im1_pp1_1.EF.clump <- eff_no_species.clump(comm_im1_pp1_1.spp.clump)
comm_im1_pp1_2.EF.clump <- eff_no_species.clump(comm_im1_pp1_2.spp.clump)
comm_im1_pp1_3.EF.clump <- eff_no_species.clump(comm_im1_pp1_3.spp.clump)
comm_im1_pp1_avg.EF.clump <- (comm_im1_pp1_1.EF.clump+comm_im1_pp1_2.EF.clump+comm_im1_pp1_3.EF.clump)/3
comm_im1_pp1_avg.EF.clump.err <- std.error(c(comm_im1_pp1_1.EF.clump,comm_im1_pp1_2.EF.clump,comm_im1_pp1_3.EF.clump))



comm_im01_pp95_1.EF.clump <- eff_no_species.clump(comm_im01_pp95_1.spp.clump)
comm_im01_pp95_2.EF.clump <- eff_no_species.clump(comm_im01_pp95_2.spp.clump)
comm_im01_pp95_3.EF.clump <- eff_no_species.clump(comm_im01_pp95_3.spp.clump)
comm_im01_pp95_avg.EF.clump <- (comm_im01_pp95_1.EF.clump+comm_im01_pp95_2.EF.clump+comm_im01_pp95_3.EF.clump)/3
comm_im01_pp95_avg.EF.clump.err <- std.error(c(comm_im01_pp95_1.EF.clump,comm_im01_pp95_2.EF.clump,comm_im01_pp95_3.EF.clump))

comm_im01_pp99_1.EF.clump <- eff_no_species.clump(comm_im01_pp99_1.spp.clump)
comm_im01_pp99_2.EF.clump <- eff_no_species.clump(comm_im01_pp99_2.spp.clump)
comm_im01_pp99_3.EF.clump <- eff_no_species.clump(comm_im01_pp99_3.spp.clump)
comm_im01_pp99_avg.EF.clump <- (comm_im01_pp99_1.EF.clump+comm_im01_pp99_2.EF.clump+comm_im01_pp99_3.EF.clump)/3
comm_im01_pp99_avg.EF.clump.err <- std.error(c(comm_im01_pp99_1.EF.clump,comm_im01_pp99_2.EF.clump,comm_im01_pp99_3.EF.clump))

comm_im01_pp1_1.EF.clump <- eff_no_species.clump(comm_im01_pp1_1.spp.clump)
comm_im01_pp1_2.EF.clump <- eff_no_species.clump(comm_im01_pp1_2.spp.clump)
comm_im01_pp1_3.EF.clump <- eff_no_species.clump(comm_im01_pp1_3.spp.clump)
comm_im01_pp1_avg.EF.clump <- (comm_im01_pp1_1.EF.clump+comm_im01_pp1_2.EF.clump+comm_im01_pp1_3.EF.clump)/3
comm_im01_pp1_avg.EF.clump.err <- std.error(c(comm_im01_pp1_1.EF.clump,comm_im01_pp1_2.EF.clump,comm_im01_pp1_3.EF.clump))



EF.no.clump<- c(comm_im1_pp95_avg.EF.clump,comm_im1_pp99_avg.EF.clump,comm_im1_pp1_avg.EF.clump,comm_im01_pp95_avg.EF.clump,comm_im01_pp99_avg.EF.clump,comm_im01_pp1_avg.EF.clump)
EF.names.clump <- c("comm_im1_pp95_avg.EF.clump", "comm_im1_pp99_avg.EF.clump", "comm_im1_pp1_avgEF", "comm_im01_pp95_avg.EF.clump", "comm_im01_pp99_avg.EF.clump", "comm_im01_pp1_avg.EF.clump" )
EF.error.clump <- c(comm_im1_pp95_avg.EF.clump.err,comm_im1_pp99_avg.EF.clump.err,comm_im1_pp1_avg.EF.clump.err,comm_im01_pp95_avg.EF.clump.err,comm_im01_pp99_avg.EF.clump.err,comm_im01_pp1_avg.EF.clump.err)
EF.clump <- rbind(EF.names.clump, EF.no.clump, EF.error.clump )

View(EF.clump)




### BP IoD ###


BP_IoD.clump <- function(df){
  df.BP <- max(df$clump)/sum(df$clump)
  return(df.BP)
}

comm_im1_pp95_1.BP.clump <- BP_IoD.clump(comm_im1_pp95_1.spp.clump)
comm_im1_pp95_2.BP.clump <- BP_IoD.clump(comm_im1_pp95_2.spp.clump)
comm_im1_pp95_3.BP.clump <- BP_IoD.clump(comm_im1_pp95_3.spp.clump)
comm_im1_pp95_avg.BP.clump <- (comm_im1_pp95_1.BP.clump+comm_im1_pp95_2.BP.clump+comm_im1_pp95_3.BP.clump)/3
comm_im1_pp95_avg.BP.clump.err <- std.error(c(comm_im1_pp95_1.BP.clump,comm_im1_pp95_2.BP.clump,comm_im1_pp95_3.BP.clump))

comm_im1_pp99_1.BP.clump <- BP_IoD.clump(comm_im1_pp99_1.spp.clump)
comm_im1_pp99_2.BP.clump <- BP_IoD.clump(comm_im1_pp99_2.spp.clump)
comm_im1_pp99_3.BP.clump <- BP_IoD.clump(comm_im1_pp99_3.spp.clump)
comm_im1_pp99_avg.BP.clump <- (comm_im1_pp99_1.BP.clump+comm_im1_pp99_2.BP.clump+comm_im1_pp99_3.BP.clump)/3
comm_im1_pp99_avg.BP.clump.err <- std.error(c(comm_im1_pp99_1.BP.clump,comm_im1_pp99_2.BP.clump,comm_im1_pp99_3.BP.clump))


comm_im1_pp1_1.BP.clump <- BP_IoD.clump(comm_im1_pp1_1.spp.clump)
comm_im1_pp1_2.BP.clump <- BP_IoD.clump(comm_im1_pp1_2.spp.clump)
comm_im1_pp1_3.BP.clump <- BP_IoD.clump(comm_im1_pp1_3.spp.clump)
comm_im1_pp1_avg.BP.clump <- (comm_im1_pp1_1.BP.clump+comm_im1_pp1_2.BP.clump+comm_im1_pp1_3.BP.clump)/3
comm_im1_pp1_avg.BP.clump.err <- std.error(c(comm_im1_pp1_1.BP.clump,comm_im1_pp1_2.BP.clump,comm_im1_pp1_3.BP.clump))



comm_im01_pp95_1.BP.clump <- BP_IoD.clump(comm_im01_pp95_1.spp.clump)
comm_im01_pp95_2.BP.clump <- BP_IoD.clump(comm_im01_pp95_2.spp.clump)
comm_im01_pp95_3.BP.clump <- BP_IoD.clump(comm_im01_pp95_3.spp.clump)
comm_im01_pp95_avg.BP.clump <- (comm_im01_pp95_1.BP.clump+comm_im01_pp95_2.BP.clump+comm_im01_pp95_3.BP.clump)/3
comm_im01_pp95_avg.BP.clump.err <- std.error(c(comm_im01_pp95_1.BP.clump,comm_im01_pp95_2.BP.clump,comm_im01_pp95_3.BP.clump))

comm_im01_pp99_1.BP.clump <- BP_IoD.clump(comm_im01_pp99_1.spp.clump)
comm_im01_pp99_2.BP.clump <- BP_IoD.clump(comm_im01_pp99_2.spp.clump)
comm_im01_pp99_3.BP.clump <- BP_IoD.clump(comm_im01_pp99_3.spp.clump)
comm_im01_pp99_avg.BP.clump <- (comm_im01_pp99_1.BP.clump+comm_im01_pp99_2.BP.clump+comm_im01_pp99_3.BP.clump)/3
comm_im01_pp99_avg.BP.clump.err <- std.error(c(comm_im01_pp99_1.BP.clump,comm_im01_pp99_2.BP.clump,comm_im01_pp99_3.BP.clump))

comm_im01_pp1_1.BP.clump <- BP_IoD.clump(comm_im01_pp1_1.spp.clump)
comm_im01_pp1_2.BP.clump <- BP_IoD.clump(comm_im01_pp1_2.spp.clump)
comm_im01_pp1_3.BP.clump <- BP_IoD.clump(comm_im01_pp1_3.spp.clump)
comm_im01_pp1_avg.BP.clump <- (comm_im01_pp1_1.BP.clump+comm_im01_pp1_2.BP.clump+comm_im01_pp1_3.BP.clump)/3
comm_im01_pp1_avg.BP.clump.err <- std.error(c(comm_im01_pp1_1.BP.clump,comm_im01_pp1_2.BP.clump,comm_im01_pp1_3.BP.clump))



BP.no.clump<- c(comm_im1_pp95_avg.BP.clump,comm_im1_pp99_avg.BP.clump,comm_im1_pp1_avg.BP.clump,comm_im01_pp95_avg.BP.clump,comm_im01_pp99_avg.BP.clump,comm_im01_pp1_avg.BP.clump)
BP.names.clump <- c("comm_im1_pp95_avg.BP.clump", "comm_im1_pp99_avg.BP.clump", "comm_im1_pp1_avgBP", "comm_im01_pp95_avg.BP.clump", "comm_im01_pp99_avg.BP.clump", "comm_im01_pp1_avg.BP.clump" )
BP.error.clump <- c(comm_im1_pp95_avg.BP.clump.err,comm_im1_pp99_avg.BP.clump.err,comm_im1_pp1_avg.BP.clump.err,comm_im01_pp95_avg.BP.clump.err,comm_im01_pp99_avg.BP.clump.err,comm_im01_pp1_avg.BP.clump.err)
BP.clump <- rbind(BP.names.clump, BP.no.clump, BP.error.clump )

View(BP.clump)



### SADs ###



comm_im01_pp1_1.spp.clump.graph <- comm_im01_pp1_1.spp.clump[1:131,1:4]
comm_im01_pp1_2.spp.clump.graph <- comm_im01_pp1_2.spp.clump[1:131,1:4]
comm_im01_pp1_3.spp.clump.graph <- comm_im01_pp1_3.spp.clump[1:131,1:4]

comm_im01_pp1_avg.spp.clump.graph <- (comm_im01_pp1_1.spp.clump.graph+comm_im01_pp1_2.spp.clump.graph+comm_im01_pp1_3.spp.clump.graph)/3

plot(comm_im01_pp1_avg.spp.clump.graph$clump ~ comm_im01_pp1_avg.spp.clump.graph$rank, type = "l")



comm_im01_pp95_1.spp.clump.graph <- comm_im01_pp95_1.spp.clump[1:121,1:4]
comm_im01_pp95_2.spp.clump.graph <- comm_im01_pp95_2.spp.clump[1:121,1:4]
comm_im01_pp95_3.spp.clump.graph <- comm_im01_pp95_3.spp.clump[1:121,1:4]

comm_im01_pp95_avg.spp.clump.graph <- (comm_im01_pp95_1.spp.clump.graph+comm_im01_pp95_2.spp.clump.graph+comm_im01_pp95_3.spp.clump.graph)/3

plot(comm_im01_pp95_avg.spp.clump.graph$clump ~ comm_im01_pp95_avg.spp.clump.graph$rank, type = "l")




comm_im01_pp99_1.spp.clump.graph <- comm_im01_pp99_1.spp.clump[1:137,1:4]
comm_im01_pp99_2.spp.clump.graph <- comm_im01_pp99_2.spp.clump[1:137,1:4]
comm_im01_pp99_3.spp.clump.graph <- comm_im01_pp99_3.spp.clump[1:137,1:4]

comm_im01_pp99_avg.spp.clump.graph <- (comm_im01_pp99_1.spp.clump.graph+comm_im01_pp99_2.spp.clump.graph+comm_im01_pp99_3.spp.clump.graph)/3

plot(comm_im01_pp99_avg.spp.clump.graph$clump ~ comm_im01_pp99_avg.spp.clump.graph$rank, type = "l")





comm_im1_pp1_1.spp.clump.graph <- comm_im1_pp1_1.spp.clump[1:139,1:4]
comm_im1_pp1_2.spp.clump.graph <- comm_im1_pp1_2.spp.clump[1:139,1:4]
comm_im1_pp1_3.spp.clump.graph <- comm_im1_pp1_3.spp.clump[1:139,1:4]

comm_im1_pp1_avg.spp.clump.graph <- (comm_im1_pp1_1.spp.clump.graph+comm_im1_pp1_2.spp.clump.graph+comm_im1_pp1_3.spp.clump.graph)/3

plot(comm_im1_pp1_avg.spp.clump.graph$clump ~ comm_im1_pp1_avg.spp.clump.graph$rank, type = "l")


comm_im1_pp95_1.spp.clump.graph <- comm_im1_pp95_1.spp.clump[1:133,1:4]
comm_im1_pp95_2.spp.clump.graph <- comm_im1_pp95_2.spp.clump[1:133,1:4]
comm_im1_pp95_3.spp.clump.graph <- comm_im1_pp95_3.spp.clump[1:133,1:4]

comm_im1_pp95_avg.spp.clump.graph <- (comm_im1_pp95_1.spp.clump.graph+comm_im1_pp95_2.spp.clump.graph+comm_im1_pp95_3.spp.clump.graph)/3

plot(comm_im1_pp95_avg.spp.clump.graph$clump ~ comm_im1_pp95_avg.spp.clump.graph$rank, type = "l")




comm_im1_pp99_1.spp.clump.graph <- comm_im1_pp99_1.spp.clump[1:141,1:4]
comm_im1_pp99_2.spp.clump.graph <- comm_im1_pp99_2.spp.clump[1:141,1:4]
comm_im1_pp99_3.spp.clump.graph <- comm_im1_pp99_3.spp.clump[1:141,1:4]

comm_im1_pp99_avg.spp.clump.graph <- (comm_im1_pp99_1.spp.clump.graph+comm_im1_pp99_2.spp.clump.graph+comm_im1_pp99_3.spp.clump.graph)/3










comm_im01_avg.spp.clump.graph <- cbind(comm_im01_pp95_avg.spp.clump.graph,comm_im01_pp99_avg.spp.clump.graph[1:121,2],comm_im01_pp1_avg.spp.clump.graph[1:121,2])



comm_im_01_graph.clump <- ggplot(comm_im01_avg.spp.clump.graph, aes(x = reorder(rank, -clump), group=1)) +
  
  geom_line(aes(y = clump), size=1.5, colour ="black") +
  geom_line(aes(y = comm_im01_pp99_avg.spp.clump.graph[1:121, 3]), size=1.5, colour ="blue") +
  geom_line(aes(y = comm_im01_pp1_avg.spp.clump.graph[1:121, 3]), size=1.5, colour ="green") +
  
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 15))+
  
  xlab("Rank") + ylab("Frequency of Meta-Individual Occurrence in Community")+
  
  scale_x_discrete(breaks=seq(0,131,10))

comm_im_01_graph.clump


comm_im1_avg.spp.clump.graph <- cbind(comm_im1_pp95_avg.spp.clump.graph,comm_im1_pp99_avg.spp.clump.graph[1:133,2],comm_im1_pp1_avg.spp.clump.graph[1:133,2])

comm_im_1_graph.clump <- ggplot(comm_im1_avg.spp.clump.graph, aes(x = reorder(rank, -clump), group=1)) +
  
  geom_line(aes(y = clump), size=1.5, colour ="black") +
  geom_line(aes(y = comm_im1_pp99_avg.spp.clump.graph[1:133, 3]), size=1.5, colour ="blue") +
  geom_line(aes(y = comm_im1_pp1_avg.spp.clump.graph[1:133, 3]), size=1.5, colour ="green") +
  
  theme(axis.text = element_text(size = 12, angle = 90), axis.title = element_text(size = 15))+
  
  xlab("Rank") + ylab("Frequency of Meta-Individual Occurence in Community")+
  
  scale_x_discrete(breaks=seq(0,121,10))

comm_im_1_graph.clump



### LANDSCAPE VISUALISATION ###
plot(raster(CommunityMatrix1))

w <- plot(raster(comm_im01_pp1_1))
e <- plot(raster(comm_im01_pp99_1))
r <- plot(raster(comm_im01_pp95_1))
t <- plot(raster(comm_im1_pp95_1))

# so compare between pp1 and pp95 for im1, then compare for im01pp95 and im1pp95

