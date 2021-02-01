# Plot 98 ----------------------------------------------
rm(list=ls())
date <- Sys.Date()
No.test<- "1_1"
lapply(c("ggplot2","mvnfast"), require, character.only = T)
# observed diameter distributions by plot, use observed.DBH[[98]]
load("observed.DBH_site_tree.RData"); source("generate_pars.R")
# mylist<-list()
# setwd("C:/Users/Alicia/Documents/GitHub/FL_Carbon")
envdata<-read.csv("slash_env_data_site_index2.csv", header=T, sep=",")

# colnames(envdata)
# set parameters for growth equation
SI<-envdata[5,16]
aridity<-envdata[5,15]*0.0001
temp<-envdata[5,10]
temp2<-envdata[5,10]^2
C_AVG<-((envdata[5,11]/10)*(1/3))+((envdata[5,12]/10)*(2/3))
N_AVG<-((envdata[5,13]/100)*(1/3))+((envdata[5,14]/100)*(2/3))
CN<-C_AVG/N_AVG
CN_SCALE<-(CN*4.941)+29.777

# set stand age and density
plot_density<-round(envdata[5,2])
observed.a<-round((10^(2.475538 - 0.3158355*temp + 0.007800543*temp2 + 0.5217529*aridity))*(envdata[5,3]^1.728844))
predict.tasb<-c() #total above stump biomass
predict.d<-c() #store average diameter for each simulation run
initial<-(10^(-1.431904 + 0.182686*temp - 0.004512*temp2 - 0.301793*aridity))*(0.5784213*1^(-0.4215787)) #initial diameter calculated from growth at year 1
initial_BA<-(pi*(initial*2.47)^2)/40000 #basal area per tree per hectare at year 1
initial_height<-10^(0.8730716)+initial^(0.8976928) #height of trees at year 1

# initialize the diameter for the first year
Diameter<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize diameter array
age<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize age array
TASB<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize the total above-stump biomass array
height<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize height (ft) array
BA<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize basal area (m^2/hectare)
RBA<-matrix(0,ncol = 100, nrow = 100*observed.a) # initialize relative basal area array
tbaha<-vector() # initialize vector for the total basal area per hectare
prob_denom.keep<-matrix(ncol=plot_density, nrow=(observed.a-1)) # initialize vector for prob of mortaility
simul.tree.d.last <- c()

#numbering rows and columns in array
Diameter[1:100,1:100]<-seq(1:10000)
BA[1:100,1:100]<-seq(1:10000)
height[1:100,1:100]<-seq(1:10000)

#selecting random numbers of array to place trees
trees <- sample.int(10000, plot_density, replace = F)

#load year 1 diameters where trees were randomly placed
Diameter <- replace(Diameter, Diameter %in% trees, initial)
Diameter <- replace(Diameter, Diameter!= initial, NA)

#load year 1 basal area where trees were randomly placed in array
BA <- replace(BA, BA %in% trees, initial_BA)
BA <- replace(BA, BA!=initial_BA, NA)

#load year 1 height where trees were randomly placed in array
height <- replace(height, height %in% trees, initial_height)
height <- replace(height, height!=initial_height, NA)
tbaha[1]<- initial_BA*plot_density

#index for trees
tree.index <- c()

for (i in 1:100) {
    tree.index1 <- cbind.data.frame("row.index"=na.omit(as.numeric(ifelse(!is.na(Diameter[i,]),i,NA))),
                                    "col.index"=na.omit(as.numeric(ifelse(!is.na(Diameter[i,]),1:100,NA))))
    tree.index <- rbind(tree.index,tree.index1)
  }
rm("tree.index1")

#introduce variable used for mortality
M <- numeric(length = 1)

#data assimilation
pmin <- c(); pmax <- c() #create empty vectors

#set prior ranges for seven parameters
#   a             b1              b2               b3           b4               b5               b6
pmin[1] <- -1.1;pmin[2] <- -0.05;pmin[3] <- -0.42;pmin[4] <- 2;pmin[5] <- -0.20;pmin[6] <- -0.12;pmin[7] <- 0.007;
pmax[1] <-  9.4;pmax[2] <-  0.45;pmax[3] <-  0.10;pmax[4] <-20;pmax[5] <- -0.03;pmax[6] <-  0.32;pmax[7] <- 0.101;

par.name <- c("a","b1","b2","b3","b4","b5","b6")
p_op <- c(4.2899, 0.1022, -0.2146, 11.54, -0.0699, -0.1112, 0.0079)
no.simu <- 100000
d <- 6
J_last <- 1000
updated <- 0

#
p_rec <- matrix(NA,length(pmin), no.simu)
p_upgraded <- matrix(NA, length(pmin), no.simu)
J_keep <- rep(NA,no.simu)
J <- c()
DJ <- c()
J_new1 <- c()
sd <- 2.38/sqrt(length(pmin))

for (simu in 1:no.simu) {
  if (simu <= 2000) { #two steps; 1-default sampling; 2-sampling from covariances of each pars
    pnew <- generate.pars(p_op, pmin, pmax, d)
  } else {
    pnew <- generate.pars.cov(p_op, pmin, pmax, covars)
  }
  
  p_rec[,simu] <- pnew #save pnew
  
  #mean of posterior parameter distribution
  # pnew <- apply(p_upgraded[,1:updated], 1, mean)

  #assign pars
  for (i in 1:length(par.name)) {
    assign(par.name[i], pnew[i])
  }

    # system.time({
    for (h in 1:(observed.a-1)){ #run simulation for h years
      k<-0
      
      for (l in 1:plot_density) { 
        i <- tree.index[l,"row.index"] # rows of array
        j <- tree.index[l,"col.index"] # columns of array
        
        Diameter[(i+h*100),j]<-Diameter[(i+(h-1)*100),j] #mortality calculations depend on previous diameter
        k<-k+1
        age[(i+h*100),j]<-age[(i+(h-1)*100),j]+1
        growth<-(10^(-1.431904 + 0.182686*temp - 0.004512*temp2 - 0.301793*aridity))*(0.5784213*age[(i+h*100),j]^(-0.4215787))
        BA[(i+h*100),j]<-BA[(i+(h-1)*100),j]
        tbaha[h+1]<-sum(BA[(1+(h-1)*100):(h*100),], na.rm=TRUE)
        height[(i+h*100),j]<- height[(i+(h-1)*100),j]
        
        # create 11x11 grid around tree, collect diameters within grid that are greater than or equal to ith jth tree
        # select 11x11 grid around tree #it's ok to remove exeeded grid as the values assigned as 0 
        guide.grid.v <- c((i-5:1),i,(i+1:5))
        guide.grid.h <- c((j-5:1),j,(j+1:5))
        guide.grid.v <- guide.grid.v[guide.grid.v > 0 & guide.grid.v < 100]
        guide.grid.h <- guide.grid.h[guide.grid.h > 0 & guide.grid.h < 100]
        
        all.grid <- BA[(guide.grid.v+h*100), guide.grid.h]
        all.grid.2 <- sum(ifelse(all.grid>=BA[(i+h*100),j], all.grid, 0), na.rm = T)

        RBA[(i+h*100),j] <- BA[(i+h*100),j]/(1/all.grid.2) #relative basal area calculated
        
        prob_denom<-(1+exp( a
                           +b1*(Diameter[(i+h*100),j]*2.54)
                           +b2*(height[(i+h*100),j]/3.281)
                           +b3*RBA[(i+h*100),j]
                           +b4*SI
                           +b5*tbaha[h+1]
                           +b6*(height[(i+h*100),j]/3.281)*tbaha[h+1]))

        M <- rbinom(1, 1,(0.125/prob_denom)) #eight year mortality probability converted to yearly mortality probability
        # M <- rbinom(1, 1,prob_denom.test[h]/2.2) # test
        prob_denom.keep[h,l] <- (.125/prob_denom) #save prob of mortality
        # Calculate the diameter for ith jth tree for h year
        Diameter[(i+h*100),j]<-Diameter[(i+(h-1)*100),j] + growth - M*(Diameter[(i+(h-1)*100),j]+growth) #if M=1, tree dies, diameter=0
        BA[(i+h*100),j]<-(pi*(Diameter[(i+h*100),j]*2.54)^2)/40000
        tbaha[h+1]<-sum(BA[(1+h*100):((h+1)*100),], na.rm=TRUE)
        height[(i+h*100),j]<-10^(0.8730716)+Diameter[(i+h*100),j]^(0.8976928)
        
        # use diameter and age to calculate total aboveground biomass of the ith jth tree in the h year
        TASB[(i+h*100),j]<-(0.041281*((Diameter[(i+h*100),j]*2.54)^2.722214))
        
        # If the tree dies, plant a new tree (age = 0)
        if (M==1){
          age[(i+h*100),j]<-0
        }
      }
      # print(paste("time step =", (h+1)))
    }
    # })
    # # save average modeled diameter
    # predict.d[o]<-mean(Diameter[(1+h*100):((h+1)*100),], na.rm = TRUE)
    # predict.tasb[o]<-sum(TASB[(1+h*100):((h+1)*100),], na.rm = TRUE)*.5*1000*(1/10000)

#simulated diameter
for (i in 1:plot_density) {
  simul.tree.d.last[i] <- Diameter[(tree.index[i,1]+h*100),tree.index[i,2]]
}


#obs
obs.d <- observed.DBH$plot_98
breaks.his <- seq(min(c(simul.tree.d.last,obs.d)),max(c(simul.tree.d.last,obs.d)), length.out = 15) 
obs.count <- hist(obs.d, breaks.his,plot = F)$counts
simul.count <- hist(simul.tree.d.last, breaks.his, plot = F)$counts

J <- (sort(simul.tree.d.last) - sort(obs.d))^2

DJ <- 2*(1*obs.d)^2
DJ[DJ == 0] <- DJ[DJ > 0][1] 
J_new <- J/DJ
delta.J <- sum(J_new) - J_last

if (min(1, exp(-delta.J)) > runif(1)) {
  p_op <- pnew
  J_last <- sum(J_new)
  updated <- updated + 1
  J_keep[updated] <- sum(J_new) 
  p_upgraded[,updated] <- p_op
  # plot(prob_denom.keep)
  if (updated %in% c(100*1:100)) {
  par(mfrow=c(2,7))
  par(mar=c(5,5,1,1))
  for (i in 1:7) {
    hist(p_upgraded[i,(updated/2):updated],xlim = c(pmin[i],pmax[i]), main = par.name[i], xlab =NA)
  }
  hist(obs.d, breaks.his, main = "Observed diameter", xlab="Diameter (in)")
  hist(simul.tree.d.last, breaks.his, main = "Simulated diameter", xlab="Diameter (in)")
  plot(sort(simul.tree.d.last), sort(obs.d), xlab = "observed DBH", ylab = "simulated DBH", xlim = c(4.5,16), ylim = c(4.5,16))
  abline(0,1)
  }
}

if (simu == 2000) {
  covars <- cov(t(p_rec[,1:simu]))
}
if (simu > 2000) {
  covars <- sd*cov(t(p_rec[,1:simu]))
}
 print(paste("simu =", simu, "updated =", updated))
}

#check posterior distribution
png(paste(date,No.test,".png"), width = 2000, height = 300, res=200)
par(mfrow=c(1,7))
par(mar=c(2,2,1.7,1.7))
for (i in 1:7) {
  hist(p_upgraded[i,(updated/2):updated],xlim = c(pmin[i],pmax[i]), main = par.name[i])
}
dev.off()

#histogram with validations
png(paste(date,No.test,"hist_val_sorted.png"), width = 3000, height = 900, res=200)
par(mfrow=c(2,7))
par(mar=c(5,5,1,1))
for (i in 1:7) {
  hist(p_upgraded[i,(updated/2):updated],xlim = c(pmin[i],pmax[i]), main = par.name[i], xlab =NA)
}
hist(obs.d, breaks.his, main = "Observed diameter", xlab="Diameter (in)")
hist(simul.tree.d.last, breaks.his, main = "Simulated diameter", xlab="Diameter (in)")
plot(sort(simul.tree.d.last), sort(obs.d), xlab = "observed DBH", ylab = "simulated DBH", xlim = c(4.5,16), ylim = c(4.5,16))
abline(0,1)
dev.off()

save.image(paste(date,No.test, "2.RData", sep = "_"))

# observed.d<-c(envdata[5,4])
# observed.tasb<-c(envdata[5,5])
# modeled.d<-mean(predict.d)
# modeled.tasb<-mean(predict.tasb)
# sd.tasb<-sd(predict.tasb)
# sd.dbh<-sd(predict.d)

# set up dataframe to store simulated data
# df<-cbind(observed.d, modeled.d, observed.a, plot_density, temp, temp2, aridity, CN_SCALE, observed.tasb, modeled.tasb, 
#           sd.tasb, sd.dbh)

# df<-data.frame(df)

# colnames(df)<-c("Observed_Diameter","Modeled_Diameter","Age", "Tree_Density", "Temperature", "Temp2", "Aridity", "Soil_CN", "Observed_Biomass", 
#                 "Modeled_Biomass", "sd.tasb", "sd.dbh")
# 
# # data will be saved as list 5
# mylist[[5]] <- df