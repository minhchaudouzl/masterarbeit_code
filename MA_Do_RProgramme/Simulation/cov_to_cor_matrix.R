library(stats)
library(rlist)

mypath <- "D:/MA_Do_RProgramme"

mvn_list_gleason <- list.load(paste0(mypath, "/Simulation/Gleason/mvn_list_gleason.rds"))
mvn_list_jorgensen <- list.load(paste0(mypath, "/Simulation/Jorgensen/mvn_list_jorgensen.rds"))
mvn_list_leon <- list.load(paste0(mypath, "/Simulation/Leon/mvn_list_leon.rds"))
mvn_list_mack <- list.load(paste0(mypath, "/Simulation/Mack/mvn_list_mack.rds"))
mvn_list_makkar2012 <- list.load(paste0(mypath, "/Simulation/Makkar2012/mvn_list_makkar2012.rds"))
mvn_list_makkar2020 <- list.load(paste0(mypath, "/Simulation/Makkar2020/mvn_list_makkar2020.rds"))
mvn_list_popma <- list.load(paste0(mypath, "/Simulation/Popma/mvn_list_popma.rds"))
mvn_list_thyregod <- list.load(paste0(mypath, "/Simulation/Thyregod/mvn_list_thyregod.rds"))
mvn_list_all <- list.load(paste0(mypath, "/Simulation/all/mvn_list_all.rds"))

cor_matrix <- function(mvn_list){
  cov_0 <- mvn_list$Sigma_0
  cov_1 <- mvn_list$Sigma_1
  
  cor_0 <- cov2cor(cov_0)
  cor_1 <- cov2cor(cov_1)
  
  return(list(cor_0 = cor_0, cor_1 = cor_1))
}

cor_gleason <- cor_matrix(mvn_list_gleason)
cor_jorgensen <- cor_matrix(mvn_list_jorgensen)
cor_leon <- cor_matrix(mvn_list_leon)
cor_mack <- cor_matrix(mvn_list_mack)
cor_makkar2012 <- cor_matrix(mvn_list_makkar2012)
cor_makkar2020 <- cor_matrix(mvn_list_makkar2020)
cor_popma <- cor_matrix(mvn_list_popma)
cor_thyregod <- cor_matrix(mvn_list_thyregod)
cor_all <- cor_matrix(mvn_list_all)


mean_vector <- function(mvn_list){
  print("SAVR:")
  print(mvn_list$mu_0)
  print("TAVR:")
  print(mvn_list$mu_1)
}
