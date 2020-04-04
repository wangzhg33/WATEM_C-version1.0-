### scenario with no AD and erosion (one run)

rm(list=ls())
graphics.off()

options(scipen = 999)

###### Set root folder

scena_no<-1

current_dir <-getwd()

setwd(paste(current_dir,'Catchments',sep='/'))

#setwd("F:/Geoscientific model development/WATEM_C/test/Catchments")

library(raster)
library(tools)
library(ggplot2)


input_dir <- "input_data"
work_dir <- "data"

###### LAND USE
parcel_fileName <- "PARCEL.tif"
PARCEL <- raster(paste(input_dir,parcel_fileName,sep='/'))
writeRaster(PARCEL, filename =paste0(work_dir,'/LAND_USE.RST'),            
            format ='IDRISI',NAflag =0, datatype ="INT2S",overwrite =TRUE)

######## C Factor
Cfactor_fileName <- "C_factor.tif"
C <- raster(paste(input_dir,Cfactor_fileName,sep='/'))
writeRaster(C, filename = paste0(work_dir,'/C_factor.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE) 

######## K Factor
Kfactor_fileName <- "K_factor.tif"
K <- raster(paste(input_dir,Kfactor_fileName,sep='/'))
writeRaster(K, filename = paste0(work_dir,'/K_factor.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE)

######## P Factor
Pfactor_fileName <- "P_factor.tif"
P <- raster(paste(input_dir,Pfactor_fileName,sep='/'))
writeRaster(P, filename = paste0(work_dir,'/P_factor.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE) 

######## R Factor
Rfactor_fileName <- "R_factor.tif"
R <- raster(paste(input_dir,Rfactor_fileName,sep='/'))
writeRaster(R, filename = paste0(work_dir,'/R_factor.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE) 

######## RKP multiplication
RKP<-R*K*P/1000

writeRaster(RKP, filename = paste0(work_dir,'/RKP.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE) 

# DEM
DEM_fileName <- "DEM.tif"
DEM <- raster(paste(input_dir,DEM_fileName,sep='/'))
writeRaster(DEM, filename = paste0(work_dir,'/DEM.RST'),
            format = "IDRISI",NAflag = 0,datatype = "FLT4S",overwrite = TRUE) 

ktc_fileName <- "ktc.tif"
ktc <- raster(paste(input_dir,ktc_fileName,sep='/'))
writeRaster(ktc, filename = paste0(work_dir,"/KTC.rst"),
            format = "IDRISI",NAflag = 0,datatype = "INT2S",overwrite = TRUE) 


### close erosion by setting start year and end year
### close diffusion and advection by setting K0 and v0


erosion_start_year <- 2100
erosion_end_year <- 2100

## variables for C cycling
depth_interval <- 5
depth <- 100

tillage_depth <- 5

deltaC13_ini_top <- -27.0
deltaC13_ini_bot <- -27.5
deltaC14_ini_top <- 50.0
deltaC14_ini_bot <- -500.0

time_equilibrium <- 100000
k1 <- 2.1
k2 <- 0.03
k3 <- 0.002
hAS <- 0.12
hAP <- 0.01
hSP <- 0.01
r0 <- 1
C_input <- 2.0 ##input from root
C_input2 <- 0.5## input from manure, residue
r_exp <- 3.30
i_exp <- 6
C13_discri <- 0.9965
C14_discri <- 0.996
deltaC13_input_default <- -29.0
deltaC14_input_default <- 100.0
Cs137_input_default <- 0.0

Sand_ini_top <- 15
Silt_ini_top <- 70
Clay_ini_top <- 15
Sand_ini_bot <- 15
Silt_ini_bot <- 70
Clay_ini_bot <- 15

K0 <- 0.0
Kfzp <- 0.00001
v0 <- 0.0
vfzp <- 0.00001

a_erer <- 1.0
b_erero <- 2000.0
b_erdepo <- -4000.0

time_step <- 5

unstable <- FALSE              

col_num<-0;
sequence <-seq(0.002,0.002, by=0.001)


for (temp_data in sequence) {  
  col_num=col_num+1
  
  cmd <- paste0("WATEM_C.exe -dir ",work_dir,
                " -o ","output",
                " -d ","DEM.rst",
                " -p ","LAND_USE.rst",
                " -u ","RKP.rst",
                " -c ","C_factor.rst",
                " -k ","KTC.rst",
                " -MC_f 1 -t 100 -r mf -b 1400 -ktil 0 -pdep 25 -frac 0 -w 1 -w_d_exp 0",
                " -C_input ", toString(format(C_input,nsmall=3)),
                " -C_input2 ", toString(format(C_input2,nsmall=3)),
#                " -deltaC13_ini_top ",toString(format(deltaC13_ini_top,nsmall=3)),
#                " -deltaC13_ini_bot ", toString(format(deltaC13_ini_bot,nsmall=3)),

                " -depth ", toString(depth),
                " -depth_interval ",toString(depth_interval),
                " -tillage_depth ",toString(tillage_depth),
                " -time_equilibrium ", toString(time_equilibrium),
                " -erosion_start_year ", toString(erosion_start_year),
                " -erosion_end_year ",toString(erosion_end_year),
                " -k1 ",toString(format(k1,nsmall=3))," -k2 ",toString(format(k2,nsmall=3)),
                " -k3 ",toString(format(temp_data,nsmall=3)),
                " -hAS ",toString(format(hAS,nsmall=3))," -hAP ",toString(format(hAP,nsmall=3)),
                " -hSP ",toString(format(hSP,nsmall=3)), 
                
                " -r_exp ",toString(format(r_exp,nsmall=3))," -i_exp ",toString(format(i_exp,nsmall=3)),
                " -a_erer ",toString(format(a_erer,nsmall=3))," -b_erero ",toString(format(b_erero,nsmall=3)),
                #                " -b_erdepo ",toString(format(b_erdepo,nsmall=3)),                  
                " -K0 ",toString(format(K0,nsmall=3))," -Kfzp ",toString(format(Kfzp,nsmall=3)),
                " -v0 ",toString(format(v0,nsmall=3))," -vfzp ",toString(format(vfzp,nsmall=3)),                
                " -C13_discri ",toString(format(C13_discri,nsmall=3)),
                " -C14_discri ",toString(format(C14_discri,nsmall=3)),
                #               " -deltaC13_input_default ",toString(format(deltaC13_input_default,nsmall=3)),
                " -deltaC14_input_default ",toString(format(deltaC14_input_default,nsmall=3)),
                " -Cs137_input_default ",toString(format(Cs137_input_default,nsmall=3)),
                " -time_step ",toString(time_step), 
                " -c14input ","C14input.txt",
                " -c13input ","C13input.txt",
                " -Cs137input ","Cs137_input.txt")
  

  # run command
  system(cmd)
  
 
  
}


setwd("..")