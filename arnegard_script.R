# Arnegard paper analysis

library(car)
library(geomorph)
library(ggplot2)
library(MuMIn)
library(psych)
library(rgl)
library(scatterplot3d) 

arnegard.all <- read.csv('data-raw/arnegard_data.csv')

### Drop the non-experiment fish
arnegard <- arnegard[-c(637:713),]
arnegard.noA <- subset(arnegard, Fifteen_Percent_Loess_Curved_Landscape_Group == (c("B", "L", "O")))
arnegard.L <- subset(arnegard, Fifteen_Percent_Loess_Curved_Landscape_Group == "L")
arnegard.wild <- arnegard.all[c(674:713),]



ng1 <- theme(aspect.ratio=1.0,panel.background = element_blank(), 
             panel.grid.major = element_blank(), 
             panel.grid.minor = element_blank(),
             panel.border=element_blank(),
             axis.line = element_line(size=1), 
             axis.line.x = element_line(color="black", size = 1),
             axis.line.y = element_line(color="black", size = 1),
             axis.ticks=element_line(color="black"), 
             axis.text=element_text(color="black"), 
             axis.title=element_text(color="black"), 
             axis.title.y=element_text(vjust=0.2, size=12),
             axis.title.x=element_text(vjust=0.1,size=12),
             axis.text.x=element_text(size=10),
             axis.text.y=element_text(size=10),
             legend.position="none",
             legend.title=element_blank(),
             plot.title = element_blank())

#Isotope_PC1 is niche score; PC2 is deviation score
#Use dredge to figure out what explains niche score
#Run test model
#
multvardata <- na.omit((arnegard.BL[, (colnames(arnegard.BL) %in% c("PC1.abs", "PC2.abs", "Loess_Predicted_SL_mm", "Opening_Inlever.res", "Protrusion.res", "Buccal_Length.res", "Total_Long_Gill_Rakers", "Total_Short_Gill_Rakers", "Total_All_Gill_Rakers", "Gape.res", "Suction_Index", "Isotope_PC1"))]))

biglm <- lm(Loess_Predicted_SL_mm ~ PC1.abs + PC2.abs + Opening_Inlever.res + Protrusion.res + Buccal_Length.res + Total_Long_Gill_Rakers + Total_Short_Gill_Rakers + Total_All_Gill_Rakers + Gape.res, data = arnegard)
summary(biglm)


#Create "morphology" PCA
#Identify terms for PCA
arnegard.morphology.terms <- c("Opening_Inlever.res", "Protrusion.res", "Buccal_Length.res", "Total_Long_Gill_Rakers", "Total_Short_Gill_Rakers", "Gape.res", "Suction_Index", "EpaxH.res")
#Run PCA using
arnegard.morph.PCA <- principal(arnegard.noA[arnegard.morphology.terms], nfactors = 3, residuals = FALSE, rotate = "varimax", scores = TRUE, missing = FALSE)
#See eigenvalues
arnegard.morph.PCA$values
#Create a 2xn matrix of PCA scores
arnegard.pca <- arnegard.morph.PCA$scores
#Attach matrix to main dataframe
arnegard.noA <- cbind(arnegard.noA, arnegard.pca)

#Use "predict" R command to figure out where everone falls along a benthic/limnetic axis. can do for LDA and also PCA.

ggplot(arnegard, aes(x=Isotope_PC1, y=Isotope_PC, colour = Fifteen_Percent_Loess_Curved_Landscape_Group)) +
  geom_point() +
  stat_ellipse()

#Rc1 & Rc3 tell them apart...
TukeyHSD(aov(RC1 ~ Fifteen_Percent_Loess_Curved_Landscape_Group, data = arnegard.noA))

#What PC separates benthic/limnetic?
plot(arnegard.w.pca$Fifteen_Percent_Loess_Curved_Landscape_Group, arnegard.w.pca$RC3)
#RC3
summary(aov(RC3 ~ Fifteen_Percent_Loess_Curved_Landscape_Group, arnegard.w.pca))


#Trying a DFA
#Can i code A and O as 'O'"

z <- lda(Fifteen_Percent_Loess_Curved_Landscape_Group ~ Opening_Inlever.res + Protrusion.res + Buccal_Length.res + Total_Long_Gill_Rakers + Total_Short_Gill_Rakers + Total_All_Gill_Rakers + Gape.res, data = arnegard)

plot(z, group=arnegard$Fifteen_Percent_Loess_Curved_Landscape_Group)


###Some figures
PC1.Fig <- ggplot(arnegard.w.pca, aes(x = RC1, y = Isotope_PC1)) +
  labs(x = "MorphPC1", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC1.Fig

PC3.Fig <- ggplot(arnegard.w.pca, aes(x = RC3, y = Isotope_PC1)) +
  labs(x = "MorphPC2", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC3.Fig

ShortGill.Fig <- ggplot(arnegard.BL, aes(x = Total_Short_Gill_Rakers, y = Isotope_PC1)) +
  labs(x = "# short gill rakers", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
ShortGill.Fig

LongGill.Fig <- ggplot(arnegard.BL, aes(x = Total_Long_Gill_Rakers, y = Isotope_PC1)) +
  labs(x = "# long gill rakers", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
LongGill.Fig

BuccalLength.Fig <- ggplot(arnegard, aes(x = Buccal_Length.res, y = Loess_Predicted_SL_mm)) +
  labs(x = "Buccal Length (residual)", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
BuccalLength.Fig

Suction_Index.Fig <- ggplot(arnegard.BL, aes(x = Suction_Index, y = Isotope_PC1)) +
  labs(x = "Suction Index", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Suction_Index.Fig

Opening_Inlever.Fig <- ggplot(arnegard.BL, aes(x = scale(Opening_Inlever.res), y = Isotope_PC1)) +
  labs(x = "Opening Inlever", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Opening_Inlever.Fig

Protrusion.Fig <- ggplot(arnegard.BL, aes(x = Protrusion.res, y = Isotope_PC1)) +
  labs(x = "Protrusion", y = "Niche score") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Protrusion.Fig


###Trying fitness
RC3Fit.Fig <- ggplot(arnegard.w.pca, aes(x = RC3, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "PC3", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
RC3Fit.Fig

LongGill.Fig <- ggplot(arnegard.BL, aes(x = Total_Long_Gill_Rakers, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "# long gill rakers", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
LongGill.Fig

BuccalLength.Fig <- ggplot(arnegard.BL, aes(x = Buccal_Length.res, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "Buccal Length (residual)", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
BuccalLength.Fig

Suction_Index.Fig <- ggplot(arnegard.BL, aes(x = Suction_Index, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "Suction Index", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Suction_Index.Fig

Opening_Inlever.Fig <- ggplot(arnegard.BL, aes(x = Opening_Inlever.res, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "Opening Inlever", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Opening_Inlever.Fig

Protrusion.Fig <- ggplot(arnegard.BL, aes(x = Protrusion.res, y = Std_Length_TPS1_to_TPS14_mm)) +
  labs(x = "Protrusion", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
Protrusion.Fig


DFE_Fig1A <- ggplot(arnegard.BL, aes(x = Opening_Inlever.res, y = Isotope_PC1, group = Fifteen_Percent_Loess_Curved_Landscape_Group)) +
  labs(x = "Inflorescence number
       (standardized residual)", y = "Relative fitness (residual)") +
  geom_point(aes(shape = factor(Fifteen_Percent_Loess_Curved_Landscape_Group), fill = factor(Fifteen_Percent_Loess_Curved_Landscape_Group))) +
  scale_shape_manual(values=c(24, 19)) +
  scale_fill_manual(values=c("white", "black")) +
  geom_smooth(colour = "black", size = 1.05, fullrange = T, aes(linetype = factor(Fifteen_Percent_Loess_Curved_Landscape_Group)))  +
  scale_linetype_manual(values=c(5, 1)) +
  ng1
DFE_Fig1A


###Testing accumulation of RI w/ deviation on PC1 vs. pc2
arnegard$PC1.abs <- abs(arnegard$Isotope_PC1)
arnegard$PC2.abs <- abs(arnegard$Isotope_PC2)

plot(arnegard$PC1.abs, arnegard$PC2.abs)

#Loess
PC1.fig <- ggplot(arnegard, aes(x = PC1.abs, y = Loess_Predicted_SL_mm)) +
  labs(x = "PC1 (absolute)", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC1.fig

PC2.fig <- ggplot(arnegard, aes(x = PC2.abs, y = Loess_Predicted_SL_mm)) +
  labs(x = "PC2 (absolute)", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC2.fig

PC1.fig <- ggplot(arnegard, aes(x = PC1.abs, y = SL_mm_Scaled_0_to_1)) +
  labs(x = "PC1 (absolute)", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC1.fig

PC2.fig <- ggplot(arnegard, aes(x = PC2.abs, y = SL_mm_Scaled_0_to_1)) +
  labs(x = "PC2 (absolute)", y = "Body Length (mm)") +
  geom_point(aes()) +
  geom_smooth(se = T, colour = "black", size = 1.05, fullrange = T)  +
  ng1
PC2.fig

PC1.lin <- lm(Loess_Predicted_SL_mm ~ PC1.abs, data = arnegard)
PC1.quad <- lm(Loess_Predicted_SL_mm ~ PC1.abs + I(PC1.abs^2), data = arnegard)

PC2.lin <- lm(Loess_Predicted_SL_mm ~ PC2.abs, data = arnegard)
PC2.quad <- lm(Loess_Predicted_SL_mm ~ PC2.abs + I(PC2.abs^2), data = arnegard)

anova(PC1.lin, PC1.quad)


scatter3d(Loess_Predicted_SL_mm~Isotope_PC1+Isotope_PC2, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

scatter3d(Loess_Predicted_SL_mm~PC1.abs+PC2.abs, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

#3D version of arnegard figure
scatter3d(Loess_Predicted_SL_mm~Isotope_PC1+Isotope_PC2, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)


# Arnegard transgression figure
summarySE(arnegard.all, measurevar="Total_Long_Gill_Rakers", groupvars=c("Group")


#Mismatch Traits 3D scatter
## opening inlever and jaw protrusion
scatter3d(Loess_Predicted_SL_mm~Opening_Inlever.res+Protrusion.res, data=arnegard, 
          surface = T, fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=F)

scatter3d(Loess_Predicted_SL_mm~Buccal_Length.res+Gape.res, data=arnegard, 
          surface = T, fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=F)

# Morphometrics
## Create 3D Array
### Create df, and array of X &
### Note: Using wild fish

arnegard.wild.XY.array <- na.omit(arnegard.wild[c(1:2, 50:51, 54:89)][-c(637:713),])

coords <- as.matrix(arnegard.XY.array[-1]) ## here we say, use all columns except the first two.

is.numeric(coords)
coords3D <- arrayspecs(coords, 19, 2) ## makes the matrix a 3D array

proc <- gpagen(coords3D)
coords2d <- two.d.array(proc$coords)

consensus <- apply(proc$coords, c(1,2), mean)
consensusvec <- apply(coords2d, 2, mean)
resids <- t(t(coords2d)-consensusvec)

P <-cov(resids)
pca.stuff <- svd(P)
eigenvalues <- pca.stuff$d
eigenvectors <- pca.stuff$u

scores <- resids%*%eigenvectors


plot(scores[,1], scores[,2])

arnegard.morph <- arnegard[-c(637:713),]

arnegard.morph$morph_PC1 <- scores[,1]


cbind(arnegard[-c(637:713),]$Fifteen_Percent_Loess_Curved_Landscape_Group, )
