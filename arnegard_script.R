# Arnegard paper analysis

library(car)
library(dplyr)
library(geomorph)
library(ggplot2)
library(MASS)
library(MuMIn)
library(psych)
library(rgl)
library(scatterplot3d) 

arnegard.all <- read.csv('data-raw/arnegard_data.csv')

### Drop the non-experiment fish
arnegard <- arnegard.all[-c(637:713),]
arnegard.noO <- subset(arnegard, Fifteen_Percent_Loess_Curved_Landscape_Group == (c("B", "L", "A")))
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

#Run a PCA
arnegard.morph.terms <- c("Std_Length_TPS1_to_TPS14_mm", "Opening_Inlever.res", "Protrusion.res", "Buccal_Length.res", "Total_Long_Gill_Rakers", "Total_Short_Gill_Rakers", "Total_All_Gill_Rakers", "Gape.res", "Suction_Index", "EpaxH.res")
arnegard.matrix <- arnegard[arnegard.morph.terms]
arnegard.pca <- principal(arnegard.matrix, nfactors = 2, rotate = "varimax")
plot(arnegard.pca$scores)


# Generate DFA to dilineate B & L groups
#First create dataset with just B and L

arnegard.BL <- arnegard %>%
  group_by(Fifteen_Percent_Loess_Curved_Landscape_Group) %>% 
  filter(Fifteen_Percent_Loess_Curved_Landscape_Group == "L" | Fifteen_Percent_Loess_Curved_Landscape_Group == "B")

#Create dataset with B, L and then 'other'
arnegard.3G <- arnegard %>%
  group_by(Fifteen_Percent_Loess_Curved_Landscape_Group) %>%
  mutate(Loess.3G = revalue(Fifteen_Percent_Loess_Curved_Landscape_Group, c("A" = "O")))

#Run DFA
z.noO <- lda(Fifteen_Percent_Loess_Curved_Landscape_Group ~ Std_Length_TPS1_to_TPS14_mm + Opening_Inlever.res + Protrusion.res + Buccal_Length.res + Total_Long_Gill_Rakers + Total_Short_Gill_Rakers + Total_All_Gill_Rakers + Gape.res + Suction_Index + EpaxH.res, data = arnegard.noO)

z.3G <- lda(Loess.3G ~ Std_Length_TPS1_to_TPS14_mm + Opening_Inlever.res + Protrusion.res + Buccal_Length.res + Total_All_Gill_Rakers + Gape.res + Suction_Index + EpaxH.res, data = arnegard.3G)

#Can do two or three groups
plot(z.3G)
print(z)
z$scaling

LD.scores <- predict(z.3G)$x
z.3G$class

#Make dataset with groups and plot
z1 <- predict(z.3G, arnegard.3G)
z1 <- as.data.frame(z1)

#Bring in 'fitness' column
z1$fitness <- arnegard.3G$Loess_Predicted_SL_mm

#Also group column
z1$group <- arnegard.3G$Loess.3G

#Plot the LDs and circle...

ggplot(z1, aes(x.LD1, x.LD2, color = group)) +
  geom_point() +
  stat_ellipse()

#See if the groups fall out on LDs
TukeyHSD(aov(x.LD1 ~ group, data = z1))

#Look at scatterplot
scatter3d(fitness~x.LD1+x.LD2, data=z1, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

#determine parameters
summary(lm(fitness ~ abs(x.LD2), data = z1))
summary(lm(fitness ~ abs(x.LD1), data = z1))


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

### Cool fig with abs pcs
scatter3d(Loess_Predicted_SL_mm~PC1.abs+PC2.abs, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

### same as above but with real length
scatter3d(Std_Length_TPS1_to_TPS14_mm~PC1.abs+PC2.abs, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

scatter3d(Std_Length_TPS1_to_TPS14_mm~Opening_Inlever.res+Protrusion.res, data=arnegard, 
          fit="quadratic", residuals=TRUE, bg="white", axis.scales=TRUE, grid=TRUE, 
          ellipsoid=FALSE)

summary(aov(Std_Length_TPS1_to_TPS14_mm~Isotope_PC1*Isotope_PC2, data=arnegard))

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

arnegard.wild.XY.array <- na.omit(arnegard.wild[c(1:2, 50:51, 54:89)][-c(637:713),], stringsAsFactors = F)

arnegard.wild.XY.array <- arnegard.wild.XY.array %>% 
  rename(ID = Individ_Num)

arnegard.wild.coords <- as.matrix(arnegard.wild.XY.array[,-(1:2)]) # make it numeric
is.numeric(coords) #check


arnegard.wild.coords.3D <- arrayspecs(arnegard.wild.coords, 19, 2) ## makes the matrix a 3D array; # p = 19 is num measures and k=2 num of dimensions

classifiers <- arnegard.wild.XY.array[,1:2]
is.factor(classifiers$Group) #Verify that group is a factor

# Can split by factor... not sure if this is what I want
arnegard.wild.benthic <- coords.subset(A = arnegard.wild.coords.3D, group = classifiers$Group)

## Generalized Procrustes Analysis 

arnegard.Proc <- gpagen(arnegard.wild.coords.3D)
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
