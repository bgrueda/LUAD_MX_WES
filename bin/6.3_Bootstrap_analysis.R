library(ggplot2)
library(tidyr)
library(dplyr)

# Pancancer data
P <- read.csv("P_18862.csv", header = F)
# Cohort data
C <- read.csv("C_18862.csv", header = F)

# Select genes per decile 
n <- 18862
r <- 10
decilP <- vector("list", length= r)
for (i in 1:r) {
  init <- (i - 1) * (n/r) + 1
  end <- i * (n/r)
  decilP[[i]] <- P$V1[init:end]
  }

decilC <- vector("list", length= r)
for (i in 1:r) {
  init <- (i - 1) * (n/r) + 1
  end <- i * (n/r)
  decilC[[i]] <- C$V1[init:end]
}

# Look for real intersections
intersections <- vector("list", length= r)
intersect_len <- vector("list", length= r)
for (i in 1:r) {
  Int <- intersect(decilC[[i]], decilP[[i]])
  len <- length(Int)
  intersections[[i]] <- c(intersections[[i]], Int)
  intersect_len[[i]] <- c(intersect_len[[i]], len)
}


# Iterations
iterations <- vector("list", length= 1)

for (i in 1:50000) {
  A <- sample(P$V1, size=(n/r))
  B <- sample(C$V1, size=(n/r))
  I <- intersect(A,B)
  value <- length(I)
  iterations[[1]] <- c(iterations[[1]], value)
}


# Distribution

data <- as.data.frame(iterations[[1]])
colnames(data) <- "Intersect"
int_plot <- unlist(intersect_len) 
dist <- ggplot(data, aes(x=Intersect)) +
  geom_histogram(stat="count", fill = "#9FA9FF", alpha = 0.7) +
  geom_vline(xintercept = 233,  color = "tomato2", size= 0.5) +
  theme_minimal()

dist

#z-scores calculation
media <- mean(iterations[[1]])
stdev <- sd(iterations[[1]])
z_scores <- vector("list", length= 1)
for (i in 1:r) {
  zs <- (intersect_len[[i]] - media) / stdev
  z_scores[[1]] <- c(z_scores[[1]], zs)
}

# p-value calculation 
pvals <- vector("list", length= r)
for (i in 1:r) {
  if (intersect_len[[i]] > media) {
    val <- sum(data > intersect_len[[i]])
  } else {
    val <- sum(data < intersect_len[[i]])  
  }
  p_val <- (val/50000)
  pvals[[i]] <- c(pvals[[i]], p_val)
}

#---------------------------------------------------------------------------------------------

# Correlations (pairs) 
corr <- read.table("join_C_vs_P_ONCO.txt", sep= " ")

plot_corr <- ggplot(corr, aes(x = V2, y = V3)) +
  geom_point(col = "#28B993", size= 0.5) +
  labs(title = "Correlation",
       x = "pglobal_cv (Pancancer oncodrive)",
       y = "pglobal_cv (Pancancer report)")
plot_corr  

p_Val_pearson <- cor.test(corr$V2, corr$V3, method = "pearson")
p_Val_pearson

p_Val_spearman <- cor.test(corr$V2, corr$V3, method = "spearman")
p_Val_spearman


# Calculation of the mean p-value in intervals of 1000
corr_sort<- corr[order(corr$V2), ]
prom <- list()
for (i in seq(1, nrow(corr_sort), by = 1000)) {
  section <- corr_sort[i:(i + 999), ]
  promA <- mean(section$V2)
  promB <- mean(section$V3)
  prom <- c(prom, list(data.frame(V1 = "Mean", V2 = promA, V3 = promB)))
}
df_prom <- do.call(rbind, prom)

# Boxplot
grp <- rep(1:ceiling(nrow(corr_sort) / 1000), each = 1000, length.out = nrow(corr_sort))
df_grp <- cbind(corr_sort, grupo = grp)
colnames(df_grp) <- c("gen", "cohort", "pancancer", "group")

df_grp2 <- pivot_longer(data= df_grp, cols='cohort':'pancancer', 
                        names_to="tool", 
                        values_to= "p_vals")

ggplot(df_grp2, aes(x=as.factor(grupo), y=p_vals, fill=tool)) + 
  geom_boxplot()


data_plot <- layer_data(plot)
medians_grp <- aggregate(cbind(cohorte, pancancer) ~ grupo, data = df_grp, FUN = median)
cor.test(medians_grp$cohorte, medians_grp$pancancer, method = "pearson")
cor.test(medians_grp$cohorte, medians_grp$pancancer, method = "spearman")


# Test intersections of first decile
DN <- read.table("D_P.txt", header = F)
PC <- read.table("P_P.txt", header = F)
MT <- read.table("M_P.txt", header = F)

int_DM <- intersect(DN, MT)
len_DM <- nrow(int_DM)

int_DMP <- intersect(int_DM, PC)
len_DMP <- nrow(int_DMP)

int_nopan <- setdiff(int_DM, PC)
len_nopan <- nrow(int_nopan)

cohort_proportion <- len_DMP/len_nopan

#Distribution 
it_pan <- vector("list", length= 1)

for (i in 1:50000) {
  D <- sample(ALL$V1, size=1700)
  E <- sample(ALL$V1, size=1700)
  int_DE <- intersect(D,E)
  int_DEP <- intersect(int_DE, PC)
  len_DE <- length(int_DE)
  len_DEP <- length(int_DEP)
  NP <- len_DE - len_DEP
  prop_dist <- len_DEG/NP
  it_pan[[1]] <- c(it_pan[[1]], len_DEP)
}

dataP <- as.data.frame(it_pan[[1]])
colnames(dataP) <- "Intersect"
ggplot(dataP, aes(x=Intersect)) +
  geom_histogram(stat="count", fill = "navy", alpha = 0.7) +
  geom_vline(xintercept = 0.1273, color= "tomato2", size= 0.5 ) 

mediaP <- mean(it_pan[[1]])
stdevP <- sd(it_pan[[1]])
zP <- (0.1273 - mediaP) / stdevP
Pval_P <- (sum(dataP < 0.1273))/50000

# Correlations per sliding windows
correlation <- read.table("/Users/bethz/Downloads/H2/correlation/join_onco_vs_mutsig.txt", sep= " ")
corrs<- correlation[order(correlation$V3), ]
n1 <- nrow(corrs)
l <- as.integer(n1/1000)

I <- vector("list", length= l)
for (i in 1:l) {
  init <- ((i - 1) * 1000) + 1
  end <- i * 1000
  I[[i]] <- corrs$V2[init:end]
}

J <- vector("list", length= n1/1000)
for (i in 1:l) {
  init <- (i - 1) * 1000 + 1
  end <- i * 1000
  J[[i]] <- corrs$V3[init:end]
}

corr_pvalues <- vector("list", length= l)
corr_Rvalues <- vector("list", length= l)
for (i in 1:l) {
  resultado_corr <- cor.test(I[[i]],  J[[i]], method = "spearman")
  PV <- resultado_corr$p.value
  RV <- resultado_corr$estimate
  corr_pvalues[[i]] <- c(corr_pvalues[[i]], PV)
  corr_Rvalues[[i]] <- c(corr_Rvalues[[i]], RV)
}

corr_df <- data.frame(no = c(1:l), rho = unlist(corr_Rvalues), pvals = unlist(corr_pvalues))

Rplot <- ggplot(data = corr_df, aes(x = no, y = rho)) +
  geom_bar(stat = "identity", fill = "#4660C1") +
  scale_x_continuous(limits = c(1, l)) +  
  theme(axis.title.x = element_blank()) 

PV_plot <- ggplot(data = corr_df, aes(x = no, y = pvals)) +
  geom_point(color = "#45C18D") +
  geom_line(color = "#45C18D") +
  theme(axis.title.x = element_blank()) +
  scale_x_continuous(limits = c(1, l)) + 
  geom_hline(yintercept = 0.05, color= "tomato2")


est <- plot_grid(Rplot, PV_plot, nrow = 2)
est


#### PROPORTION DISTRIBUTIONS 
proportions <- numeric(70)
DMP_data <- numeric(70)
No_pan <- numeric(70)
PR <- read.table("PR.txt", header = F, nrows = 1700)

for (i in 1:70) {
  tabD <- paste0(i,"_DR.txt")
  tabM <- paste0(i,"_MR.txt")
  DR <- read.table(tabD, header = F, nrows = 1700)
  MR <- read.table(tabM, header = F, nrows = 1700)
  int_DM <- intersect(DR$V1, MR$V1)
  len_DM <- length(int_DM)
  int_DMP <- intersect(int_DM, PR$V1)
  len_DMP <- length(int_DMP)
  int_nopan <- setdiff(int_DM, PR$V1)
  len_nopan <- length(int_nopan)
  proportion_d <- len_nopan/(len_nopan+len_DMP)
  DMP_data[i] <- len_DMP
  No_pan[i] <- len_nopan
  proportions[i] <- proportion_d
}

df_proportions <- data.frame(valor = proportions)
ggplot(df_proportions, aes(x = valor)) +
  geom_histogram(binwidth = 0.01, fill = "#4660C1") + 
  geom_vline(xintercept = c(0.424, 0.887), color= "tomato2", size= 0.5 ) +
  theme_minimal()


### New distributions 
CDrivers <- read.table(file= "/Users/bethz/Downloads/Venn/luad_drivers.txt")
P1485 <- read.table(file = "/Users/bethz/Downloads/Venn/pancancer_3tools/pan_common.txt")
C1107 <- read.table(file = "/Users/bethz/Downloads/Venn/SCommon.txt")

dist_drivers <- numeric(50000)
T1 <- C1107 

# For cycle
for (i in 1:50000) {
  sample <- sample(T1$V1, size = round(0.1 * nrow(T1)))
  int_d <- intersect(sample, CDrivers$V1)
  dist_drivers[i] <- length(int_d)
}

mean(dist_drivers)
sd(dist_drivers)
(8 - 5.69042) / 2.247158
(sum(dist_drivers > 8))/50000
df_driv <- data.frame(valor = dist_drivers)
ggplot(df_driv, aes(x = valor)) +
  geom_histogram(binwidth = 1, fill = "#4660C1") + 
  geom_vline(xintercept = c(0,4,5,8), color= "tomato2", size= 0.5 ) 
