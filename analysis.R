# Multivariate Analysis - PCA and EFA
# Holzinger-Swineford (1939) dataset
# Author: Alexis Alberto Zuniga Alonso


# 1. Libraries -----------------------------------------------------------------

libs <- c("lavaan", "psych", "GPArotation", "FactoMineR",
  "factoextra", "corrplot", "skimr", "tidyverse")

# Install any missing packages, then load them all
to_install <- setdiff(libs, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install)

invisible(lapply(libs, library, character.only = TRUE))


# 2. Data (a) ------------------------------------------------------------------

# The dataset ships with lavaan
data <- HolzingerSwineford1939

# Structure and descriptive statistics of all variables
glimpse(data)
skim(data)

# Keep only the nine test scores for the analysis
df <- data |> select(x1:x9)

corrplot(cor(df), method = "number")


# 3. Principal Component Analysis (b) ------------------------------------------

# Suitability: KMO > 0.6 and a significant Bartlett test
KMO(df)
cortest.bartlett(cor(df), n = nrow(df))

# PCA on standardised variables
pca <- PCA(df, scale.unit = TRUE, graph = FALSE)


# 4. Number of components (c) --------------------------------------------------

eig <- get_eigenvalue(pca)
eig

# Kaiser: keep eigenvalues above 1
print(
  fviz_eig(pca, addlabels = TRUE, choice = "eigenvalue") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red")
)

# 60% rule: keep the fewest components reaching 60% cumulative variance
print(
  ggplot(eig, aes(x = seq_len(nrow(eig)), y = cumulative.variance.percent)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 60, linetype = "dashed", colour = "red") +
    scale_x_continuous(breaks = seq_len(nrow(eig))) +
    labs(x = "Components", y = "Cumulative variance (%)")
)

# Horn's parallel analysis: keep eigenvalues above those of random data
set.seed(123)
fa.parallel(df, fa = "pc", n.iter = 500)


# 5. Interpreting the components -----------------------------------------------

# Correlations between variables and the first three components
pca$var$cor[, 1:3] |> round(2)

# Quality of representation (cos2) and contributions per component
corrplot(pca$var$cos2[, 1:3], is.corr = FALSE)
print(fviz_contrib(pca, choice = "var", axes = 1))
print(fviz_contrib(pca, choice = "var", axes = 2))
print(fviz_contrib(pca, choice = "var", axes = 3))

# Correlation circles and biplot coloured by school
print(fviz_pca_var(pca, axes = c(1, 2), col.var = "contrib", repel = TRUE,
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
print(fviz_pca_var(pca, axes = c(1, 3), col.var = "contrib", repel = TRUE,
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")))
print(fviz_pca_biplot(pca, repel = TRUE, label = "var", col.ind = "grey70",
  habillage = data$school, addEllipses = TRUE))


# 6. Selection for the EFA (d) -------------------------------------------------

# Parallel analysis on the factor model (common variance only)
set.seed(123)
fa.parallel(df, fm = "ml", fa = "fa", n.iter = 500)

# Individual KMO per variable (should be above 0.6)
KMO(df)$MSAi |> round(2)


# 7. Exploratory Factor Analysis (e) -------------------------------------------

# 3 factors, maximum likelihood, oblique rotation (factors may correlate)
efa <- fa(df, nfactors = 3, fm = "ml", rotate = "oblimin")

# Loadings (hiding those below 0.3), communalities and factor correlations
print(efa$loadings, cutoff = 0.3, digits = 2)
efa$communality |> round(2)
efa$Phi |> round(2)

# Fit indices: TLI > 0.95 and RMSEA < 0.06 indicate good fit
efa$TLI
efa$RMSEA

# Diagram of which variables load on each factor
fa.diagram(efa)
