# Redundancy Analysis */
# 2012-10-20 CJS First edition */

# This is taken from Legendre and Legendre, Numerical Ecology, Second Edition, Table 11.3 */

# Lines starting in ##--part001b; and ##---part001e; are used
# to bracket portions of the R code for inclusion into my
# course notes and are not normally coded.

options(useFancyQuotes=FALSE) # renders summary output correctly


# Get the raw data
fish <- read.csv("fish.csv", header=TRUE)
fish

# Load the appropriate package
library(vegan)



# Do the redundancy analysis
# It is not necessary to specify the 'other class' as it is an indicator
# variable that is not coral or sand.

# Extract the data into separate matrices
Species <- as.matrix(fish[,paste('sp',1:6, sep="")])
Species
# We will center the species with mean 0 but no change to the std deviation
# We don't need to standardize to std=1 because measures are comensurate
Species <- Species - matrix(apply(Species,2,"mean"), nrow=nrow(Species),
                         ncol=ncol(Species),byrow=TRUE)
cat("\n\nStandardized Species scores\n")
Species

Env <- as.matrix(fish[,c("depth", "coral", "sand")])
Env
# Standardize to mean 0 and variance 1
Env <- Env - matrix(apply(Env,2,"mean"), nrow=nrow(Env),
                         ncol=ncol(Env),byrow=TRUE)
Env <- Env / matrix(apply(Env,2,"sd"), nrow=nrow(Env),
                         ncol=ncol(Env),byrow=TRUE)
cat("\n\nStandardized Environmental Covariates\n")
Env


result <- rda ( Species ~ Env)

# Get the results
summary(result)

# Extract the species and site scores on the ordination axes 
#    choice=c(1,2) indicates the first two axes
#    scaling=1 (site), 2 (species), or 3 (both) are scales by eigenvalues
# Refer to vignette("decision-vegan") for a discussion of the choice of
# scaling.
# This gives the scores in Table 11.4 of Legendre and Legendre, Numerical Ecology
const.sites <- sqrt((nrow(Species)-1)*sum(eigenvals(result)))
scores(result, choice=c(1,2), scaling=1, const=c(1,const.sites))



# Do a biplot
# The Biplot function doesn't work because RDA is a constrained ordination
biplot(result) # fails

# .. so we construct our own plot
# See  https://stat.ethz.ch/pipermail/r-sig-ecology/2011-May/002125.html
plot(result, dis=c("sp","cn"))
plot(envfit(result, Species, display="lc"), add = TRUE, col="red")




# Use the rdaTest library from Legendre directly.
# Available from https://stat.ethz.ch/pipermail/r-sig-ecology/2011-May/002125.html
# Get the tar ball 
#     http://www.bio.umontreal.ca/legendre/indexEn.html#RFunctions
# Unpack the tar ball (likely need to do this outside of Windows)
# Start the CMD option of windows and navigate to the directory containing rdaTest
#     use R CMD INSTALL rdaTest
# to install the package.
library(rdaTest)

result2 <- rdaTest(Species, Env, test.F=FALSE, nperm=999)
summary(result2)

# List the various scores

# Species scores, scaling=1 and 2
result2$U.sc1
result2$U.sc2

# Site scores, scaling=1 and 2
result2$F.sc1
result2$F.sc2

# Fitted objects, scaline=1 and 2
result2$Z.sc1
result2$Z.sc2

# Biplot scores, scaling=1 and 2 (#2 is actual correlation with Canonical variables)
result2$biplotScores1
result2$biplotScores2

# distance biplot	
plot(result2, graph.type="Z",scaling=1)

# correlation biplot
plot(result2, graph.type="Z",scaling=2)