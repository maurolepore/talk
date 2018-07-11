##################################################################################################################################
# Introduction to phylogenies in R: data structures, visualization, and analysis of comparative data (July 11, 2018)
# by Davey Wright 
# Peter Buck Postdoctoral Fellow in Paleobiology, NMNH (Smithsonian Institution)
# Contact info: wrightda@si.edu (or davey.f.wright@gmail.com if you happend email me in the distant future...) and/or follow me on Twitter: @Davey_F_Wright
# "I salute the echinoderms as a noble group especially designed to puzzle the zoologist." --Libbie H. Hyman (1955: p. vi)
##################################################################################################################################

# Set your working directy:

#setwd("Path_to_your_working_directory_here") # Be sure to be in the same directory as the files I sent!

# R packages we'll be using:

R_packages <- c("ape", "geiger", "phytools", "strap", "paleotree", "nlme","paleoTS")

# Uncomment the following code to install all of the R packages we'll be needing today (if you don't already have them): 

#lapply(R_packages, install.packages(R_packages), character.only = T)

# If you want, you can use the following to automatically load all of the installed packages we need. Alternatively, you can follow along below and load the packages as needed. (I did this to try to make it clear what packages are useful for different purposes.)

lapply(R_packages, library, character.only = T)

# Before we begin I should remind everyone that if you use an R package in your analyses, please be nice and cite the author's R package in your manuscript! People who write and share their code deserve credit for their hard work! Citing them is the least we can do for their contributions.


###########################################################################

# Phylogenetics in R

###########################################################################


# PRELIMINARY QUESTION: What is a "phylogeny"?

# Tree anatomy: nodes (vertices), branches (edges), tips

# What do the branch lengths represent? Is there a difference between "phylogenies", "cladograms", "phylograms", etc.?

# What, if anything, is the difference between a cladogram and a phylogeny?

# Ultrametric or non-ultrametric?

# These may sound like silly questions for learning about R packages dealing with phylogenies, but knowing in advance what kind of tree you have, the meaning of its branches (if any), and what questions you're asking of it is crucial to know before beginning any analyses. It's impossible to know where to begin unless you've thought about this.

# One thing we're not going to address is tree inference. Most analyses inferring phylogenies use software written in other languages (PAUP*, MrBayes, RevBayes, etc.). However, the R package phangorn can perform maximum parsimony and likelihood analyses. Phangorn also has a lot of useful function for dealing with phylogenies and I recommend people check out.

# There are many ways to represent tree-like structures in a programming environment (Newick, NEXUS, etc.).

# In R, most packages involving phylogenies are built using functions in the package APE. 

# APE treats evolutionary trees as objects of class "phylo"

library(ape)

# Let's build a tree using Newick Format. Many other software programs (e.g., McClade) can be used to build Newick files, which are basicaly build trees using strings of characters

# To see how it works, let's assemble a cladogram for the five extant classes of echinoderms. First, let's build the character string:

Echinodermata <-"(Crinoidea,((Asteroidea,Ophiuroidea),(Echinoidea,Holothuroidea)));"

# Now, we can use the read.tree function in APE to import the character string as a "phylo" object

E.tree <- read.tree(file = NA, text = Echinodermata)

class(E.tree)

# Now that we have a tree stored as a "phylo" object, let's plot it:

plot(E.tree)

# since E.tree is of class "phylo", R is calling the plot.phylo from APE. Check out ?plot.phylo to see arguments for changing the appearance of a plot

plot(E.tree, type = "fan")
plot(E.tree, type = "radial")
plot(E.tree, type = "cladogram", direction = "up") # Aside: Notice the correspondence between the Newick string and the tree structure?

# You can output Newick trees in R using the function write.tree:

write.tree(E.tree, file = "Echinodermata.Newick")

# This can be read in using the read.tree function as above. Alternatively, you can use the write.nexus function:

write.nexus(E.tree, file = "Echinodermata.nex")

# Many programs (PAUP*, MrBayes, BEAST, etc.) use NEXUS files. You can read them in using read.nexus:

tree <- read.nexus("Echinodermata.nex")


###########################################################################

# Inside a "phylo" object

###########################################################################

# What *is* an R object of class "phylo" anyway?

typeof(tree)

str(tree)

# Let's have a look at the components. What are they?

tree$edge
tree$tip.label
tree$Nnode

# Now let's try to make sense of the edge matrix. First let's look at the node, tip, and edge labels:

nodelabels()
tiplabels()
edgelabels()

# Aside: Sometimes you may want to write your own function that traverses the tree (Independent contrasts, Simulate BM evolution, etc.). 
# However, you'll nearly always find a useful function someone else wrote to do what you want, but how do you do this with "phylo" objects if you needed to?

# Depending on whether you're traversing from the tips to the root (or root to the tips), sometimes it's necessary to reorder the edge matrix of the tree

tree <- reorder(tree, order = "postorder") # alternative is "cladewise"

plot(tree); nodelabels(); tiplabels();

# How would you find the internal nodes? [Hint: think of the edge matrix]

InternalNodes <- unique(tree$edge[,1])

InternalNodes # For paleotrees, this assumes no taxa are related by ancestor-descendant relationships

# How to you get the branch descendants of an internal node?

d <- which(tree$edge[, 1] == InternalNodes[3]) # I picked the 3rd element here, which is actually node 7
d
# Now check the edge matrix

tree$edge 

# Cool. You can also see this is the case by plotting the edgelables:

edgelabels()

# What about descendant nodes?

tree$edge[d,2] 


# Tree manipulation and other utility functions:

Ntip(tree)

Nedge(tree)

Nnode(tree)

# A number of functions are useful for getting ancestors/descendants. For example, see the functions Ancestors, Descendants, Children, and other related functions in the R package phangorn. Also, the FindDescendants function in strap.

# one useful function for finding the node with the MRCA between two tips is mrca:

mrca(tree)

# There are many ways you can modify, prune, or bind trees:

tr <- drop.tip(tree, "Holothuroidea")   # removes a taxon or character vector of taxa
plot(tr)

Eleutherozoa <- extract.clade(tree, 7)
plot(Eleutherozoa)

# You can combine trees using the function bind.tree. Also you can re-root trees using the root function. (Many useful functions in the packages APE and paleotree (e.g, dropZLB, dropExtant, etc...))


###########################################################################

# Branch lengths and time trees

###########################################################################

# Thus far, we have only been dealing with a cladogram. The branches here have no defined length. We can assign branch lengths to
# a vector of edge.length inside the "phylo" object. For example:

tree$edge.length <- rep(1,Nedge(tree)); plot(tree)

# Now we have a vector of branch lengths
tree$edge.length

# Obviously, branches don't have to be unit length! 
tree$edge.length <-runif(Nedge(tree)); plot(tree)

# But these branch lengths are just made up. Let's simulate a birth-death tree (with many more tips than the echinoderm tree from earlier). Note that other R packages are useful for simulating birth-death and birth-death-sampling (FBD) models. See function in TreeSim, paleotree.

plot(rlineage(0.1, 0, Tmax = 30)) # Yule process (birth-death with death.rate == 0) birth.rate = 0.1

tr <- rlineage(0.1, 0.05, Tmax = 100) # constant rate birth-death process
plot(tr)

tr <- pbtree(b = .1, d = 0.05, n = 100,scale = 1) # another function to simulate birth-

# We can make trees like nice by using the function ladderize.  

plot(ladderize(tr), show.tip.label = FALSE) 

# Since we simulated a birth-death process over time, we can plot a time axis:

axisPhylo()

# # Note the tree is non-ultrametric, what if we made it ultrametric

no.fossils <- drop.fossil(tr)
plot(ladderize(no.fossils), show.tip.label=FALSE); axisPhylo()

# Compare the lineage diversity trends through time with and without fossils

ltt.plot(tr, log = "y")
ltt.plot(no.fossils, log = "y")

# Birth-death models allow one to calculate the exact probability of a tree given certain conditions (e.g., time, # of taxa, or both) 
# fit a constant rate birth-death model to the tree with no fossils (using likelihood [i.e., Pr(data | model)]:

birthdeath(no.fossils)

# What are the parameter estimates? Note the models fit parameters of Net Diversification and Turnover, which are:

#turnover extinction.rate / branching.rate
0.05 / 0.1 # == simulated value

# Net diversification (branching.rate - extinction.rate)
0.1 - 0.05 # == simulated value

# How close are the parameter estimates to their simulated values? What if you increase the extinction rate?

# More complicated processes can be simulated, such as a birth-death-sampling process with rate shifts using the R package TreeSim.

# If you're interested in this subject, be sure to check out simFossilRecord and related (and well documented) functions in the  R package paleotree. If your a paleontologist not regularly using *some* function(s) in paleotree you're missing out!


# Okay, enough (for now) with simulated trees. Let's use an empirical tree of fossil crinoids taken from a MrBayes tip-dating analysis.

# Data from : Wright, D.F. 2017. Bayesian estimation of fossil phylogenies and the evolution of early to middle Paleozoic crinoids (Echinodermata). Journal of Paleontology, 91(4): 799-814

tree <- read.nexus("Crinoid.tre")

plot(ladderize(tree)); axisPhylo()

# Cool. But look at the time scale? What's up with that?

# The units are wrong! (What are the units?)

# Correct them by dividing by the clock rate for that tree (the SI on dryad also have instructions for getting this):

clock.rate <- 0.03517385 # clock rate for the MCC tree

tree$edge.length <- tree$edge.length / clock.rate 

plot(ladderize(tree)); axisPhylo()

# What's wrong now? R doesn't know these crinoids are extinct! By default, the age of the youngest tip is zero (i.e., the present). You can fix this by setting a root.time:

tree$root.time <- 485.4

plot(ladderize(tree), cex = 0.5); axisPhylo()

# Much better! (NB: For fossil-only data like this, you have to fix the age of a tip in order to set the root.time)

# To make nice plots with geologic time scales, you can use the R package strap:

library(strap)

geoscalePhylo(ladderize(tree),direction="rightwards", units = c("Period", "Age"), boxes = "Age", arotate = 90, cex.tip = 0.6, cex.ts = 1.1, cex.age = .9, width = 2, tick.scale = "no")

# Rough example plot, but you get the idea. Play around with the arguments until it looks nice.

# This fossil phylogeny was already time-scaled using fosslil Bayesian tip-dating methods in MrBayes. Trees of extant taxa can be timescaled using either node or tip-dating methods (typically Bayesian approaches are used). For fossil-only datasets, other approaches exist, such as cal3, which is a model-based approach similar in spirit to tip-dating.
# I strongly recommend all paleo users time-scale trees using some model-based (not arbitrary) approaches, such as tip-dating or cal3 (others exist). The function bin_cal3TimePaleoPhy in the R package paleotree (written by Dave Bapst) is one of the only a posteriori methods I suggest paleontologists use at the moment (of course, this could change). Other approaches do not account for incomplete sampling and/or make indefensible assumptions (e.g., adding an arbitrary constant to minimum branch lengths). Just don't do it. 

# Seriously...don't.

# NB: All analyses done here only look at results from a single phylogeny and doesn't account for uncertainty in tree topology or divergence dates. In a real analysis, you'll want to perform analyses to summarize results over multiple trees and/or dates to accomodate these sources of uncertainty.




###########################################################################

# Trait evolution basics, visualization, and intro to comparative analysis

###########################################################################

# Brownian motion is one of the simplest models of character change. Mathematically, it's just a continuous time version of a random walk.
# e.g.,
# Trait(i) = Trait(i-1) + rnorm(n = 1, mean = 0, sd = 1)

# Brownian motion without trees:

# non-phylogenetic BM:
changes <- sapply(1:50, function(x) c(0, rnorm(999, mean = 0, sd = 0.5)))  # sd = 0.5 --> BM rate == ~0.25
plot(1:1000, cumsum(changes[1:1000, 1]), 'l', main = "Brownian motion", xlab = "Time", ylab = "Body size")

# Note that the variance increases through time. We can also see the variance increase this visually:
changes <- sapply(1:50, function(x) c(0, rnorm(999, mean = 0, sd = 0.5)))  # sd = 0.5 --> BM rate == ~0.25
plot(1:1000, cumsum(changes[1:1000, 1]), 'l', main = "Brownian motion", xlab = "Time", ylab = "Body size", ylim = c(-45,45))
for (i in 1:50){
  lines(1:1000, cumsum(changes[1:1000, i]), 'l', main = "Brownian motion", xlab = "Time", ylab = "Body size")
}

# How does do you simulate a BM process with phylogenetic autocorrelation this over a tree? I'm not going to go into detail here, but the basic idea is to have this process traverse the tree. 
# Luckily, there are already functions to simulate these kinds of processes and fit models (e.g., using AIC). Many functions in the R packages geiger and phytools are great for this kind of analysis.

# Let's use the fastBM function in phytools to simulate and visualize a character evolving under a Brownian motion process:

library(phytools)

tr <- pbtree(b = .1, d = 0.05, n = 100,scale = 1) 

#simulate data with BM
y <- fastBM(tr, internal = TRUE, sig2 = 0.5)
phenogram(tr,y,ftype="off",spread.labels=FALSE)

# Let's plot the evolution of a character across the tree (as inferred using likelihood):
y <- fastBM(tr, internal = FALSE, sig2 = 0.5)
figure <- contMap(tr, y, fsize = c(0.5,1), outline = FALSE)

# change the color map or invert it:
figure <- setMap(figure, fsize = c(0.5,1), outline = FALSE, invert = TRUE)

# with color
plot(figure, fsize = c(0.5,1), outline = FALSE)

# in black and white
figure <- setMap(figure, colors = c("black", "white"))
plot(figure, fsize = c(0.5,1), outline = FALSE)

# Before of overinterpretting your results! Look at the confidence intervals using the R package paleotree:
library(paleotree)
y <- fastBM(tr, internal = FALSE, sig2 = 0.5)
plotTraitgram(y, tr) 

# Yikes! This should worry you! 
# There's a lot that could be discussed here but no time. However, think about this anytime you're presented with a study where the interpretations rely on ancestral state reconstructions! 

# This is fun, but...
# Now let's look at some trait data for the crinoid dataset we looked at earlier do a quick comparative analysis

# reminder of the tree:
plot(ladderize(tree)); axisPhylo()

# Now let's load some trait data:

data <- read.table("Trait.data.txt")

# Two traits: (1) Body size [calyx volume], (2) Stem length [measures the height above seafloor]. 

plot(data$Body_size ~ data$Stem_length, xlab = "Stem length", ylab = "Body size")

# Is there a correlation between these traits?

cor.test(data$Body_size, data$Stem_length, method = "pearson")

# Okay, let's use Ordinary Least Squares (OLS) to test to see if there's a relationship between these two trait:

summary(lm(data$Body_size ~ data$Stem_length))

# But OLS assumes the variance in the residuals is normally distributed. Comparative data violate this because species are related (i.e., evolution)!

# Phylomorphospace plot:

phylomorphospace(tree, data, label = "off", xlab = "Body size", ylab = "Stem length")
text(data, rownames(data), pos = 2, cex = 0.5, offset = 0.5)

# Let's look at the variance-covariance matrix:

VCV <- vcv.phylo(tree)
plot(ladderize(tree)); axisPhylo()
head(VCV) # what do the elements in this matrix mean? 

tip.heights <- diag(VCV)
tip.heights # what do the diagonal elements mean?


# Phylogenetic Generalized Least Squares: relaxes assumption about normally distributed residuals, but assumes a model of trait evolution like BM.

# Phylogenetic GLS was first introduced as a way to account for spurious autocorrelation between traits arising from shared evolutionary history. Often though, the strength of the relationship is underestimated and PGLS can sometimes *improve* correlations among traits. If you have comparative data and a time-scaled tree, you really should use PGLS because otherwise you're not accounting for the non-independent nature of your observations.

# One of the easiest ways to do PGLS is to take advantage of the R package nlme. 
# Note that GLS is, of course, very general. It has many other uses in science, especially time series and spatial analysis. For example, if you had a fossil time series of body size and temerpature dataset, you could use GLS to incorporate a correlation structure for an autoregressive model. 

library(nlme)

cor.BM <- corBrownian(phy = tree) # assumes a Brownian motion model of trait evolution
data.f <- data.frame(data,row.names=row.names(data))	

result <- gls(Body_size ~ Stem_length, correlation = cor.BM, weights = varFixed(~ tip.heights), data = data.f, method = "ML") #tip.heights are used as weights because the tree is non-ultrametric! see ?gls for use
summary(result)

# What's the relationships between body size and stem length?


###########################################################################

# Beyond BM: Fitting models of trait evolution

###########################################################################

# BM isn't the only game in town. Other models are also available. Most models are BM with added complexity, so they simplify to BM:

tr <- pbtree(b = .2, d = 0.05, n = 100,scale = 1) # another function to simulate birth-death trees. This one is in Phytools.

# Simulate data with a trend (i.e., the mean of the rnom() draw is != zero. Requires a non-ultrametric tree)
x <- fastBM(tr, internal = TRUE, mu = 3)
phenogram(tr, x , ftype="off", spread.labels=FALSE)

# Simulate data with Ornstein-Uhlenbeck process with a single, stationary peak. It's like BM with a "rubber band" parameter alpha that "pulls" traits to "optimal" value.
# In a macroevolutionary context, it's used to model the evolution of a clade within a Simpsonian adaptive zone.
z <- fastBM(tr, internal = TRUE, a = 10, sig2 = 1, mu = 0)
phenogram(tr, z, ftype = "off", spread.labels = FALSE)

library(geiger) # needed for rescale function, we'll discuss more about geiger later...

# Simulate data with Eary Burst diversification. It's basically BM with an exponentially declining rate paramter. Clades undergoing adaptive radiation are predicted to fit this pattern.
z <- fastBM(rescale(tr, "EB", a = -2), sig2 = 1)
phenogram(tr, z, ftype = "off", spread.labels = FALSE)

# Simulate data with Late Burst. Models a burst of trait diversification occurring sometime *after* the early evolution of a clade. It's like an EB process but with the opposite sign. 
z <- fastBM(rescale(tr, "EB", a = 2), sig2 = 1)
phenogram(tr, z, ftype = "off", spread.labels = FALSE)


# Let's revisit the crinoid body size data:

size <- data$Body_size
names(size) <- rownames(data)

# Is BM the best model of evolution for these data?

# We can use the fitContinuous function in the R package geiger to estimate likelihoods and compare AICs across alternative models.
# There's a lot of "philosophy" about model selection to discuss here, but I don't have enough time to get to it...

BM <- fitContinuous(phy = tree, dat = size, model = "BM")
BM # let's see parameter estimates, likelihood, AIC, etc.

# Let's compare Brownian motion, Orstein-Uhlenbeck (single peak), Early Burst, and White Noise.
# White noise assumes no covariance among species, which is sort of equivalent to saying it's a non-phylogenetic model. 
# Another way of saying it is that under white noise, species traits evolve so fast they preserve no record of evolutionary history in their traits.

# For the moment, let's worry only about the "best" model fit using AICc (bad practice, look at parameter estimates!).

# So let's write a quick script to get AICc for several models:

models <- c("BM", "OU", "EB", "white") 
AICc <- vector(mode = "numeric", length = length(models)); names(AICc) <- models

for (i in 1:length(models)){
  temp <- suppressWarnings(fitContinuous(phy = tree, dat = size, model = models[i])) # the suppressWarnings is used to get rid of a warning about how fitContinuous is using the VCV to compute likelihoods for the OU model with a non-ultrametric tree
  AICc[i] <- temp$opt$aicc
}

# Now let's plot the Akaike weights for each model. These are loosely interpreted (by some) as model probabilites, which tacitly assumes the "True" model is one of the models you compared. Not realistic, but you have to deal with these concerns one way or another (need more time for philosophy here!).

# transform AICc to Akaike weights using the akaike.wts function in paleoTS:
library(paleoTS)
weights <- akaike.wts(AICc)

# make a barplot of results
barplot(weights, ylim = c(0,1), main = "Akaike Weights")

# Which model best fits body size evolution? How good is "best"? Good enough?

# Which model has the worst fit? Why?


# Other packages with addional models of trait evolution include OUwei, ouch, mvMORPH, surface


################################################

# Discrete models of trait evolution

################################################


# You can also simulate character data for tips using the sim.char function in geiger. It can be used to simulate either continuous
# time BM or discrete traits

tr <- pbtree(b = .1, d = 0.05, n = 100,scale = 1) # get another tree because why not?

# Let's simulate a discrete character:
x <- sim.char(tr, par = list(rbind(c(-1, 1), c(1, -1))), model = "discrete") # discrete characters (1's and 2's). Needs a Q (rate) matrix. See ?sim.char

# We can fit models of discrete character evolution using geiger just like we did for the body size data. Let's fit and compare an
# "Equal Rates" (ER) model with a "All Rates Different" (ARD) model:

# Fit under equal rates model
ER <- fitDiscrete(tr, dat = x[,,1], model = "ER")
ER
# Get the sample size corrected AIC:
ER.aicc <- ER[[4]]$aicc
ER.aicc

# Git under "all rates different" model"
ARD <- fitDiscrete(tr, dat = x[,,1], model = "ARD")

ARD.aicc <- ARD[[4]]$aicc
ARD.aicc

# Which model is the better fit? What one should we expect to be a better fit?

# Let's get ancestral states of this character and plot a pie chart of the likelihood of a given state 

anc.char <- ace(x, tr, type = "discrete")
plot(ladderize(tr), show.tip.label = FALSE)
nodelabels(pie = anc.char$lik.anc, piecol = c("lightblue","salmon"), cex = 0.5)

# The R package phangorn is also a great tool for this kind of analysis (e.g., it can get parsimony-based estimates of ancestral states, etc.).


# Lastly, we can use the R package phytools to conduct stochastic character mapping analysis. It's an alternative approach to maximum parsimony or maximum likelihood techniques for ASR that can account for uncertainty by summarizing simulations over  many possible character evolution histories.

# get data in correct format
x <- as.vector(x)
names(x) <- tr$tip.label

# simulate character evolution:
char.histories <- make.simmap(tr, x, model = "ER", nsim = 100)

# plot

plotSimmap(char.histories, pts = TRUE) 

# To summarize the total distribution of simulated character histories, use the densityMap function:

densityMap(char.histories, fsize = c(0.5,1))

# Cool! Right?

# Another package very useful for simulating and analyzing discrete (morphological) characters is dispaRity. Again, paleotree has a number of functions that might be useful in this context.


##########################################################################################################################
# NB:
# The crinoid body size and stem length data are NOT real and have nothing to do with the empirical values for those taxa.
# They were simulated only for use in this workshop and biologically nonsensical. 
##########################################################################################################################