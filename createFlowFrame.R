library(Biobase)
library(flowCore)
library(flowViz)

# simply replace the dummy dta below with your CSV data
# dta <- matrix(rnorm(10000),ncol=10)
# dimnames(dta)[[2]] <- LETTERS[1:10]
dta <- read.csv(dat)

head(dta)

# you need to prepare some metadata
meta <- data.frame(name=dimnames(dta)[[2]],
		   desc=paste('this is column',dimnames(dta)[[2]],'from your CSV')
)
meta$range <- apply(apply(dta,2,range),2,diff)
meta$minRange <- apply(dta,2,min)
meta$maxRange <- apply(dta,2,max)

head(meta)

mat.dta <- as.matrix(dta)

# all these are required for the following steps to work

# a flowFrame is the internal representation of a FCS file
ff <- new("flowFrame",
	exprs=mat.dta,
	parameters=AnnotatedDataFrame(meta)
)

# a simple plot to check that it worked
xyplot(SSC.A~FSC.A,ff)

# now you can save it back to the filesystem
write.FCS(ff,paste(tempfile(),'FCS',sep='.'))

