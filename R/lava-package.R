
##' Estimation and simulation of latent variable models
##' 
##' Framwork for estimating parameters and simulate data from Latent Variable
##' Models.
##' 
##' @name lava-package
##' @import utils stats guts
##' @aliases lava-package lava
##' @docType package
##' @author Klaus K. Holst Maintainer: <k.k.holst@@biostat.ku.dk>
##' @keywords package
##' @examples
##' 
##' lava()
##' 
NULL

##' Longitudinal Bone Mineral Density Data
##' 
##' Bone Mineral Density Data consisting of 112 girls randomized to receive
##' calcium og placebo. Longitudinal measurements of bone mineral density
##' (g/cm^2) measured approximately every 6th month in 3 years.
##' 
##' 
##' @name calcium
##' @docType data
##' @format A data.frame containing 560 (incomplete) observations. The 'person'
##' column defines the individual girls of the study with measurements at
##' visiting times 'visit', and age in years 'age' at the time of visit. The
##' bone mineral density variable is 'bmd' (g/cm^2).
##' @source Vonesh & Chinchilli (1997), Table 5.4.1 on page 228.
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name bmd
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name brisa
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name bmidata 
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name hubble
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name hubble2
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name hubble2
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name indoorenv
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name missingdata
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name nldata
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name nsem
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name semdata
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Serotonin data
##'
##' \tabular{rll}{
##'         day    \tab numeric \tab Scan day of the year \cr
##'         age    \tab numeric \tab Age at baseline scan \cr
##'         mem    \tab numeric \tab Memory performance score \cr
##'         depr   \tab numeric \tab Depression (mild) status 500 days after baseline \cr
##'         gene1  \tab numeric \tab Gene marker 1 (HTR2A) \cr
##'         gene2  \tab numeric \tab Gene marker 2 (HTTTLPR) \cr
##'         cau \tab numeric \tab SERT binding, Caudate Nucleus \cr
##'         th  \tab numeric \tab SERT binding, Thalamus \cr
##'         put \tab numeric \tab SERT binding, Putamen \cr
##'         mid \tab numeric \tab SERT binding, Midbrain \cr
##'         aci \tab numeric \tab 5-HT2A binding, Anterior cingulate gyrus \cr
##'         pci  \tab numeric \tab 5-HT2A binding, Posterior cingulate gyrus \cr
##'         sfc \tab numeric \tab 5-HT2A binding, Superior frontal cortex \cr
##'         par \tab numeric \tab 5-HT2A binding, Parietal cortex \cr
##' }
##' @name serotonin
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name serotonin2
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL

##' Data
##'
##' Description
##' @name twindata
##' @docType data
##' @format data.frame containing ...
##' @keywords datasets
NULL


##' For internal use
##' 
##' @title For internal use
##' @name startvalues
##' @rdname internal
##' @author Klaus K. Holst
##' @keywords utilities
##' @export
##' @aliases addattr startvalues2 startvalues3 starter.multigroup
##' modelPar modelVar matrices pars pars.lvm pars.lvmfit pars.glm
##' score.glm procdata.lvmfit modelPar modelVar matrices reorderdata
##' graph2lvm igraph.lvm subgraph finalize index.lvm index.lvmfit
##' index reindex index<- survival survival<- randomslope
##' randomslope<- lisrel variances offdiags describecoef parlabels
##' stdcoef CoefMat CoefMat.multigroupfit deriv updatelvm
##' checkmultigroup profci estimate.MAR missingModel Inverse
##' gaussian_logLik.lvm addhook gethook multigroup Weight fixsome
##' parfix parfix<- merge IV parameter
##' Specials decomp.specials getoutcome
NULL
