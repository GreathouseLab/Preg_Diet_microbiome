# =================================================== #
# =================================================== #
#              			                                  
#            Microbiome Analysis Functions
#                    Jun Chen, PhD
#                                         
# =================================================== #
# =================================================== #
# Created: 06/10/2018                                 #
# Last Editted: 06/10/2018                            #
# By: R. Noah Padgett                                 #
# =================================================== #
# =================================================== #
# Copyright R. Noah Padgett, 2018
#
# This script is not guaranteed to be free of bugs and/or errors.
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
#
#
# You can edit script files (for example, this file)
# and either cut and paste lines from file into R command line
# (in Windows you can use ctrl-R to do this)
# or in the R command line type:
#
#    source("main.R")
#
# You may need to use full path name in the filename, or alternatively in the R console
# window change to the directory containing the file by using the command:
#
#    setwd("<path of your directory>") 
#
# Or this file can be sourced through a different file that needs
# the functions that are listed in this file. 
# =================================================== #
# =================================================== #
# This file contains the following functions
# 
#     Function           
#   subset_data()  
#   is.na.null() 
#   perform_differential_analysis_para()
#   perform_differential_analysis_para_single_FE()
#       - FE = Fixed Effects
#   
# =================================================== #
# =================================================== #


# =================================================== #
#                   subset_data()
# =================================================== #
# Inputs
#   data.obj = a dataframe
#   samIDs = a vector of IDs to subset to
# Returns subsetted dataset
subset_data <- function (data.obj, samIDs) {
    
    data.obj$meta.dat <- data.obj$meta.dat[samIDs, , drop=FALSE]
    
    if (!is.na.null(data.obj$otu.tab)) {
        data.obj$otu.tab <- data.obj$otu.tab[, samIDs, drop=FALSE]
        data.obj$otu.tab <- data.obj$otu.tab[rowSums(data.obj$otu.tab) != 0, , drop=FALSE]
        data.obj$otu.name <- data.obj$otu.name[rownames(data.obj$otu.tab), , drop=FALSE]
        if (!is.na.null(data.obj$otu.name.full)) {
            data.obj$otu.name.full <- data.obj$otu.name.full[rownames(data.obj$otu.tab), , drop=FALSE]
        }
    }
    
    if (!is.na.null(data.obj$abund.list)) {
        data.obj$abund.list <- lapply(data.obj$abund.list, function(x) {
            xx <- x[, samIDs, drop=FALSE]
            xx <- xx[rowSums(xx) != 0, , drop=FALSE]
        })
    }
    
    if (!is.na.null(data.obj$size.factor)) {
        data.obj$size.factor <- data.obj$size.factor[samIDs]
    }
    
    if (!is.na.null(data.obj$ko.list)) {
        data.obj$ko.list <- lapply(data.obj$ko.list, function(x) {
            xx <- x[, samIDs, drop=FALSE]
            xx <- xx[rowSums(xx) != 0, , drop=FALSE]
        })
    }
    if (!is.na.null(data.obj$cog.list)) {
        data.obj$cog.list <- lapply(data.obj$cog.list, function(x) {
            xx <- x[, samIDs, drop=FALSE]
            xx <- xx[rowSums(xx) != 0, , drop=FALSE]
        })
    }
    data.obj
}

# =================================================== #
#                   is.na.null()
# =================================================== #
# This function is called by other functions as a logical check
# Do not mess with this function!
is.na.null <- function (x) {
    if (is.null(x)) {
        return(TRUE)
    } else {
        if (is.na(x)[1]) {
            return(TRUE)
        }  else {
            return(FALSE)
        }
    }
    
}


# =================================================== #
#           perform_differential_analysis_para()
# =================================================== #
# The following function performs the taxonomic differential analysis
# this is a monster of a function that does a LOT with lots of options.
# 
# An example call:
#   perform_differential_analysis_para(data.obj0,  grp.name='Status.cat', 
#               adj.name=c('Gender', 'Age'),  RE=FALSE, method='NB', 
#               taxa.levels=c('Genus', 'Species'), winsor=TRUE, 
#               winsor.qt=0.97, norm='TSS', norm.level='Species', 
#               intersect.no=4, prev=0.1, minp=0.002, medianp=NULL, 
#               mt.method='raw', cutoff=0.05, ann=paste0(df, '.BMI.TSSNB'))
#
# Arguments/Inputs
#   data.obj            Data for analysis
#   grp.name            
#   adj.name
#   subject
#   RE                  Random Effects? - Logical
#   method 
#   zerop.cutoff
#   ZINB
#   LRT
#   taxa.levels
#   winsor
#   winsor.qt
#   norm
#   norm.level
#   intersect.no
#   prev
#   minp
#   medianp
#   mt.method
#   cutoff
#   ann
#   ...
perform_differential_analysis_para <- function (data.obj,  grp.name, adj.name=NULL, 
                                                subject=NULL, RE=FALSE, method='Adaptive0', 
                                                zerop.cutoff=0.25, ZINB='ZINB1', LRT=FALSE, 
                                                taxa.levels=c('Phylum', 'Order', 'Class', 
                                                              'Family', 'Genus'), 
                                                winsor=TRUE, winsor.qt=0.97, norm='GMPR', 
                                                norm.level='Genus', intersect.no=4,prev=0.1, 
                                                minp=0.002, medianp=NULL, mt.method='fdr', 
                                                cutoff=0.15, ann='', ...) {
    # To be completed
    # subject holds the random effects formula
    if (!RE) {
        if (!(method %in% c('ZINB', 'B', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model is not supported!\n')
        perform_differential_analysis_para_single <- perform_differential_analysis_para_single_FE
        if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
    } else {
        if (!(method %in% c('ZINB', 'B', 'B0', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model does not have random effects implementation!\n')
        if (ZINB != 'ZINB1') stop('Currently only ZINB1 is supported!\n')
        if (is.null(subject)) warning('subject is not supplied. Fixed effects model will be used instead!\n')
        perform_differential_analysis_para_single <- perform_differential_analysis_para_single_RE
    }
    
    df <- data.obj$meta.dat
    grp <- df[, grp.name]
    
    ind <- !is.na(grp)
    #data.obj <- subset_data(data.obj, ind)
    grp <- grp[ind]
    df <- df[ind, ]
    
    if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
        data.obj$abund.list[['Species']] <- data.obj$otu.tab
        rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
                                                             data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
    }
    
    dep <- colSums(data.obj$otu.tab)
    diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
    
    if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
        cat("Signficant sequencing depth confounding!\n")
        cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
        cat("There may be potential residual sequence depth confounding!\n")
    }
    
    pv.list <- qv.list <-  fc.list <- fc.lc.list <- fc.uc.list <- met.list <- list()
    res.final <- NULL
    
    if (norm == 'Precalculated') {
        dep <- data.obj$size.factor
    }
    if (norm == 'GMPR') {
        dep <- GMPR(data.obj$abund.list[[norm.level]], intersect.no)
    }
    
    if (norm == 'TSS') {
        dep <- colSums(data.obj$abund.list[[norm.level]])
    }
    
    ldep <- log(dep)
    
    for (LOI in taxa.levels) {
        cat(LOI, "\n")
        
        taxon.ct <- data.obj$abund.list[[LOI]]
        
        
        if (winsor == TRUE) {
            # Addressing the outlier (97% percent) or at least one outlier
            
            taxon.ct.p <- t(t(taxon.ct) / dep)
            taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
                cutoff <- quantile(x, winsor.qt)
                x[x >= cutoff] <- cutoff
                x
            }
            )
            # column/row switch
            taxon.ct <- t(round(taxon.ct.p * dep))
            
        }
        
        prop <- t(t(taxon.ct) / colSums(taxon.ct))
        if (!is.null(minp)) {
            prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
            taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
        }
        
        if (!is.null(medianp)) {
            nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
            prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
            taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
        }
        
        pv.vec <-  fc.vec <- fc.lc.vec <- fc.uc.vec <- met.vec <- conv.vec <- NULL
        obj <- NULL
        for (taxon in rownames(taxon.ct)) {
            cat('.')
            taxon.abund <- taxon.ct[taxon, ]
            
            ######## Logistic regression ###############
            if (method == 'B0') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B0', LRT, ...)) 
            if (method == 'B') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B', LRT, ...)) 
            if (method == 'QB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
            ######## Overdispersed Poisson regression #########	
            if (method == 'OP') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
            ######## Negative binomial regression #########
            if (method == 'NB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
            ######## Zeroinflated negbinomial regression 1 ########
            if (method == 'ZINB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
            
            # Adpative 0 selects OP and QB based on the zero proportion (Not optimal)
            if (method == 'Adaptive0') {
                temp <- mean(as.numeric(taxon.abund != 0))
                
                if (temp > zerop.cutoff) {
                    error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
                } else {
                    error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
                }
            }
            
            # Adpative 1 selects NB and ZIB based on AIC
            if (method == 'Adaptive1') {
                error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
                error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
                if (class(error1) != 'try-error' & class(error2) != 'try-error') {
                    if (obj1$aic < obj2$aic) {
                        obj <- obj1
                    } else {
                        obj <- obj2
                    }
                    error <- error1
                } else {
                    # pv == 0 indicates some problems in fitting
                    if (class(error1) != 'try-error' & obj1$pv != 0) {
                        obj <- obj1
                        error <- error1
                    } else {
                        if (class(error2) != 'try-error' & obj1$pv != 0) {
                            obj <- obj2
                            error <- error2
                        } else {
                            error <- error2
                        }
                    }
                }	
            }
            
            # Adaptive 2 starts with NB model, if it fails, it switches ZINB
            if (method == 'Adaptive2') {
                error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
                
                if (class(error1) == 'try-error' | obj1$pv == 0) {
                    error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
                    if (class(error2) != 'try-error') {
                        obj <- obj2
                        error <- error2
                    } else {
                        error <- error2
                    }
                } else {
                    error <- error1
                    obj <- obj1
                }
                
            }
            
            # Random Effects model 
            # ZINB, B, NB, Adpative1 is implemented based on glmmADMB
            
            # Set P value NA for those not makes sense
            if (class(error) == "try-error" | abs(obj$lfc) > 100) {
                obj$pv <- obj$lfc <- obj$lfc.lci <- obj$lfc.uci <- obj$method <- NA
            }
            
            
            pv.vec <- rbind(pv.vec, obj$pv)
            fc.vec <- rbind(fc.vec, obj$lfc)
            fc.lc.vec <- rbind(fc.lc.vec, obj$lfc.lci)
            fc.uc.vec <- rbind(fc.uc.vec, obj$lfc.uci)
            met.vec <- rbind(met.vec, obj$method)
            
        }
        cat('\n')
        
        qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
        
        rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(fc.uc.vec) <- rownames(fc.lc.vec) <- rownames(met.vec) <- rownames(prop)
        colnames(pv.vec) <- 'Pvalue'
        colnames(qv.vec) <- 'Qvalue'
        colnames(met.vec) <- 'Method'
        
        
        pv.list[[LOI]] <- pv.vec
        qv.list[[LOI]] <- qv.vec
        fc.list[[LOI]] <- fc.vec
        fc.lc.list[[LOI]] <- fc.lc.vec
        fc.uc.list[[LOI]] <- fc.uc.vec
        met.list[[LOI]] <- met.vec
        
        res <- cbind(pv.vec, qv.vec, fc.vec, fc.lc.vec, fc.uc.vec, met.vec)
        rownames(res) <- rownames(prop)
        #write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
        
        if (mt.method == 'fdr') {
            res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
        }
        if (mt.method == 'raw') {
            res.final <- rbind(res.final, res[ as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
        }
    }
    
    if (!is.null(res.final)) {
        colnames(res.final) <- colnames(res)
        res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
        #write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
    }
    return(list(pv.list=pv.list, qv.list=qv.list, fc.list=fc.list, fc.uc.list=fc.uc.list, fc.lc.list=fc.lc.list, met.list=met.list))
}


# =================================================== #
#           perform_differential_analysis_para_FE()
# =================================================== #
# The following function performs the taxonomic differential analysis
# this is a monster of a function that does a LOT with lots of options.
# 
perform_differential_analysis_para_single_FE <- function (taxon.abund, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
    # ldep: log depth (size factor)
    if (!is.null(adj.name)) {
        if (sum(grepl(grp.name, c(adj.name)))) {
            stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
        }
    }
    
    if (!is.null(subject)) {
        warnings('Fixed effects model will ignore the subject variable! Please use randome effects model!\n')
    }
    if (LRT & method == 'OP') warning('Overdispersed Poisson does not support LRT! Wald test used!\n')
    
    if (is.null(adj.name)) {
        grp.name.adj.name <- grp.name
    } else {
        grp.name.adj.name <- paste(grp.name, '+', adj.name)
    }
    if (method == 'NB') {
        m1.nb <- glm.nb(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')), data = df)
        if (LRT) {
            m0.nb <- update(m1.nb, as.formula(paste('. ~ . -',  grp.name)))
            code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary, m0.conv=m0.nb$converged, m0.bound=m0.nb$boundary)
            pv.nb <- anova(m1.nb, m0.nb)['Pr(Chi)'][2, ]
            method <- paste(method, 'LRT')
        } else {
            code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary)	
            pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
            method <- paste(method, 'Wald')
        }
        
        aic.nb <- summary(m1.nb)$aic
        
        coef.nb <- coef(m1.nb)
        fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
        
        ci.nb <- confint.default(m1.nb)
        obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
        
        if (is.vector(obj)) {
            fc.lc.nb <- obj[1]
            fc.uc.nb <- obj[2]
        } else {
            fc.lc.nb <- obj[, 1]
            fc.uc.nb <- obj[, 2]
            names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
            names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
        }
        return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
    }
    if (method == 'B') {
        taxon.abund2 <- as.numeric(taxon.abund != 0)
        m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=binomial)
        if (LRT) {
            m0.b <- update(m1.b, as.formula(paste('. ~ . -',  grp.name)))
            code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary, m0.conv=m0.b$converged, m0.bound=m0.b$boundary)
            pv.b <- pchisq(2 * (logLik(m1.b) - logLik(m0.b)), df = df.residual(m0.b) - df.residual(m1.b), lower.tail=FALSE)
            method <- paste(method, 'LRT')
        } else {
            code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
            pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
            method <- paste(method, 'Wald')
        }
        
        aic.b <- summary(m1.b)$aic
        coef.b <- coef(m1.b)
        fc.b <- coef.b[grep(grp.name, names(coef.b))]
        
        ci.b <- confint.default(m1.b)
        obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
        
        if (is.vector(obj)) {
            fc.lc.b <- obj[1]
            fc.uc.b <- obj[2]
        } else {
            fc.lc.b <- obj[, 1]
            fc.uc.b <- obj[, 2]
            names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
            names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
        }
        return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
    }
    
    # Rev: 2016_09_13 add 'QB', No likelihood ratio test
    if (method == 'QB') {
        taxon.abund2 <- as.numeric(taxon.abund != 0)
        m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=quasibinomial)
        
        code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
        pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
        method <- paste(method, 'Wald')
        
        
        aic.b <- summary(m1.b)$aic
        coef.b <- coef(m1.b)
        fc.b <- coef.b[grep(grp.name, names(coef.b))]
        
        ci.b <- confint.default(m1.b)
        obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
        
        if (is.vector(obj)) {
            fc.lc.b <- obj[1]
            fc.uc.b <- obj[2]
        } else {
            fc.lc.b <- obj[, 1]
            fc.uc.b <- obj[, 2]
            names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
            names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
        }
        return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
    }
    
    if (method == 'OP') {
        # No LRT
        m1.op <- glm(as.formula(paste('taxon.abund ~', grp.name.adj.name)), offset=ldep, data = df, family=quasipoisson)
        code <- list(m1.conv=m1.op$converged, m1.bound=m1.op$boundary)
        
        # pv.op <- pchisq(2 * (logLik(m1.op) - logLik(m0.op)), df = df.residual(m0.op) - df.residual(m1.op), lower.tail=FALSE) # LRT not applicable
        coef.op <- coef(m1.op)		
        pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = grep(grp.name, names(coef.op)))$result$chi2['P']
        method <- paste(method, 'Wald')
        fc.op <- coef.op[grep(grp.name, names(coef.op))]
        
        ci.op <- confint.default(m1.op)
        obj <- ci.op[grep(grp.name, rownames(ci.op)), ]
        
        if (is.vector(obj)) {
            fc.lc.op <- obj[1]
            fc.uc.op <- obj[2]
        } else {
            fc.lc.op <- obj[, 1]
            fc.uc.op <- obj[, 2]
            names(fc.lc.op) <- paste(names(fc.lc.op), '2.5%')
            names(fc.uc.op) <- paste(names(fc.uc.op), '97.5%')
        }
        return(list(method=method, pv=pv.op, lfc=fc.op, lfc.lci=fc.lc.op, lfc.uci=fc.uc.op, aic=NULL, code=code))
    }
    
    if (method == 'ZINB0') {
        m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')),
                            data = df, dist = "negbin", EM = TRUE)
        if (LRT) {
            if (is.null(adj.name)) {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep)')),
                                    data = df, dist = "negbin", EM = TRUE)
            } else {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep)')),
                                    data = df, dist = "negbin", EM = TRUE)
            }
            code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
            # LRT
            pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
            method <- paste(method, 'LRT')
        } else {
            code <- list(m1.conv=m1.zinb$converged)	
            pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
            method <- paste(method, 'Wald')
        }
        
        aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
        
        coef.zinb <- coef(m1.zinb)
        fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
        
        ci.zinb <- confint.default(m1.zinb)
        obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
        
        if (is.vector(obj)) {
            fc.lc.zinb <- obj[1]
            fc.uc.zinb <- obj[2]
        } else {
            fc.lc.zinb <- obj[, 1]
            fc.uc.zinb <- obj[, 2]
            names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
            names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
        }
        return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
    }
    
    
    if (method == 'ZINB1') {
        m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) | ldep')),
                            data = df, dist = "negbin", EM = TRUE)
        if (LRT) {
            if (is.null(adj.name)) {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
                                    data = df, dist = "negbin", EM = TRUE)
            } else {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) | ldep')),
                                    data = df, dist = "negbin", EM = TRUE)
            }
            code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
            # LRT
            pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
            method <- paste(method, 'LRT')
        } else {
            code <- list(m1.conv=m1.zinb$converged)	
            pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
            method <- paste(method, 'Wald')
        }
        
        aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
        
        coef.zinb <- coef(m1.zinb)
        fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
        
        ci.zinb <- confint.default(m1.zinb)
        obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
        
        if (is.vector(obj)) {
            fc.lc.zinb <- obj[1]
            fc.uc.zinb <- obj[2]
        } else {
            fc.lc.zinb <- obj[, 1]
            fc.uc.zinb <- obj[, 2]
            names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
            names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
        }
        return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
    }
    
    
    if (method == 'ZINB2') {
        m2.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) |', grp.name.adj.name, '+ ldep')),
                            data = df, dist = "negbin", EM = TRUE)
        if (LRT) {
            if (is.null(adj.name)) {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
                                    data = df, dist = "negbin", EM = TRUE)
            } else {
                m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) |', adj.name, ' + ldep')),
                                    data = df, dist = "negbin", EM = TRUE)
            }
            
            code <- list(m1.conv=m2.zinb$converged, m0.conv=m0.zinb$converged)
            # LRT
            pv2.zinb <-  pchisq(2 * (logLik(m2.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m2.zinb), lower.tail=FALSE)
            method <- paste(method, 'LRT')
        } else {
            code <- list(m2.conv=m2.zinb$converged)	
            pv2.zinb <- wald.test(b = coef(m2.zinb), Sigma = vcov(m2.zinb), Terms = grep(grp.name, names(coef(m2.zinb))))$result$chi2['P']
            method <- paste(method, 'Wald')
        }
        
        aic2.zinb <- -2 * logLik(m2.zinb) + 2 * (m2.zinb$n - m2.zinb$df.residual)
        
        coef.zinb <- coef(m2.zinb)
        fc2.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
        
        ci.zinb <- confint.default(m2.zinb)
        obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
        
        if (is.vector(obj)) {
            fc2.lc.zinb <- obj[1]
            fc2.uc.zinb <- obj[2]
        } else {
            fc2.lc.zinb <- obj[, 1]
            fc2.uc.zinb <- obj[, 2]
            names(fc2.lc.zinb) <- paste(names(fc2.lc.zinb), '2.5%')
            names(fc2.uc.zinb) <- paste(names(fc2.uc.zinb), '97.5%')
        }
        return(list(method=method, pv=pv2.zinb, lfc=fc2.zinb, lfc.lci=fc2.lc.zinb, lfc.uci=fc2.uc.zinb, aic=aic2.zinb, code=code))
    }
    
}

# =================================================== #
#           abundance_list_create()
# =================================================== #
# this function is designed to help transform the raw data into the
# form usable in the Jun Chen functions above.
# mydata is an object of raw counts and a grouping variable
# group.Var is a grouping variable.
#
# return a matrix of the observed counts in the groups
# across all individuals.
abundance_list_create <- function(mydata,Group.Var)
{
  N <- ncol(mydata) - 1
  output <- matrix(nrow=length(unique(Group.Var)), ncol = N)
  i <- 1
  for(i in 1:N){
    person.counts <- aggregate(mydata[,(1+i)], by=list(Group.Var), FUN=sum)
    output[,i] <- person.counts[,2]
  }
  rownames(output) <- person.counts[,1]
  colnames(output) <- colnames(mydata)[2:(N+1)]
  return(output)
}

