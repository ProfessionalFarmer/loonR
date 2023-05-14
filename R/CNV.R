#' Categorical variable test
#'
#' @param sample1
#' @param sample2
#' @param f Default chisq.test, can be fisher.test
#'
#' @return
#' @export
#'
#' @examples
#' g1 = sample(c("Amp","N","Del"), 10, replace = T)
#' g2 = c( sample(c("Amp","N","Del"), 15, replace = T), sample(c("Del"), 5, replace = T) )
#' res = loonR::testCategoricalVariableFreq(g1, g2)
#' res$P
testCategoricalVariableFreq <- function(sample1, sample2,
                                        s1.name = "Sample1", s2.name = "Sample2",
                                        f = "chisq.test"){

        twoXN.table = list()
        twoXN.p = list()
        raw.test.res = list()
        ######### Prepare
        sample1 = sample1[!is.na(sample1)]
        sample2 = sample2[!is.na(sample2)]
        s1 = unlist(table(sample1))
        s2 = unlist(table(sample2))

        # in case the same variable in two independant samples
        all.names = unique(c(names(s1), names(s2)))

        if(length(all.names)==1){
                twoXN.p = data.frame(Round = all.names, P = 1, check.names = F)
                row.names(twoXN.p) = all.names

                res = list(p = twoXN.p)
                return( res )
        }

        tmp = rep(0, length(all.names))
        names(tmp) = all.names
        tmp[names(s1)] = as.vector(s1)
        s1 = tmp

        tmp = rep(0, length(all.names))
        names(tmp) = all.names
        tmp[names(s2)] = as.vector(s2)
        s2 = tmp

        chisq.table = rbind(s1, s2)
        rownames(chisq.table) = c(s1.name, s2.name)



        ###########
        overall.test = get(f)(chisq.table)
        cat("Using ", f, "\n --------------------\n")
        twoXN.p[["Overall"]] = overall.test$p.value
        raw.test.res[["Overall"]] = overall.test
        twoXN.table[["Overall"]] = chisq.table

        for(variable in all.names){

                n11 = sum(sample1 == variable)
                n21 = sum(sample2 == variable)

                n12 = length(sample1) - n11
                n22 = length(sample2) - n21

                chisq.table <- matrix(c(n11, n21, n12, n22),
                                      nrow = 2,
                                      dimnames = list(Sample = c(s1.name, s2.name),
                                                      Group = c(variable, paste0("Non-",variable))
                                      ))

                single.test = get(f)(chisq.table)

                twoXN.p[[variable]] = single.test$p.value
                raw.test.res[[variable]] = single.test
                twoXN.table[[variable]] = chisq.table

        }

        twoXN.p = unlist(twoXN.p)
        twoXN.p = data.frame(Round = names(twoXN.p), P = twoXN.p, check.names = F)
        row.names(twoXN.p) = twoXN.p$Round
        twoXN.p$P = as.numeric(twoXN.p$P)

        res = list(p = twoXN.p,
                   Count.table = twoXN.table,
                   RawTestResult = raw.test.res)
        res
}


