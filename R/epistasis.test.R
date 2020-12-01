#' A function to run a test of the null hypothesis that a collection of SNPs do not exhibit epistasis, conditional
#' upon the observed marginal associations.
#'
#' This function runs a permutation based test of the null hypothesis that a collection of SNPs do not exhibit epistasis,
#' conditional upon the observed marginal associations.
#'
#' @param snp.cols An integer vector specifying the columns in the input data containing the SNPs to be tested.
#' @param preprocessed.list The initial list produced by function \code{preprocess.genetic.data}.
#' @param n.permutes The number of permutations on which to base the test. Defaults to 1000.
#' @param n.different.snps.weight The number by which the number of different SNPs between a case and complement or unaffected sibling
#'  is multiplied in computing the family weights. Defaults to 2.
#' @param n.both.one.weight The number by which the number of SNPs equal to 1 in both the case and complement or unaffected sibling
#' is multiplied in computing the family weights. Defaults to 1.
#' @param weight.function.int An integer used to assign family weights. Specifically, we use \code{weight.function.int} in a function that takes the weighted sum
#' of the number of different SNPs and SNPs both equal to one as an argument, denoted as x, and
#' returns a family weight equal to \code{weight.function.int}^x. Defaults to 2.
#' @param recode.threshold For a given SNP, the minimum test statistic required to recode and recompute the fitness score using recessive coding. Defaults to 3.
#' See the GADGET paper for specific details.
#' @param bp.param The BPPARAM argument to be passed to bplapply when estimating marginal disease associations for each SNP.
#'  If using a cluster computer, this parameter needs to be set with care. See \code{BiocParallel::bplapply} for more details
#' @return A list of thee elements:
#' \describe{
#'  \item{pval}{The p-value of the test.}
#'  \item{obs.fitness.score}{The fitness score from the observed data}
#'  \item{perm.fitness.scores}{A vector of fitness scores for the permuted datasets.}
#' }
#' @examples
#'
#' data(case)
#' data(dad)
#' data(mom)
#' data(snp.annotations)
#' library(Matrix)
#' block.ld.mat <- as.matrix(bdiag(list(matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25),
#'                               matrix(rep(TRUE, 25^2), nrow = 25))))
#'
#' pp.list <- preprocess.genetic.data(case, father.genetic.data = dad,
#'                                mother.genetic.data = mom,
#'                                block.ld.mat = block.ld.mat)
#'
#' run.gadgets(pp.list, n.chromosomes = 5, chromosome.size = 3,
#'        results.dir = "tmp", cluster.type = "interactive",
#'        registryargs = list(file.dir = "tmp_reg", seed = 1300),
#'        n.top.chroms = 5, n.islands = 8, island.cluster.size = 4,
#'        n.migrations = 2)
#'
#' combined.res <- combine.islands('tmp', snp.annotations, pp.list)
#'
#' top.snps <- as.vector(t(combined.res$unique.results[1, 1:3]))
#' set.seed(10)
#' epi.test.res <- epistasis.test(top.snps, pp.list)
#'
#' unlink('tmp', recursive = TRUE)
#' unlink('tmp_reg', recursive = TRUE)
#'
#' @importFrom data.table rbindlist setkey setorder `:=`
#' @importFrom BiocParallel bplapply bpparam
#' @export

epistasis.test <- function(snp.cols, preprocessed.list, n.permutes = 1000,
                     n.different.snps.weight = 2, n.both.one.weight = 1,
                     weight.function.int = 2, recode.threshold = 3,
                     bp.param = bpparam()) {

    ### pick out the target columns in the pre-processed data ###
    original.col.numbers <- preprocessed.list$original.col.numbers
    target.snps <- which(original.col.numbers %in% snp.cols)

    ### grab LD mat ###
    block.ld.mat <- preprocessed.list$block.ld.mat
    target.block.ld.mat <- block.ld.mat[target.snps, target.snps]

    # stop if all SNPs are in LD
    if (all(target.block.ld.mat)){

        stop("cannot run test, all SNPs are in LD")

    }

    ### compute the weight lookup table ###
    max.sum <- max(n.different.snps.weight, n.both.one.weight)*length(snp.cols)
    weight.lookup <- vapply(seq_len(max.sum), function(x) weight.function.int^x, 1)

    ### grab pre-processed case and complement data ###
    case.genetic.data <- preprocessed.list$case.genetic.data
    complement.genetic.data <- preprocessed.list$complement.genetic.data

    ### Compute matrices of differences between cases and complements ###
    case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
    case.comp.different <- case.minus.comp != 0

    ### Compute matrix indicating whether both the case and control have the same ###
    ### number of copies of the minor allele ###
    both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1

    ### compute matrices indicating whether case carries two copies of the risk allele
    case2.mat <- case.genetic.data == 2
    case0.mat <- case.genetic.data == 0

    ### compute fitness score and find the families with the full risk set ###
    fitness.score <- chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different,
                                    target.snps, case.minus.comp, both.one.mat,
                                    block.ld.mat, weight.lookup, case2.mat, case0.mat,
                                    n.different.snps.weight, n.both.one.weight,
                                    recode.threshold, epi.test = TRUE)

    ### uses informative families already, so no need to recompute the observed fitness score here ###
    obs.fitness.score <- fitness.score$fitness.score*min(abs(fitness.score$sum.dif.vecs))

    ### restrict the data to informative families ###
    case.inf <- case.genetic.data[fitness.score$inf.families, target.snps]
    comp.inf <- complement.genetic.data[fitness.score$inf.families, target.snps]
    n.families <- length(fitness.score$inf.families)

    ### identify ld blocks for the target snps ###
    remaining.snps <- seq_len(length(target.snps))
    ld.blocks <- list()
    while (length(remaining.snps) > 0){

        old.ld.blocks <- unlist(ld.blocks)
        starting.snp <- min(remaining.snps)
        snp.ld.block <- target.block.ld.mat[ , starting.snp]
        ld.blocks <- c(ld.blocks, list(which(snp.ld.block)))
        if (length(old.ld.blocks > 0)){

            remaining.snps <- remaining.snps[!snp.ld.block[-old.ld.blocks]]

        } else {

            remaining.snps <- remaining.snps[!snp.ld.block]

        }

    }

    ### loop over ld blocks and shuffle rows of the case/complement data ###
    permutes <- seq_len(n.permutes)
    permuted.data <- lapply(permutes, function(permute){

        case.permuted <- case.inf
        comp.permuted <- comp.inf
        for (block.snps in ld.blocks){

            this.order <- sample(seq_len(n.families), n.families)
            case.permuted[ , block.snps] <- case.permuted[this.order , block.snps]
            comp.permuted[ , block.snps] <- comp.permuted[this.order , block.snps]

        }
    return(list(case = case.permuted, complement = comp.permuted))

    })

    ### loop over permuted datasets and compute fitness scores
    perm.fitness.scores <- unlist(bplapply(permuted.data,

        function(permute, target.snps, target.block.ld.mat, weight.lookup, n.different.snps.weight,
                 n.both.one.weight, recode.threshold){

        ### grab case and complement data ###
        case.genetic.data <- permute$case
        complement.genetic.data <- permute$complement

        ### Compute matrices of differences between cases and complements ###
        case.minus.comp <- sign(as.matrix(case.genetic.data - complement.genetic.data))
        case.comp.different <- case.minus.comp != 0

        ### Compute matrix indicating whether both the case and control have the same ###
        ### number of copies of the minor allele ###
        both.one.mat <- complement.genetic.data == 1 & case.genetic.data == 1

        ### compute matrices indicating whether case carries two copies of the risk allele
        case2.mat <- case.genetic.data == 2
        case0.mat <- case.genetic.data == 0

        ### compute fitness score ###
        fitness.score <- chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different,
                                             seq_len(length(target.snps)), case.minus.comp, both.one.mat,
                                             target.block.ld.mat, weight.lookup, case2.mat, case0.mat,
                                             n.different.snps.weight, n.both.one.weight,
                                             recode.threshold)
        fs <- fitness.score$fitness.score*min(abs(fitness.score$sum.dif.vecs))
        return(fs)

    }, target.snps  = target.snps, target.block.ld.mat = target.block.ld.mat, weight.lookup = weight.lookup,
    n.different.snps.weight =  n.different.snps.weight, n.both.one.weight = n.both.one.weight,
    recode.threshold = recode.threshold, BPPARAM = bp.param))

    ### compute p-value ###
    pval <- sum(perm.fitness.scores >= obs.fitness.score)/(n.permutes + 1)

    ### return result ###
    return(list(pval = pval, obs.fitness.score = obs.fitness.score, perm.fitness.scores = perm.fitness.scores))

}