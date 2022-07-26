test_that("chrom.fitness.score works", {
    data(case)
    data(dad)
    data(mom)
    case <- as.matrix(case)
    dad <- as.matrix(dad)
    mom <- as.matrix(mom)
    comp <- mom + dad - case
    case.comp.diff <- case != comp
    case.minus.comp <- case - comp
    storage.mode(case.minus.comp) <- "integer"
    both.one.mat <- case == 1 & comp == 1
    case2.mat <- case == 2
    case0.mat <- case == 0
    comp2.mat <- comp == 2
    comp0.mat <- comp == 0
    library(Matrix)
    block.ld.mat <- as.matrix(bdiag(list(
        matrix(rep(TRUE, 25^2), nrow = 25),
        matrix(rep(TRUE, 25^2), nrow = 25),
        matrix(rep(TRUE, 25^2), nrow = 25),
        matrix(rep(TRUE, 25^2), nrow = 25)
    )))
    weight.lookup <- vapply(seq_len(8), function(x) 2^x, 1)
    storage.mode(weight.lookup) <- "integer"
    f <- chrom.fitness.score(
        case, comp, case.comp.diff, c(51, 52, 76, 77),
        case.minus.comp, both.one.mat,
        block.ld.mat, weight.lookup,
        case2.mat, case0.mat, comp2.mat,
        comp0.mat
    )
    expect_equal(f$fitness_score, 18.189961)
})
