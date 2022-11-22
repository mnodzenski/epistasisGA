test_that("chrom.fitness.score works", {
    data(case)
    data(dad)
    data(mom)
    case <- as.matrix(case)
    dad <- as.matrix(dad)
    mom <- as.matrix(mom)
    comp <- mom + dad - case
    block.ld.vec <- cumsum(rep(25, 4))
    weight.lookup <- vapply(seq_len(8), function(x) 2^x, 1)
    storage.mode(weight.lookup) <- "integer"
    f <- chrom.fitness.score(case, comp,
                             c(51, 52, 76, 77),
                             block.ld.vec,
                             weight.lookup)
    expect_equal(f$fitness_score, 18.189961)
})
