Version 1.1.1
+ Updated functions to improve memory usage. 
+ Related to point 1, changed how linkage information is specified in 'preprocess.genetic.data'.
  See the GADGETS vignette for details. 
+ Also in 'preprocess.genetic.data', the function no longer re-codes input genotypes 
  based on the apparent minor allele. If, for a particular SNP, an analyst wants the 
  genotype coded '1' to represent 1 copy of the minor allele and '2' to represent 
  two copies of the minor allele, they need to code the genotypes in that way 
  before inputting the data. With that said, the functions in this package are 
  invariant to genotype coding choice (and always have been) and results 
  will be the same regardless. It is therefore completely unnecessary to do any 
  re-coding, which is the reason we have removed that step from the function. 
+ Added in experimental functionality for GxGxE search and maternal-fetal effects

Version 1.1.0
+ accepted to Bioconductor

Version 0.99.0
 + pushed to github
