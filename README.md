# cvtree
Predictive score for phylogenetic trees

Installation:
```
devtools::install_github("rizbicki/cvtree")  
```

Usage:
```
set.seed(0)
my_tree<- ape::rtree(n=26,tip.label=LETTERS)
chars_test <- replicate(2,ape::rTraitDisc(my_tree,rate = 0.5))
cvtree(my_tree,chars_test,root="A")
```
