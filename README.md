# Kemeny-and-Conquer
Computation of Kemeny's constant using censored Markov Chains

## Hutch++

Some of the tests make use of the stochastic trace estimator called Hutch++. The repository with related code is added as a git submodule. If you use that version of the algorithm for estimating the Kemeny constant, you should also cite the original article:

```{bibtex}
@incollection {MR4537958,
    AUTHOR = {Meyer, Raphael A. and Musco, Cameron and Musco, Christopher
              and Woodruff, David P.},
     TITLE = {Hutch++: optimal stochastic trace estimation},
 BOOKTITLE = {Symposium on {S}implicity in {A}lgorithms ({SOSA})},
     PAGES = {142--155},
 PUBLISHER = {[Society for Industrial and Applied Mathematics (SIAM)],
              Philadelphia, PA},
      YEAR = {2021},
      ISBN = {978-1-61197-649-6},
   MRCLASS = {68W20 (15A15)},
  MRNUMBER = {4537958},
}
```

## hm-toolbox

One of the implemented version of the recursive algorithm for the approximation
of the Kemeny constant uses the hm-toolbox by Massei, Robol and Kressner. If you
use that version of the algorithm for estmating the Kemeny constant, you should
also cite the original article:

```{bibtex}
@article {MR4082285,
    AUTHOR = {Massei, Stefano and Robol, Leonardo and Kressner, Daniel},
     TITLE = {hm-toolbox: {MATLAB} software for {HODLR} and {HSS} matrices},
   JOURNAL = {SIAM J. Sci. Comput.},
  FJOURNAL = {SIAM Journal on Scientific Computing},
    VOLUME = {42},
      YEAR = {2020},
    NUMBER = {2},
     PAGES = {C43--C68},
      ISSN = {1064-8275,1095-7197},
   MRCLASS = {65F99 (15A03 65F05 65F08 65F10 65F50 65Y15)},
  MRNUMBER = {4082285},
MRREVIEWER = {Bruno\ Carpentieri},
       DOI = {10.1137/19M1288048},
       URL = {https://doi.org/10.1137/19M1288048},
}
```
