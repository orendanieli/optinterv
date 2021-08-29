
Thoughts for the next version:
==============================

1. Change the metric for ordered/unordered categorical variables (CDF distance is more suitable for continuous variables).

2. Validate that there is no multicollinearity problem with the controls (this throws an error).

3. Remove unnecessary documentation and update the package on CRAN

Nearest Neighbors Issues:

1. The confidence intervals are not expected to work when there is a small number of neighbors. That's because they are calculated with bootstrap and are only true asymptotically. We should report somehow the effective number of neighbors. E.g. 1/HHI for every row
2. potentially lambda should be set in a way that is proportional to sigma. If sigma is very small, you start putting very high weight on distance and less on Y. 
3. For some reason the package doesn't support very small sigmas. 
