# ncImpute
Matrix Completion via Non-Convex Regularization

This R package provides iterative algorithms for matrix completion based on non-convex regularization of the singular values. The main approach uses iterative MC+ thresholded singular value decompositions to impute the missing values, and has an "EM" flavor, in that at each iteration the matrix is completed with the current estimate. For large matrices there is a special sparse-matrix class named "Incomplete" that efficiently handles all computations. The package includes procedures for centering and scaling rows, columns or both.
