#include <iostream>
#include "matrix.hpp"
#include <vector>
#include <numeric>
#include <string>
 
 mat_Tridiagonal::mat_Tridiagonal(size_t M): M(M)
{
    mat = new double *[M]; 
    for (std::size_t i=0; i<M; i++) 
    {
        mat[i] = new double[M]; 
        std::fill(&mat[i][0], &mat[i][M], 0);
}