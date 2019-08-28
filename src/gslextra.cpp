#include "../include/gslextra.h"


void gsl_matrix_flip(gsl_matrix * dest,int height,int width)
{
    int i,j;
    
    for(i=0;i<height;i++)
    {
        for(j=0;j<i;j++)
        {
            gsl_matrix_set(dest,i,j,gsl_matrix_get(dest,j,i));
        }
    }
}