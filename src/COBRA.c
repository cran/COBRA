#include <stdio.h>
#include <stdlib.h>
#include <R.h>

void COBRA(int *n1,
	   int *n2,
	   double *train,
	   double *test,
	   int *nmachines,
	   double *alpha,
	   double *poids,
	   double *e)
{
  int i,j,k;
  int seuil = *alpha*(*nmachines);
  int *vres = NULL;
  vres = calloc(*n1,sizeof(int));
  for(i=0;i<*n2;i++)
    {
      for(j=0;j<*n1;j++)
	{
	  vres[j] = 0;
	  for(k=0;k<*nmachines;k++)
	    {
	      if(fabs(train[*n1*k+j] - test[*n2*k+i])<=*e)
		{
		  vres[j]++;
		}
	    }
	  if(vres[j]>=seuil)
	    {
	      poids[*n2*j+i] = 1; 
	    }
	}
    }
  free(vres);
}
