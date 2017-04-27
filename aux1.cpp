#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// [[Rcpp::export]]
NumericMatrix convertSBtoNormal(NumericMatrix vmat, 
                                int ncol, int nrow, 
                                NumericVector prod) {
  NumericMatrix res(nrow,ncol);

  for(int j=0; j<ncol;j++){
    res(_,j)=vmat(_,j)*prod;    
    prod=prod*(1-vmat(_,j));
  }

  return (res);
}

