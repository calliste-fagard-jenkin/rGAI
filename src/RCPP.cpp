#include <Rcpp.h>
using namespace Rcpp;

// @export
//[[Rcpp::export]]
NumericVector stopover_for_loop(NumericVector betas, NumericVector phi){
  // purpose : Performs the evaluation of the lambda_j equation in Matechou et 
  //           al (2014), given the beta estimates and the constant probability
  //           phi that an individual remains from one occasion to the next
  // inputs  : betas - The vector of betas, one for each occasion 
  //           phi   - The retention probability of an indivdual.
  // output  : A vector of the same length as betas, with the density expected
  //           at each occasion.
  
  int nT = betas.size();
  NumericVector a_func(nT);
  a_func[0] = betas[0];
  // copy the object with a loop, or beta value updates with a_func:
  for (int i = 0; i < nT; i++){
    a_func[i] = betas[i];
  }
  

  // same for loop as supplementary code, rewritten in C++:
  for(int j = 2; j <= nT; j++){
    for(int b = 1; b < j; b++){
      double update =  a_func[j - 1] + betas[b - 1]*pow(phi[0], j - b);
      a_func[j - 1] = update;
    }
  }
  
  return a_func;
}// stopover_for_loop

// @export
//[[Rcpp::export]]
NumericVector weightedSelection(NumericVector x, NumericVector probs){
  // purpose : Takes a vector of values, and a vector of probs and returns
  //           a vector which contains each element of x with the probability
  //           specific by probs
  // inputs  : x      - the vector of values which we want to select from
  //           breaks - an ordered vector where the ith entry is the probability
  //                    of keeping entry i of x.
  // output  : A numeric vector of selected values from x
  
  // declare vars:
  int xLen;
  xLen = x.size();
  NumericVector storage(xLen);
  NumericVector deviate = runif(xLen);
  int counter;
  counter = 0;
  
  for (int i=0; i<xLen; i++){
    if (deviate[i]<probs[i]){
      storage[counter] = x[i];
      counter++;
    }
  }
  
  NumericVector output(counter);
  for (int i=0; i < counter; i++){
    output[i] = storage[i];
  }
  
  return output;
}//weightedSelection
