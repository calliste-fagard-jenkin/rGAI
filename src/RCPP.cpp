#include <Rcpp.h>
using namespace Rcpp;

// @export
// [[Rcpp::export]]
NumericVector probs_link(NumericVector x){
  // purpose : Takes a vector of parameter guesses on the real line and maps
  //           them to the probability of belonging to any given component
  // input   : x - The vector of guesses
  // output  : A vector of probabilities that sum to 1
  
  // Apply the link function to all the values:
  int xLen = x.size();
  NumericVector linked(xLen);
  linked = plogis(x);
  
  // Create a variable to store the results:
  NumericVector output(xLen + 1);
  output[0] = linked[0];
  
  // Store the remaining proportion of the prob left to assign:
  double left(1);
  left = (1 - linked[0]);
  
  // Loop to determine the probability of being in each group:
  double prob(1);
  for (int i = 1; i < xLen; i++){
    prob = left*linked[i];
    output[i] = prob;
    left = left - prob;
  }
  
  output[xLen] = left;
  return output;
}//probs_link

// @export
//[[Rcpp::export]]
NumericVector vector_to_counts(NumericVector x, NumericVector breaks){
  // purpose : Takes a vector of values, and a vector of breaks and returns
  //           a vector which indicates how many of the observations are in 
  //           each of the intervals defined by breaks.
  // inputs  : x      - the vector of values which we want to count through
  //           breaks - an ordered (increasing) vector of break points, defining
  //                    the boundaries of the length(breaks) - 1 intervals.
  // output  : A vector of length breaks - 1, where each entry is the count in 
  //           x of values between the ith and ith + 1 breakpoint.
  
  // initialise output vector:
  int intervalNum;
  int xLen;
  intervalNum = breaks.size() - 1;
  xLen = x.size();
  NumericVector output(intervalNum);
  
  for (int i = 0; i < xLen; i++){
    // for each element in x:
    for (int j = 0; j < intervalNum; j++){
      // for each break:
      if ((x[i] > breaks[j]) & (x[i] <= breaks[j+1])){
        output[j]++;
      }
    }
  }
  
  return output;
}//vector_to_counts

//' Rescale a numeric vector so that its entries sum to one
//' 
//' @param x A numeric vector.
//' @export
//[[Rcpp::export]]
NumericVector sum_to_one(NumericVector x){
  // purpose : Takes a vector of counts and returns a vector of the same
  //           length where each entry is the proportion of the total counts
  //           of that entry in x
  
  int xLen = x.size();
  
  double sum = 0;
  for (int i=0; i < xLen; i++){sum = sum + x[i];}
  
  NumericVector output(xLen);
  for (int i=0; i < xLen; i++){output[i] = x[i]/sum;}
  return output;
}// sum_to_one

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
NumericVector means_link(NumericVector means){
  // purpose : Takes a vector of the type (mean1, diff1, diff2, ...) and turns
  //           it into a vector of the type (mean1, mean2, means3, ...)
  // inputs  : means - The vector of inputs
  // output  : A vector of the same length as means
  
  int meansLen = means.size();
  NumericVector output(meansLen);
  output[0] = exp(means[0]);
  
  for (int i = 1; i < meansLen; i++){
    output[i] = exp(means[i]) + output[i - 1];
  }
  
  return output;
}

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
