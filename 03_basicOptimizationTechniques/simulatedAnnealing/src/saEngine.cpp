#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Simulated Annealing Core Algorithm
//' 
// [[Rcpp::export]]
List SA (NumericVector initState, 
         Function objFunc, 
         Function perturbFunc, 
         double initialTemp, 
         double alpha, 
         double minTemp, 
         int markovIter)
{
  // Initilization states and energy
  NumericVector currentPoint = clone(initState); 
  NumericVector bestPoint = clone(initState); 
  
  // Source the functions from R
  double currentEnergy = as<double>(objFunc(currentPoint)); 
  double bestEnergy = currentEnergy; 
  
  double currentTemp = initialTemp; 
  
  // Record the global best history
  NumericVector historyBest; 
  historyBest.push_back(bestEnergy); 
  
  // Cooling Schedules
  while (currentTemp > minTemp)
  {
    // Markov Chain
    for (int i = 0; i < markovIter; ++i)
    {
      // 1. Generate the new state
      NumericVector newPoint = perturbFunc(currentPoint); 
      
      // 2. Calculated the new energy and diff. of energy
      double newEnergy = as<double>(objFunc(newPoint)); 
      double deltaEnergy = newEnergy - currentEnergy; 
      
      // 3. Metropolis criterion
      if (deltaEnergy <= 0)
      {
        currentPoint = clone(newPoint); 
        currentEnergy = newEnergy; 
        
        if (currentEnergy < bestEnergy)
        {
          bestPoint = clone(currentPoint); 
          bestEnergy = currentEnergy; 
        }
      }
      else
      {
        double probAccept = std::exp(-deltaEnergy / currentTemp); 
        if (R::runif(0.0, 1.0) < probAccept)
        {
          currentPoint = clone(newPoint); 
          currentEnergy = newEnergy; 
        }
      }
    }
    
    historyBest.push_back(bestEnergy); 
    currentTemp = currentTemp / alpha; 
  }
  
  return List::create(
    Named("bestState") = bestPoint,
    Named("bestEnergy") = bestEnergy,
    Named("history") = historyBest
  );
}