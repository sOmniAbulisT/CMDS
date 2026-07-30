#include <Rcpp.h>
#include <vector>
#include <random>
#include <algorithm>

using namespace Rcpp;
using namespace std;

// ==============================================================================
// 1. Define Particle Data Structure 
// ==============================================================================
struct Particle 
{
  vector<double> position; 
  vector<double> velocity; 
  
  vector<double> pBestPosition; 
  double pBestFitness; 
};

// ==============================================================================
// 2. Core algorithm for particle swarm
// ==============================================================================
// [[Rcpp::export]]
List PSO (Function objFunc, int popSize, int numVars, 
          double lowerBound, double upperBound, 
          int maxGenerations, double wMax, double wMin, 
          double c1, double c2, double vMax)
{
  random_device rd; 
  mt19937 gen(rd()); 
  uniform_real_distribution<double> distPos(lowerBound, upperBound); 
  uniform_real_distribution<double> distProb(0.0, 1.0); 
  
  vector<Particle> swarm(popSize); 
  vector<double> gBestPosition(numVars); 
  double gBestFitness = numeric_limits<double>::infinity(); // Minimize
  
  NumericVector historyBest(maxGenerations); 
  
  // Initialization
  for (int i = 0; i < popSize; ++i)
  {
    swarm[i].position.resize(numVars); 
    swarm[i].velocity.resize(numVars, 0.0); // V_0 = 0
    swarm[i].pBestPosition.resize(numVars); 
    
    for (int j = 0; j < numVars; ++j)
    {
      swarm[i].position[j] = distPos(gen);
      swarm[i].pBestPosition[j] = swarm[i].position[j]; 
    }
    
    //
    NumericVector rcppPos = wrap(swarm[i].position); 
    double currentFitness = as<double>(objFunc(rcppPos)); 
    swarm[i].pBestFitness = currentFitness; 
    
    // gBest
    if (currentFitness < gBestFitness)
    {
      gBestFitness = currentFitness; 
      gBestPosition = swarm[i].position; 
    }
  }
  
  // Iteration and Kinematics update
  for (int iter = 0; iter < maxGenerations; ++iter)
  {
    
    // Linear decay of inertia weight
    double w = wMax - ((wMax - wMin) / maxGenerations) * iter; 
    
    for (int i = 0; i < popSize; ++i)
    {
      // A. calculated new velocity and position
      for (int j = 0; j < numVars; ++j)
      {
        double r1 = distProb(gen); 
        double r2 = distProb(gen); 
        
        // final velocity
        swarm[i].velocity[j] = w * swarm[i].velocity[j] + 
          c1 * r1 * (swarm[i].gBestPosition[j] - swarm[i].position[j]) + 
          c2 * r2 * (gBestPosition[j] - swarm[i].position[j]); 
        
        // Velocity clamping
        swarm[i].velocity[j] = clamp(swarm[i].velocity[j], -vMax, vMax); 
        
        // Position update
        swarm[i].velocity[j] += swarm[i].velocity[j]; 
        
        // Boundary defense
        if (swarm[i].position[j] > upperBound || swarm[i].position[j] < lowerBound)
        {
          swarm[i].position[j] = clamp(swarm[i].position[j], lowerBound, upperBound); 
          swarm[i].position[j] = 0.0; 
        }
      }
      
      // B. evaluation
      NumericVector rcppPos = wrap(swarm[i].position); 
      double currentFitness = as<double>(objFunc(rcppPos)); 
      
      // C. Update pBest and gBest
      if (currentFitness < swarm[i].pBestFitness)
      {
        swarm[i].pBestFitness = currentFitness; 
        swarm[i].pBestPosition = swarm[i].position; 
        
        if (currentFitness < gBestFitness)
        {
          gBestFitness = currentFitness; 
          gBestPosition = swarm[i].position; 
        }
      }
    }
    
    historyBesy[iter] = gBestFitness;
  }
  
  return List::create(
    Named("bestPosition") = gBestPosition,
    Named("bestFitness")  = gBestFitness,
    Named("history")      = historyBest 
  );
}
