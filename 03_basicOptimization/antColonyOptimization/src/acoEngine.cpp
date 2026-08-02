#include <Rcpp.h>
#include <vector>
#include <random>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace Rcpp;
using namespace std;

struct Ant 
{
  vector<int> path; 
  double fitness;   
};

// [[Rcpp::export]]
List ACO (NumericMatrix rDistMatrix, int numAnts, double evaporationRate, int maxGenerations)
{
  int numNodes = rDistMatrix.nrow();
  
  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> distProb(0.0, 1.0);
  uniform_int_distribution<int> distNode(0, numNodes - 1);
  
  vector<vector<double>> distanceMatrix(numNodes, vector<double>(numNodes));
  for (int i = 0; i < numNodes; ++i) {
    for (int j = 0; j < numNodes; ++j) {
      distanceMatrix[i][j] = rDistMatrix(i, j);
    }
  }
  
  // Step 1: Initialization
  vector<vector<double>> pheromoneMatrix(numNodes, vector<double>(numNodes, 1.0));
  vector<int> globalBestPath;
  double globalBestFitness = numeric_limits<double>::infinity();
  NumericVector historyBest(maxGenerations);
  
  // Main Generations Loop
  for (int iter = 0; iter < maxGenerations; ++iter)
  {
    vector<Ant> colony(numAnts);
    
    // Step 2 & 3: Path Construction
    for (int i = 0; i < numAnts; ++i)
    {
      vector<bool> visited(numNodes, false); // Tabu List
      
      int currentNode = distNode(gen);
      colony[i].path.push_back(currentNode);
      visited[currentNode] = true;
      
      for (int step = 1; step < numNodes; ++step)
      {
        vector<double> probabilities;
        vector<int> candidates;
        double probabilitySum = 0.0;
        
        for (int nextNode = 0; nextNode < numNodes; ++nextNode)
        {
          if (!visited[nextNode])
          {
            candidates.push_back(nextNode);
            double tau = pheromoneMatrix[currentNode][nextNode];
            probabilities.push_back(tau);
            probabilitySum += tau;
          }
        }
        
        // Roulette-Wheel Selection
        double randomPoint = distProb(gen) * probabilitySum;
        double cumulativeProbability = 0.0;
        int selectedNode = candidates.back(); 
        
        for (size_t c = 0; c < probabilities.size(); ++c)
        {
          cumulativeProbability += probabilities[c];
          if (randomPoint <= cumulativeProbability)
          {
            selectedNode = candidates[c];
            break;
          }
        }
        
        colony[i].path.push_back(selectedNode);
        visited[selectedNode] = true;
        currentNode = selectedNode;
      }
      
      colony[i].path.push_back(colony[i].path[0]); // Return to start
      
      double totalDistance = 0.0;
      for (size_t s = 0; s < colony[i].path.size() - 1; ++s) {
        totalDistance += distanceMatrix[colony[i].path[s]][colony[i].path[s + 1]];
      }
      colony[i].fitness = totalDistance;
      
      if (colony[i].fitness < globalBestFitness)
      {
        globalBestFitness = colony[i].fitness;
        globalBestPath = colony[i].path;
      }
    }
    
    // Step 4: Pheromone Update
    for (int r = 0; r < numNodes; ++r) {
      for (int c = 0; c < numNodes; ++c) {
        pheromoneMatrix[r][c] *= (1.0 - evaporationRate);
      }
    }
    
    double currentBestFitness = numeric_limits<double>::infinity();
    double currentWorstFitness = -numeric_limits<double>::infinity();
    int bestAntIndex = -1;
    
    for (int i = 0; i < numAnts; ++i)
    {
      if (colony[i].fitness < currentBestFitness) {
        currentBestFitness = colony[i].fitness;
        bestAntIndex = i;
      }
      if (colony[i].fitness > currentWorstFitness) {
        currentWorstFitness = colony[i].fitness;
      }
    }
    
    if (bestAntIndex != -1 && currentWorstFitness > 0)
    {
      double deltaTau = (currentWorstFitness - currentBestFitness) / currentWorstFitness;
      const auto& bestPath = colony[bestAntIndex].path;
      
      for (size_t s = 0; s < bestPath.size() - 1; ++s)
      {
        int fromNode = bestPath[s];
        int toNode = bestPath[s + 1];
        pheromoneMatrix[fromNode][toNode] += deltaTau;
        pheromoneMatrix[toNode][fromNode] += deltaTau;
      }
    }
    
    historyBest[iter] = globalBestFitness;
  }
  
  // Convert 0-based to 1-based index for R
  IntegerVector rBestPath(globalBestPath.size());
  for (size_t i = 0; i < globalBestPath.size(); ++i) {
    rBestPath[i] = globalBestPath[i] + 1; 
  }
  
  return List::create(
    Named("bestPath")     = rBestPath,
    Named("bestFitness")  = globalBestFitness,
    Named("history")      = historyBest
  );
}