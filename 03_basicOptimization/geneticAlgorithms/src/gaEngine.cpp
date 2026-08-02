#include <Rcpp.h>
#include <vector>
#include <random>
#include <numeric>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ==============================================================================
// 1. Data Structure Definitions
// ==============================================================================
struct Individual {
  std::vector<int> genes;
  double fitness;
};

// Calculated the fitness of population
double calculateStdDev (const std::vector<Individual>& pop)
{
  double sum = 0.0; 
  for (const auto& ind : pop) sum += ind.fitness; 
  double mean = sum / pop.size(); 
  
  double accum = 0.0;
  for (const auto& ind : pop) accum += (ind.fitness - mean) * (ind.fitness - mean); 
  return std::sqrt(accum / (pop.size() - 1)); 
}

// ==============================================================================
// 2. GA operators
// ==============================================================================

// Selection
inline Individual selection (const std::vector<Individual>& population, 
                             std::mt19937& gen, 
                             std::uniform_int_distribution<>& d_tour)
{
  Individual winner = population[d_tour(gen)]; 
  for (int k = 0; k < 2; ++k)
  {
    Individual competitor = population[d_tour(gen)]; 
    if (competitor.fitness < winner.fitness)
    {
      winner = competitor; 
    }
  }
  
  return winner; 
}

// Crossover (single-point)
inline void crossover (Individual& child1, Individual& child2, 
                       int chromLength, std::mt19937& gen, 
                       std::bernoulli_distribution& d_cross, 
                       std::uniform_int_distribution<>& d_cut)
{
  if (d_cross(gen))
  {
    int cutPoint = d_cut(gen); 
    for (int j = cutPoint; j < chromLength; ++j)
    {
      std::swap(child1.genes[j], child2.genes[j]); 
    }
  }
}

// Mutation
inline void mutation (Individual& child, int chromLength, 
                      std::mt19937& gen, std::bernoulli_distribution& d_mut)
{
  for (int j = 0; j < chromLength; ++j)
  {
    if (d_mut(gen))
    {
      child.genes[j] ^= 1; 
    }
  }
}

// ==============================================================================
// 3. Genetic Algorithm main loop
// ==============================================================================

//' @export
// [[Rcpp::export]]
List GA (Function objFunc, 
         int popSize, 
         int chromLength, 
         double pCrossover, 
         double pMutation, 
         int maxGenerations, 
         double minStdDev)
{
  std::random_device rd; 
  std::mt19937 gen(rd()); 
  std::bernoulli_distribution d_init(0.5); 
  std::bernoulli_distribution d_cross(pCrossover); 
  std::bernoulli_distribution d_mut(pMutation); 
  std::uniform_int_distribution<> d_tour(0, popSize - 1); 
  std::uniform_int_distribution<> d_cut(1, chromLength - 1); 
  
  // Initialization
  std::vector<Individual> population(popSize); 
  Individual globalBest; 
  globalBest.fitness = R_PosInf;
  
  for (int i = 0; i < popSize; ++i)
  {
    population[i].genes.resize(chromLength); 
    for (int j = 0; j < chromLength; ++j)
    {
      population[i].genes[j] = d_init(gen) ? 1 : 0; 
    }
    population[i].fitness = as<double>(objFunc(wrap(population[i].genes))); 
    
    if (population[i].fitness < globalBest.fitness)
    {
      globalBest = population[i]; 
    }
  }
  
  NumericVector historyBest(maxGenerations); 
  int currentGen = 0; 
  double currentStdDev = R_PosInf; 
  
  // Main Loop
  while (currentGen < maxGenerations && currentStdDev > minStdDev)
  {
    std::vector<Individual> nextPopulation; 
    nextPopulation.reserve(popSize); 
    
    while (nextPopulation.size() < static_cast<size_t>(popSize))
    {
      // selection
      Individual p1 = selection(population, gen, d_tour); 
      Individual p2 = selection(population, gen, d_tour); 
      
      Individual child1 = p1; 
      Individual child2 = p2; 
      
      // crossover
      crossover(child1, child2, chromLength, gen, d_cross, d_cut); 
      
      // mutation
      mutation(child1, chromLength, gen, d_mut); 
      mutation(child2, chromLength, gen, d_mut); 
      
      nextPopulation.push_back(child1); 
      if (nextPopulation.size() < static_cast<size_t>(popSize))
      {
        nextPopulation.push_back(child2); 
      }
    }
    
    // Evaluation
    for (int i = 0; i < popSize; ++i)
    {
      nextPopulation[i].fitness = as<double>(objFunc(wrap(nextPopulation[i].genes))); 
      if (nextPopulation[i].fitness < globalBest.fitness)
      {
        globalBest = nextPopulation[i]; 
      }
    }
    
    // alternative of generations
    population = std::move(nextPopulation); 
    currentStdDev = calculateStdDev(population); 
    historyBest[currentGen] = globalBest.fitness; 
    currentGen++; 
  }
  
  // mistake prevention
  if (currentGen < maxGenerations)
  {
    historyBest = historyBest[Range(0, currentGen - 1)]; 
  }
  
  return List::create(
    Named("bestGenes") = globalBest.genes,
    Named("bestFitness") = globalBest.fitness,
    Named("history") = historyBest,
    Named("generationsRun") = currentGen
  );
}
