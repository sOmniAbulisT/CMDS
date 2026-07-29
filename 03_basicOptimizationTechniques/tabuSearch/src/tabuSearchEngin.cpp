#include <Rcpp.h>
#include <deque>
#include <algorithm>

using namespace Rcpp; 

// Move Attribute

struct MoveAttribute 
{
  int i; 
  int j; 
  int moveType; 
  
  // reload
  bool operator == (const MoveAttribute& other) const
  {
    return (i == other.i && j == other.j && moveType == other.moveType) ||
           (i == other.j && j == other.i && moveType == other.moveType)
  }
};

//' @export
//[[Rcpp::export]]
List tsEngine(IntegerVector initState, 
              Function objFunc, 
              int maxIterations, 
              int tabuTenure) {
  
  // --------------------------------------------------------------------------
  // [Initialization Phase]
  // --------------------------------------------------------------------------
  IntegerVector currentPoint = clone(initState);
  double currentEnergy = as<double>(objFunc(currentPoint));
  
  IntegerVector bestPoint = clone(currentPoint);
  double bestEnergy = currentEnergy;
  
  std::deque<MoveAttribute> tabuList;
  
  NumericVector energyHistory(maxIterations + 1);
  energyHistory[0] = currentEnergy;
  
  int vectorLength = currentPoint.size();
  
  // --------------------------------------------------------------------------
  // [Main Optimization Meta-Loop]
  // --------------------------------------------------------------------------
  for (int iter = 1; iter <= maxIterations; ++iter) {
    
    IntegerVector bestNeighbor;
    double bestNeighborEnergy = R_PosInf; 
    MoveAttribute bestMoveOfIteration = {-1, -1, -1};
    bool foundValidNeighbor = false;
    
    for (int i = 0; i < vectorLength - 1; ++i) {
      for (int j = i + 1; j < vectorLength; ++j) {
        
        // 
        IntegerVector ptest = clone(currentPoint);
        std::swap(ptest[i], ptest[j]);
        
        // 
        double ptestEnergy = as<double>(objFunc(ptest));
        
        // 
        MoveAttribute currentMove = {i, j, 0}; 
        
        //
        auto it = std::find(tabuList.begin(), tabuList.end(), currentMove);
        bool isTabu = (it != tabuList.end());
        
        bool isStandardPath = (!isTabu && (ptestEnergy < bestNeighborEnergy || !foundValidNeighbor));
        bool isAspirationPath = (ptestEnergy < bestEnergy); 
        
        if (isStandardPath || isAspirationPath) {
          bestNeighbor = clone(ptest);
          bestNeighborEnergy = ptestEnergy;
          bestMoveOfIteration = currentMove;
          foundValidNeighbor = true;
        }
      }
    }
    
    // --------------------------------------------------------------------------
    // [State Guard & Deadlock Protection]
    // --------------------------------------------------------------------------
    if (!foundValidNeighbor) {
      Rcout << "Iteration " << iter << ": Trapped in neighborhood deadlock (Empty Set). Graceful halt triggered.\n";
      energyHistory = energyHistory[Range(0, iter - 1)];
      break;
    }
    
    // --------------------------------------------------------------------------
    // [Forced State Transition]
    // --------------------------------------------------------------------------
    currentPoint = clone(bestNeighbor);
    currentEnergy = bestNeighborEnergy;
    energyHistory[iter] = currentEnergy;
    
    if (currentEnergy < bestEnergy) {
      bestPoint = clone(currentPoint);
      bestEnergy = currentEnergy;
    }
    
    // --------------------------------------------------------------------------
    // [Memory Struct Dynamics (FIFO Update)]
    // --------------------------------------------------------------------------
    tabuList.push_back(bestMoveOfIteration);
    
    if (tabuList.size() > static_cast<size_t>(tabuTenure)) {
      tabuList.pop_front();
    }
  }
  
  // --------------------------------------------------------------------------
  // [Output Packaging]
  // --------------------------------------------------------------------------
  return List::create(
    Named("bestPoint")  = bestPoint,
    Named("bestEnergy") = bestEnergy,
    Named("history")    = energyHistory
  );
}