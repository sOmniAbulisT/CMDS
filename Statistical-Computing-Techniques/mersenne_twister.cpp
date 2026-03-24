// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cstdint>
#include <vector>

using namespace Rcpp; 

// [[Rcpp::export]]
NumericVector gen_unif(int n, uint32_t seed = 12345){
  NumericVector res(n); 
  
  std::vector<uint32_t> mt(624); 
  int index = 624;
   
  mt[0] = seed; 
  for(int i = 1; i < 624; ++i){
    mt[i] = (1812433253UL * (mt[i-1]^(mt[i - 1] >> 30)) + i); 
  }
  
  for (int k = 0; k < n; ++k) {
    
    if (index >= 624) {
      for (int i = 0; i < 624; ++i) {
        uint32_t y = (mt[i] & 0x80000000UL) + (mt[(i + 1) % 624] & 0x7fffffffUL);
        mt[i] = mt[(i + 397) % 624] ^ (y >> 1);
        
        if (y % 2 != 0) {
          mt[i] = mt[i] ^ 0x9908b0dfUL;
        }
      }
      index = 0; 
    }
    
    uint32_t y = mt[index++];
    
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    
    res[k] = y / 4294967296.0; 
  }
  
  return res;
}