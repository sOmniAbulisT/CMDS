# Particle Swarm Optimization 

## Overview
Particle Swarm Optimization (PSO) is a population-based, stochastic optimization algorithm inspired by the swarm intelligence observed in biological phenomena such as bird flocking and fish schooling. Unlike traditional genetic algorithms that rely on discrete crossover and mutation, PSO is inherently designed for **continuous real-valued search spaces**. Through the sharing of information among microscopic particles, the algorithm exhibits powerful "Emergence" and "Self-Organization," allowing it to efficiently locate global optima in high-dimensional spaces with a minimal memory footprint and small population size.

This module implements the standard PSO algorithm, featuring a linearly decreasing inertia weight mechanism and strict spatial boundary defense constraints.

---

## Mathematical Model

The core of the PSO engine is driven by two linear update equations for "Velocity" and "Position".

### 1. Velocity Update Equation
Determines the flight direction and step size of a particle for the next iteration. It is a superposition of three physical vectors: inertia, individual cognition, and social influence:
$$V_{j,t} = w_t V_{j,t-1} + c_1 r_1(P_j - X_{j,t-1}) + c_2 r_2(P_G - X_{j,t-1})$$
* **$w_t$**: Inertia Weight. Provides the momentum to maintain the current trajectory, ensuring global exploration capabilities.
* **$c_1, c_2$**: Cognitive & Social learning factors. These control the relative pull between a particle trusting its own memory versus blindly following the swarm.
* **$r_1, r_2$**: Independent uniform random numbers generated between $0 \sim 1$, ensuring the stochasticity and diversity of the optimization trajectory.
* **$P_j$**: The historical personal best position of particle $j$ (pBest).
* **$P_G$**: The historical global best position discovered by the entire swarm (gBest).

### 2. Position Update Equation
After the new velocity is calculated, a linear spatial translation is applied:
$$X_{j,t} = X_{j,t-1} + V_{j,t}$$

### 3. Linear Decay of Inertia Weight
To achieve a dynamic balance between early-stage "global exploration" and late-stage "local exploitation," the inertia weight $w$ decays linearly as the iteration count $i$ increases:
$$w_t = w_{\text{max}} - \left(\frac{w_{\text{max}} - w_{\text{min}}}{i_{\text{max}}}\right) i$$

---

## Repository Structure

This branch (`03_basicOptimizationTechniques/particleSwarmOptimization`) follows a modular architecture, 
separating the high-performance C++ computational core from the high-level R scripting environment.

```text
.
├── src/
│   ├── psoEngine.cpp      # The core C++ engine: defines the Particle struct, memory allocations, and kinematics update equations.
│   └── Makevars           # (Optional) Compiler flags for C++ optimization (e.g., -O3, -march=native).
├── main.R                 # Master controller: handles Rcpp compilation, objective function definition, hyperparameter initialization, and executes the C++ solver.
│ 
├── visualize.R            # Analytics script: reads the serialized output to generate convergence curves and trajectory plots using ggplot2.  
│           
└── README.md              # Algorithm documentation, mathematical models, and architectural overview.
```

---

## Algorithm Workflow

1. **Initialization:**
   * Define a lightweight swarm size $N$ (typically 20-30 particles).
   * Randomly scatter particle coordinates within the specified design space boundaries $[X^{(l)}, X^{(u)}]$.
   * Initialize all particle velocities to absolute zero ($V_{j,0} = 0$).
   * Evaluate the initial objective function fitness and establish the first generation of $P_j$ and $P_G$ memories.

2. **Evaluation & Memory Update:**
   * Calculate the objective function fitness $f(X)$ at the current position.
   * If a particle surpasses its personal historical record, update the cognitive memory $P_j$.
   * If any particle outperforms the swarm's current global record, a throne transition occurs; update the global memory $P_G$.

3. **Kinematics Update:**
   * Compute the new velocity $V_{j,t}$ based on the velocity equation.
   * **Velocity Clamping Defense:** If the calculated absolute velocity exceeds the predefined limit, clamp it to $\pm V_{\text{max}}$ to prevent severe overshooting and oscillation.
   * Update the spatial position $X_{j,t}$.
   * **Boundary Collision Defense:** If the new position violates the design space boundaries, force the coordinates back onto the boundary limits and neutralize the momentum (velocity) in that specific dimension.

4. **Termination:**
   * Continue the iteration loop until the maximum number of generations ($i_{\text{max}}$) is reached, or the spatial standard deviation of the swarm falls below a strict convergence tolerance (indicating all particles have condensed into a single optimum point).

---

## Implementation Architecture (C++ / Rcpp)
* **High-Performance Core:** The underlying physical kinematics are fully implemented in C++, utilizing `std::vector` for contiguous memory allocation to maximize CPU Cache hit rates. Boundary and velocity defenses are securely handled via `std::clamp`.
* **Stochastic Engine:** Random number generation is powered by the `std::mt19937` engine from the standard `<random>` library, guaranteeing high-quality statistical Monte Carlo simulation properties.
* **Seamless Interface:** Encapsulated as an Rcpp module, allowing the higher-level R environment to dynamically inject custom objective functions, define search boundaries, and orchestrate core hyperparameters without sacrificing compiled execution speed.