# Protein Folding by Ant Colony Optimization (ACO) algorithm
The protein folding problem is one of the most fundamental and challenging problems in molecular biology. In this project, I used Python to solve this problem by Ant Colony Optimization (ACO) algorithm. This algorithm was implemented after reading many research papers. References to those are mentioned.

## Requirements

- Python 3.x
- 
- Numpy

- Matplotlib

### Protein Folding (PF)
Protein folding is the process by which a linear chain of amino acids folds into a three-dimensional structure that determines its function. The 2D HP model is a simplified version of protein folding that represents the protein as a sequence of amino acids in a two-dimensional lattice, where each amino acid is either hydrophobic (H) or polar (P). In this model, the hydrophobic amino acids tend to cluster together in the interior of the protein, while the polar amino acids tend to be on the exterior.

The goal of the 2D HP model is to find the lowest energy conformation of the protein, which corresponds to the most stable and biologically functional structure. The energy of a conformation is determined by the number of hydrophobic contacts between neighboring amino acids in the lattice. Two hydrophobic amino acids are considered to be in contact if they are adjacent on the lattice and their side chains face each other.

The problem of finding the lowest energy conformation of a protein in the 2D HP model is known to be NP-hard, which means that it is computationally intractable for large proteins. However, several heuristic algorithms have been developed to approximate the solution, such as Monte Carlo simulations, genetic algorithms, and ant colony optimization.

The 2D HP model is useful for studying the general principles of protein folding and for exploring the effects of different factors, such as temperature, pressure, and mutations, on the stability and function of proteins. It is also a useful tool for designing artificial proteins with specific structures and functions for applications in medicine, biotechnology, and nanotechnology.

### Ant Colony Optimization (ACO)
Ant Colony Optimization (ACO) is a metaheuristic algorithm inspired by the behavior of ant colonies. It is used to solve complex optimization problems that involve finding the optimal path in a graph or network. The algorithm is based on the observation that ants are able to find the shortest path between their nest and food sources by leaving a trail of pheromones that guides other ants towards the food.

In the ACO algorithm, a colony of artificial ants is created and allowed to explore the solution space by constructing and evaluating candidate solutions. Each ant follows a probabilistic decision rule that takes into account the amount of pheromone deposited on each edge and the distance to the destination. The ants deposit pheromone on the edges they traverse, with the amount of pheromone being proportional to the quality of the solution. As a result, good solutions are reinforced by the deposition of more pheromone, making them more likely to be chosen in subsequent iterations.

Over time, the pheromone trail evaporates, which means that solutions that are no longer optimal are gradually eliminated from consideration. This helps the algorithm to avoid getting stuck in local optima and encourages exploration of the solution space.


## Output
Here, cyan represents Polar amino acids, and, black represents Hydrophobic amino acids.

![Output](https://github.com/DivyanshPandey99/Protein-Folding-ACO/blob/main/Screenshot.png)
