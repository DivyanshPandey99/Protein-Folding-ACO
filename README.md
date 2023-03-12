# Protein Folding by Ant Colony Optimization (ACO) algorithm


### Ant Colony Optimization (ACO)
Ant Colony Optimization (ACO) is a metaheuristic algorithm inspired by the behavior of ant colonies. It is used to solve complex optimization problems that involve finding the optimal path in a graph or network. The algorithm is based on the observation that ants are able to find the shortest path between their nest and food sources by leaving a trail of pheromones that guides other ants towards the food.

In the ACO algorithm, a colony of artificial ants is created and allowed to explore the solution space by constructing and evaluating candidate solutions. Each ant follows a probabilistic decision rule that takes into account the amount of pheromone deposited on each edge and the distance to the destination. The ants deposit pheromone on the edges they traverse, with the amount of pheromone being proportional to the quality of the solution. As a result, good solutions are reinforced by the deposition of more pheromone, making them more likely to be chosen in subsequent iterations.

Over time, the pheromone trail evaporates, which means that solutions that are no longer optimal are gradually eliminated from consideration. This helps the algorithm to avoid getting stuck in local optima and encourages exploration of the solution space.


## Output
Here, cyan represents Polar amino acids, and, black represents Hydrophobic amino acids.

![Output](https://github.com/DivyanshPandey99/Protein-Folding-ACO/blob/main/Screenshot.png)
