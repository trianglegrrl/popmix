# DATES Ground Truth Reference Data

## Simulation Parameters

The example data is from a simulated admixture event with:
- **Source 1**: French (CEU proxy) - 20% ancestry
- **Source 2**: Yoruba (YRI) - 80% ancestry
- **Admixture time (lambda)**: 50 generations
- **Sample size**: 1 admixed individual (Simulation) and 4 admixed individuals (Simulation.1)

## DATES Output (Simulation_French_Yoruba)

Reference from original DATES software v753:
- **Estimated generations**: 35.738
- **Standard error**: 6.346
- **True value**: 50 generations

Note: The estimate is approximately 2 SE below the true value, which is within expected variation for a single admixed individual.

## Output Files

- `Simulation.out`: LD decay data (distance_cM, covariance, z_score, correlation, n_pairs)
- `Simulation.jin`: Jackknife estimates per chromosome
- `Simulation.jout`: Final result (mean, SE)
- `Simulation.fit`: Exponential fit (distance_cM, observed, fitted, residual)

## Usage in Testing

The `Simulation.out` file provides reference LD decay values that our implementation should match.
The jackknife and fitting outputs can be used to verify our exponential fitting and SE estimation.
