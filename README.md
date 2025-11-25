# Nuclear Pair Approximation Code (paircode-for-shell)

A C++ implementation of the BCS (Bardeen-Cooper-Schrieffer) pairing approximation for nuclear structure calculations.

## Overview

This code implements the nuclear pair approximation, which describes pairing correlations between nucleons (protons or neutrons) in atomic nuclei. The BCS approximation is fundamental to understanding nuclear structure, particularly for even-even nuclei where pairs of like nucleons tend to couple to angular momentum zero.

## Features

- BCS gap equation solver with iterative convergence
- Calculation of pairing gap (Δ), chemical potential (λ), and pairing energy
- Occupation probabilities and BCS u,v factors
- Harmonic oscillator single-particle state generator
- Woods-Saxon potential state generator (simplified)
- Support for both proton and neutron calculations

## Building

### Requirements

- CMake 3.10 or higher
- C++17 compatible compiler (GCC, Clang, MSVC)

### Build Instructions

```bash
mkdir build
cd build
cmake ..
make
```

### Running Tests

```bash
cd build
ctest --output-on-failure
```

## Usage

```bash
./paircode [options]
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `-A <number>` | Mass number | 120 |
| `-Z <number>` | Proton number | 50 |
| `-N <number>` | Neutron number | A-Z |
| `-G <value>` | Pairing strength (MeV) | 0.25 |
| `-hw <value>` | Oscillator parameter (MeV) | auto |
| `-shells <num>` | Maximum shell number | 5 |
| `-proton` | Calculate for protons | |
| `-neutron` | Calculate for neutrons | default |
| `-v` | Verbose output | |
| `-h` | Show help | |

### Examples

```bash
# Calculate neutron pairing for Sn-120
./paircode -A 120 -Z 50 -G 0.25 -neutron

# Calculate proton pairing for Pb-208 with verbose output
./paircode -A 208 -Z 82 -G 0.1 -proton -shells 6 -v

# Custom calculation with specific parameters
./paircode -A 100 -Z 44 -G 0.3 -hw 10.0 -shells 4 -neutron -v
```

## Theory Background

### BCS Approximation

The BCS (Bardeen-Cooper-Schrieffer) theory, originally developed for superconductivity, is applied to nuclear physics to describe pairing correlations. The key equations are:

**Gap Equation:**
```
1/G = Σₖ Ωₖ / (2Eₖ)
```

where:
- G is the pairing strength
- Ωₖ is the degeneracy of state k
- Eₖ = √[(εₖ - λ)² + Δ²] is the quasiparticle energy
- εₖ is the single-particle energy
- λ is the chemical potential
- Δ is the pairing gap

**Number Equation:**
```
N = Σₖ Ωₖ vₖ²
```

where vₖ² = ½(1 - (εₖ - λ)/Eₖ) is the occupation probability.

### Physical Interpretation

- **Pairing gap (Δ):** Energy required to break a nucleon pair
- **Occupation probability (v²):** Probability that a state is occupied
- **BCS factors (u, v):** Related to the wave function of the paired system

## API Reference

### Classes

#### `NuclearPairing::PairApproximation`

Main class for performing BCS calculations.

```cpp
PairApproximation(const std::vector<SingleParticleState>& states,
                  int particle_number,
                  double pairing_strength);

BCSResult solve(double initial_gap = 1.0,
                double tolerance = 1e-8,
                int max_iterations = 1000);
```

#### `NuclearPairing::SingleParticleState`

Structure representing a single-particle state.

```cpp
struct SingleParticleState {
    double energy;      // Single-particle energy (MeV)
    double degeneracy;  // Degeneracy (2j+1)
    int quantum_n;      // Principal quantum number
    int quantum_l;      // Orbital angular momentum
    double quantum_j;   // Total angular momentum
};
```

#### `NuclearPairing::BCSResult`

Structure containing the results of a BCS calculation.

```cpp
struct BCSResult {
    double gap;                         // Pairing gap (MeV)
    double chemical_potential;          // Chemical potential (MeV)
    double pairing_energy;              // Total pairing energy (MeV)
    std::vector<double> occupation;     // Occupation probabilities
    std::vector<double> u_factors;      // BCS u-factors
    std::vector<double> v_factors;      // BCS v-factors
    int iterations;                     // Convergence iterations
    bool converged;                     // Convergence status
};
```

### Utility Functions

```cpp
namespace NuclearPairing::Utils {
    // Generate harmonic oscillator single-particle states
    std::vector<SingleParticleState> generateHarmonicOscillatorStates(
        double hw, int max_shell);

    // Generate Woods-Saxon single-particle states
    std::vector<SingleParticleState> generateWoodsSaxonStates(
        int A, int Z, bool is_proton, int max_n = 5);

    // Calculate Fermi energy
    double calculateFermiEnergy(const std::vector<SingleParticleState>& states,
                                int particle_number);
}
```

## License

This project is provided for educational and research purposes.

## References

1. Bohr, A., & Mottelson, B. R. (1969). *Nuclear Structure, Vol. 1*. Benjamin.
2. Ring, P., & Schuck, P. (1980). *The Nuclear Many-Body Problem*. Springer.
3. Bardeen, J., Cooper, L. N., & Schrieffer, J. R. (1957). Theory of Superconductivity. *Physical Review*, 108(5), 1175.