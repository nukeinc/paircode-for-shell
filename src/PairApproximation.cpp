#include "PairApproximation.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

namespace NuclearPairing {

PairApproximation::PairApproximation(const std::vector<SingleParticleState>& states,
                                     int particle_number,
                                     double pairing_strength)
    : states_(states)
    , particle_number_(particle_number)
    , pairing_strength_(pairing_strength)
{
    if (states.empty()) {
        throw std::invalid_argument("States vector cannot be empty");
    }
    if (particle_number <= 0) {
        throw std::invalid_argument("Particle number must be positive");
    }
    if (particle_number % 2 != 0) {
        throw std::invalid_argument("Particle number must be even for BCS pairing");
    }
    if (pairing_strength < 0) {
        throw std::invalid_argument("Pairing strength must be non-negative");
    }
}

double PairApproximation::quasiparticleEnergy(double epsilon, double lambda, double gap) {
    double xi = epsilon - lambda;
    return std::sqrt(xi * xi + gap * gap);
}

double PairApproximation::vSquared(double epsilon, double lambda, double gap) {
    if (gap < 1e-15) {
        // In the limit of zero gap, use step function
        return (epsilon < lambda) ? 1.0 : 0.0;
    }
    double xi = epsilon - lambda;
    double E_k = quasiparticleEnergy(epsilon, lambda, gap);
    return 0.5 * (1.0 - xi / E_k);
}

double PairApproximation::gapEquation(double gap, double lambda) const {
    if (gap < 1e-15) {
        return -gap;  // Trivial solution check
    }

    double sum = 0.0;
    for (const auto& state : states_) {
        double E_k = quasiparticleEnergy(state.energy, lambda, gap);
        // Sum over pair states (degeneracy/2 pairs per level)
        sum += state.degeneracy / (2.0 * E_k);
    }

    // Gap equation: 1/G = sum_k (degeneracy_k / (2 * E_k))
    // Rearranged: gap - G * gap * sum = 0
    return 1.0 - pairing_strength_ * sum;
}

double PairApproximation::particleNumberEquation(double gap, double lambda) const {
    double sum = 0.0;
    for (const auto& state : states_) {
        double v2 = vSquared(state.energy, lambda, gap);
        sum += state.degeneracy * v2;
    }
    return sum - particle_number_;
}

double PairApproximation::findChemicalPotential(double gap, double tol) const {
    // Find bounds for chemical potential
    double e_min = states_.front().energy;
    double e_max = states_.back().energy;
    for (const auto& state : states_) {
        e_min = std::min(e_min, state.energy);
        e_max = std::max(e_max, state.energy);
    }

    // Expand bounds slightly
    double margin = 0.2 * (e_max - e_min) + 10.0;
    double lambda_min = e_min - margin;
    double lambda_max = e_max + margin;

    // Bisection to find lambda where N(lambda) = particle_number
    int max_iter = 200;
    for (int iter = 0; iter < max_iter; ++iter) {
        double lambda_mid = 0.5 * (lambda_min + lambda_max);
        double N_residual = particleNumberEquation(gap, lambda_mid);

        if (std::abs(N_residual) < tol) {
            return lambda_mid;
        }

        if (N_residual > 0) {
            // Too many particles, lower chemical potential
            lambda_max = lambda_mid;
        } else {
            // Too few particles, raise chemical potential
            lambda_min = lambda_mid;
        }
    }

    return 0.5 * (lambda_min + lambda_max);
}

BCSResult PairApproximation::solve(double initial_gap,
                                    double tolerance,
                                    int max_iterations) {
    BCSResult result;
    result.converged = false;
    result.iterations = 0;

    double gap = initial_gap;
    double lambda = Utils::calculateFermiEnergy(states_, particle_number_);

    // Iterative solution of coupled gap and particle number equations
    for (int iter = 0; iter < max_iterations; ++iter) {
        result.iterations = iter + 1;

        // Find chemical potential for current gap
        lambda = findChemicalPotential(gap, tolerance * 0.01);

        // Calculate new gap from gap equation
        double gap_sum = 0.0;
        for (const auto& state : states_) {
            double E_k = quasiparticleEnergy(state.energy, lambda, gap);
            if (E_k > 1e-15) {
                gap_sum += state.degeneracy / (2.0 * E_k);
            }
        }

        double new_gap = (gap_sum > 1e-15) ? pairing_strength_ * gap * gap_sum : 0.0;

        // Check convergence
        if (std::abs(new_gap - gap) < tolerance) {
            result.converged = true;
            gap = new_gap;
            break;
        }

        // Damped update for stability
        gap = 0.5 * gap + 0.5 * new_gap;

        // Prevent gap from becoming too small (collapse to normal state)
        if (gap < tolerance * 0.1) {
            gap = 0.0;
            result.converged = true;
            break;
        }
    }

    // Final chemical potential
    lambda = findChemicalPotential(gap, tolerance * 0.01);

    // Store results
    result.gap = gap;
    result.chemical_potential = lambda;

    // Calculate occupation probabilities and BCS factors
    result.occupation.resize(states_.size());
    result.u_factors.resize(states_.size());
    result.v_factors.resize(states_.size());

    for (size_t i = 0; i < states_.size(); ++i) {
        double v2 = vSquared(states_[i].energy, lambda, gap);
        result.occupation[i] = v2;
        result.v_factors[i] = std::sqrt(v2);
        result.u_factors[i] = std::sqrt(1.0 - v2);
    }

    // Calculate pairing energy
    // E_pair = -Delta^2 / G (simplified form)
    if (pairing_strength_ > 1e-15 && gap > 1e-15) {
        result.pairing_energy = -gap * gap / pairing_strength_;
    } else {
        result.pairing_energy = 0.0;
    }

    return result;
}

namespace Utils {

std::vector<SingleParticleState> generateHarmonicOscillatorStates(double hw, int max_shell) {
    std::vector<SingleParticleState> states;

    for (int N = 0; N <= max_shell; ++N) {
        // N = 2n + l
        for (int l = N % 2; l <= N; l += 2) {
            int n = (N - l) / 2;

            // j = l + 1/2 and j = l - 1/2 (for l > 0)
            std::vector<double> j_values;
            if (l > 0) {
                j_values.push_back(l + 0.5);
                j_values.push_back(l - 0.5);
            } else {
                j_values.push_back(0.5);
            }

            for (double j : j_values) {
                // Energy without spin-orbit: E = hw * (N + 3/2)
                // Simple spin-orbit correction: -kappa * (l, s) term
                double kappa = 0.05 * hw;  // Typical spin-orbit strength
                double spin_orbit = (l > 0) ? 
                    kappa * ((j > l) ? l / 2.0 : -(l + 1) / 2.0) : 0.0;

                double energy = hw * (N + 1.5) + spin_orbit;
                double degeneracy = 2.0 * j + 1.0;

                states.emplace_back(energy, degeneracy, n, l, j);
            }
        }
    }

    // Sort by energy
    std::sort(states.begin(), states.end(),
              [](const SingleParticleState& a, const SingleParticleState& b) {
                  return a.energy < b.energy;
              });

    return states;
}

std::vector<SingleParticleState> generateWoodsSaxonStates(int A, int Z, bool is_proton, int max_n) {
    // Simplified Woods-Saxon energies based on empirical fits
    // These are approximate values for illustration
    std::vector<SingleParticleState> states;

    // Woods-Saxon parameters (for reference and potential future enhancement)
    // V0: potential depth, hw: oscillator frequency derived from mass number
    double V0 = 51.0;  // MeV

    // Coulomb correction for protons
    double coulomb = is_proton ? 0.7 * Z / std::pow(A, 1.0/3.0) : 0.0;

    // Generate states for each shell
    // Using simplified formula: E_nlj ~ hw * (2n + l) + spin-orbit
    double hw = 41.0 / std::pow(A, 1.0/3.0);  // MeV

    for (int n = 0; n <= max_n; ++n) {
        for (int l = 0; l <= max_n - n; ++l) {
            std::vector<double> j_values;
            if (l > 0) {
                j_values.push_back(l + 0.5);
                j_values.push_back(l - 0.5);
            } else {
                j_values.push_back(0.5);
            }

            for (double j : j_values) {
                // Approximate energy formula
                int N = 2 * n + l;
                double energy_base = hw * (N + 1.5) - V0;

                // Spin-orbit splitting
                double spin_orbit_strength = 0.4 * hw;  // Typical value
                double spin_orbit = (l > 0) ?
                    spin_orbit_strength * ((j > l) ? l : -(l + 1)) / (2.0 * l + 1.0) : 0.0;

                double energy = energy_base + spin_orbit + coulomb;
                double degeneracy = 2.0 * j + 1.0;

                states.emplace_back(energy, degeneracy, n, l, j);
            }
        }
    }

    // Sort by energy
    std::sort(states.begin(), states.end(),
              [](const SingleParticleState& a, const SingleParticleState& b) {
                  return a.energy < b.energy;
              });

    return states;
}

double calculateFermiEnergy(const std::vector<SingleParticleState>& states,
                            int particle_number) {
    // Sort states by energy
    std::vector<SingleParticleState> sorted_states = states;
    std::sort(sorted_states.begin(), sorted_states.end(),
              [](const SingleParticleState& a, const SingleParticleState& b) {
                  return a.energy < b.energy;
              });

    // Fill states up to particle number
    int filled = 0;
    double fermi_energy = sorted_states.front().energy;

    for (const auto& state : sorted_states) {
        int can_fill = static_cast<int>(state.degeneracy);
        if (filled + can_fill >= particle_number) {
            fermi_energy = state.energy;
            break;
        }
        filled += can_fill;
        fermi_energy = state.energy;
    }

    return fermi_energy;
}

} // namespace Utils

} // namespace NuclearPairing
