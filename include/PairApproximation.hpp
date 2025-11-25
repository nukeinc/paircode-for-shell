#ifndef PAIR_APPROXIMATION_HPP
#define PAIR_APPROXIMATION_HPP

#include <vector>
#include <cmath>
#include <stdexcept>
#include <functional>

namespace NuclearPairing {

/**
 * Structure representing a single-particle state in the nuclear shell model
 */
struct SingleParticleState {
    double energy;      // Single-particle energy in MeV
    double degeneracy;  // Degeneracy (2j+1) of the state
    int quantum_n;      // Principal quantum number
    int quantum_l;      // Orbital angular momentum
    double quantum_j;   // Total angular momentum

    SingleParticleState(double e = 0.0, double deg = 1.0, 
                        int n = 0, int l = 0, double j = 0.5)
        : energy(e), degeneracy(deg), quantum_n(n), quantum_l(l), quantum_j(j) {}
};

/**
 * Results from the BCS calculation
 */
struct BCSResult {
    double gap;                         // Pairing gap (Delta) in MeV
    double chemical_potential;          // Chemical potential (lambda) in MeV
    double pairing_energy;              // Total pairing energy in MeV
    std::vector<double> occupation;     // Occupation probabilities (v^2)
    std::vector<double> u_factors;      // BCS u-factors
    std::vector<double> v_factors;      // BCS v-factors
    int iterations;                     // Number of iterations to converge
    bool converged;                     // Whether the solution converged
};

/**
 * Class implementing the BCS (Bardeen-Cooper-Schrieffer) pairing approximation
 * for nuclear structure calculations.
 * 
 * The BCS approximation describes the pairing correlations between nucleons
 * (protons or neutrons) in the nuclear shell model.
 */
class PairApproximation {
public:
    /**
     * Constructor
     * @param states Vector of single-particle states
     * @param particle_number Number of particles (must be even for pairing)
     * @param pairing_strength Pairing interaction strength G (typically 0.1-0.5 MeV)
     */
    PairApproximation(const std::vector<SingleParticleState>& states,
                      int particle_number,
                      double pairing_strength);

    /**
     * Solve the BCS gap equation
     * @param initial_gap Initial guess for the gap parameter (MeV)
     * @param tolerance Convergence tolerance
     * @param max_iterations Maximum number of iterations
     * @return BCSResult containing the solution
     */
    BCSResult solve(double initial_gap = 1.0,
                    double tolerance = 1e-8,
                    int max_iterations = 1000);

    /**
     * Calculate the pairing gap equation residual
     * @param gap Current gap parameter
     * @param lambda Chemical potential
     * @return Gap equation residual
     */
    double gapEquation(double gap, double lambda) const;

    /**
     * Calculate the particle number equation residual
     * @param gap Current gap parameter
     * @param lambda Chemical potential
     * @return Particle number equation residual
     */
    double particleNumberEquation(double gap, double lambda) const;

    /**
     * Get the single-particle states
     */
    const std::vector<SingleParticleState>& getStates() const { return states_; }

    /**
     * Get the particle number
     */
    int getParticleNumber() const { return particle_number_; }

    /**
     * Get the pairing strength
     */
    double getPairingStrength() const { return pairing_strength_; }

    /**
     * Calculate quasiparticle energy
     * @param epsilon Single-particle energy
     * @param lambda Chemical potential
     * @param gap Pairing gap
     * @return Quasiparticle energy E_k
     */
    static double quasiparticleEnergy(double epsilon, double lambda, double gap);

    /**
     * Calculate BCS v^2 factor (occupation probability)
     * @param epsilon Single-particle energy
     * @param lambda Chemical potential
     * @param gap Pairing gap
     * @return v^2 factor
     */
    static double vSquared(double epsilon, double lambda, double gap);

private:
    std::vector<SingleParticleState> states_;
    int particle_number_;
    double pairing_strength_;

    /**
     * Find chemical potential for given gap using bisection
     */
    double findChemicalPotential(double gap, double tol = 1e-10) const;
};

/**
 * Utility functions for nuclear structure calculations
 */
namespace Utils {
    /**
     * Generate Woods-Saxon single-particle states
     * @param A Mass number
     * @param Z Proton number (for Coulomb correction)
     * @param is_proton Whether these are proton states
     * @param max_n Maximum principal quantum number
     * @return Vector of single-particle states
     */
    std::vector<SingleParticleState> generateWoodsSaxonStates(
        int A, int Z, bool is_proton, int max_n = 5);

    /**
     * Generate harmonic oscillator single-particle states
     * @param hw Oscillator parameter (in MeV)
     * @param max_shell Maximum major shell number
     * @return Vector of single-particle states
     */
    std::vector<SingleParticleState> generateHarmonicOscillatorStates(
        double hw, int max_shell);

    /**
     * Calculate the Fermi energy for a given set of states and particle number
     */
    double calculateFermiEnergy(const std::vector<SingleParticleState>& states,
                                 int particle_number);
}

} // namespace NuclearPairing

#endif // PAIR_APPROXIMATION_HPP
