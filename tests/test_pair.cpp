#include "PairApproximation.hpp"
#include <iostream>
#include <cmath>
#include <cassert>
#include <string>

using namespace NuclearPairing;

// Simple test framework
int tests_passed = 0;
int tests_failed = 0;

void test_result(bool condition, const std::string& test_name) {
    if (condition) {
        std::cout << "[PASS] " << test_name << "\n";
        tests_passed++;
    } else {
        std::cout << "[FAIL] " << test_name << "\n";
        tests_failed++;
    }
}

// Test quasiparticle energy calculation
void test_quasiparticle_energy() {
    // E_k = sqrt((e - lambda)^2 + Delta^2)
    double epsilon = 5.0;
    double lambda = 4.0;
    double gap = 3.0;
    
    double expected = std::sqrt(1.0 + 9.0);  // sqrt(10)
    double result = PairApproximation::quasiparticleEnergy(epsilon, lambda, gap);
    
    test_result(std::abs(result - expected) < 1e-10, 
                "Quasiparticle energy calculation");
}

// Test v^2 factor calculation
void test_v_squared() {
    // v^2 = 0.5 * (1 - (e - lambda) / E_k)
    double epsilon = 5.0;
    double lambda = 4.0;
    double gap = 3.0;
    
    double E_k = std::sqrt(1.0 + 9.0);
    double expected = 0.5 * (1.0 - 1.0 / E_k);
    double result = PairApproximation::vSquared(epsilon, lambda, gap);
    
    test_result(std::abs(result - expected) < 1e-10,
                "v^2 factor calculation");
}

// Test v^2 at Fermi surface (epsilon = lambda)
void test_v_squared_at_fermi() {
    double epsilon = 5.0;
    double lambda = 5.0;  // At Fermi surface
    double gap = 2.0;
    
    // At Fermi surface, v^2 = 0.5
    double result = PairApproximation::vSquared(epsilon, lambda, gap);
    
    test_result(std::abs(result - 0.5) < 1e-10,
                "v^2 = 0.5 at Fermi surface");
}

// Test v^2 limits for zero gap
void test_v_squared_zero_gap() {
    double gap = 0.0;
    
    // Below Fermi: v^2 = 1
    double v2_below = PairApproximation::vSquared(3.0, 5.0, gap);
    test_result(std::abs(v2_below - 1.0) < 1e-10,
                "v^2 = 1 below Fermi (zero gap)");
    
    // Above Fermi: v^2 = 0
    double v2_above = PairApproximation::vSquared(7.0, 5.0, gap);
    test_result(std::abs(v2_above) < 1e-10,
                "v^2 = 0 above Fermi (zero gap)");
}

// Test harmonic oscillator state generation
void test_ho_states() {
    double hw = 10.0;
    int max_shell = 2;
    
    auto states = Utils::generateHarmonicOscillatorStates(hw, max_shell);
    
    // Count total degeneracy
    int total_deg = 0;
    for (const auto& state : states) {
        total_deg += static_cast<int>(state.degeneracy);
    }
    
    // For shells 0, 1, 2: should have specific degeneracies
    // Shell 0: 1s1/2 (deg=2)
    // Shell 1: 1p3/2 (deg=4), 1p1/2 (deg=2)
    // Shell 2: 1d5/2 (deg=6), 1d3/2 (deg=4), 2s1/2 (deg=2)
    // Total: 2 + 4 + 2 + 6 + 4 + 2 = 20
    test_result(total_deg == 20,
                "Harmonic oscillator states total degeneracy");
    
    // Check states are sorted by energy
    bool sorted = true;
    for (size_t i = 1; i < states.size(); ++i) {
        if (states[i].energy < states[i-1].energy) {
            sorted = false;
            break;
        }
    }
    test_result(sorted, "HO states sorted by energy");
}

// Test Fermi energy calculation
void test_fermi_energy() {
    std::vector<SingleParticleState> states;
    states.emplace_back(1.0, 2.0);  // Can hold 2 particles
    states.emplace_back(2.0, 4.0);  // Can hold 4 particles
    states.emplace_back(3.0, 6.0);  // Can hold 6 particles
    
    // For 4 particles: fill first level (2), partially fill second
    double fermi = Utils::calculateFermiEnergy(states, 4);
    test_result(std::abs(fermi - 2.0) < 1e-10,
                "Fermi energy for 4 particles");
    
    // For 6 particles: fill first two levels
    fermi = Utils::calculateFermiEnergy(states, 6);
    test_result(std::abs(fermi - 2.0) < 1e-10,
                "Fermi energy for 6 particles");
}

// Test BCS solver with simple case
void test_bcs_solver_simple() {
    // Create simple model space
    std::vector<SingleParticleState> states;
    states.emplace_back(-5.0, 2.0);  // Deep level
    states.emplace_back(0.0, 4.0);   // Near Fermi
    states.emplace_back(5.0, 4.0);   // Above Fermi
    states.emplace_back(10.0, 2.0);  // High level
    
    int particle_number = 6;  // Fill first two levels in normal state
    double G = 0.5;  // Strong pairing
    
    PairApproximation solver(states, particle_number, G);
    BCSResult result = solver.solve(1.0, 1e-8, 1000);
    
    test_result(result.converged, "BCS solver converges");
    test_result(result.gap > 0.0, "BCS gap is positive");
    
    // Check particle number conservation
    double total_particles = 0.0;
    for (size_t i = 0; i < states.size(); ++i) {
        total_particles += states[i].degeneracy * result.occupation[i];
    }
    test_result(std::abs(total_particles - particle_number) < 0.1,
                "Particle number conserved");
}

// Test BCS with weak pairing (should give small or zero gap)
void test_bcs_weak_pairing() {
    std::vector<SingleParticleState> states;
    states.emplace_back(0.0, 4.0);
    states.emplace_back(10.0, 4.0);  // Large gap between levels
    
    int particle_number = 4;
    double G = 0.01;  // Very weak pairing
    
    PairApproximation solver(states, particle_number, G);
    BCSResult result = solver.solve(0.5, 1e-8, 1000);
    
    test_result(result.converged, "Weak pairing converges");
    // Weak pairing should give small gap
    test_result(result.gap < 1.0, "Weak pairing gives small gap");
}

// Test error handling
void test_error_handling() {
    std::vector<SingleParticleState> states;
    states.emplace_back(1.0, 4.0);
    
    bool caught_empty = false;
    try {
        std::vector<SingleParticleState> empty;
        PairApproximation solver(empty, 2, 0.5);
    } catch (const std::invalid_argument&) {
        caught_empty = true;
    }
    test_result(caught_empty, "Empty states throws exception");
    
    bool caught_odd = false;
    try {
        PairApproximation solver(states, 3, 0.5);  // Odd number
    } catch (const std::invalid_argument&) {
        caught_odd = true;
    }
    test_result(caught_odd, "Odd particle number throws exception");
    
    bool caught_negative_g = false;
    try {
        PairApproximation solver(states, 2, -0.5);  // Negative G
    } catch (const std::invalid_argument&) {
        caught_negative_g = true;
    }
    test_result(caught_negative_g, "Negative pairing strength throws exception");
}

// Test realistic nucleus (Sn-120)
void test_sn120() {
    double hw = 41.0 / std::pow(120.0, 1.0/3.0);  // ~8.3 MeV
    int max_shell = 5;
    
    auto states = Utils::generateHarmonicOscillatorStates(hw, max_shell);
    
    int N = 70;  // Neutron number for Sn-120
    // Use smaller pairing strength for this simplified model
    // Typical G ~ 25/A MeV ~ 0.2 MeV, but our model space is limited
    double G = 0.05;  // Adjusted for simplified harmonic oscillator model
    
    PairApproximation solver(states, N, G);
    BCSResult result = solver.solve(1.0, 1e-8, 1000);
    
    test_result(result.converged, "Sn-120 neutron calculation converges");
    
    // Gap should be positive and finite in this model
    test_result(result.gap > 0.0 && result.gap < 10.0,
                "Sn-120 gap is positive and finite");
    
    // Pairing energy should be negative
    test_result(result.pairing_energy < 0.0,
                "Pairing energy is negative");
}

int main() {
    std::cout << "=================================\n"
              << "Nuclear Pair Approximation Tests\n"
              << "=================================\n\n";
    
    test_quasiparticle_energy();
    test_v_squared();
    test_v_squared_at_fermi();
    test_v_squared_zero_gap();
    test_ho_states();
    test_fermi_energy();
    test_bcs_solver_simple();
    test_bcs_weak_pairing();
    test_error_handling();
    test_sn120();
    
    std::cout << "\n=================================\n"
              << "Results: " << tests_passed << " passed, " 
              << tests_failed << " failed\n"
              << "=================================\n";
    
    return tests_failed > 0 ? 1 : 0;
}
