#include "PairApproximation.hpp"
#include <iostream>
#include <iomanip>
#include <string>

using namespace NuclearPairing;

void printUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n"
              << "\nNuclear Pair Approximation Calculator\n"
              << "Performs BCS pairing calculations for nuclear structure.\n"
              << "\nOptions:\n"
              << "  -A <number>     Mass number (default: 120)\n"
              << "  -Z <number>     Proton number (default: 50)\n"
              << "  -N <number>     Neutron number (default: A-Z)\n"
              << "  -G <value>      Pairing strength in MeV (default: 0.25)\n"
              << "  -hw <value>     Oscillator parameter in MeV (default: auto)\n"
              << "  -shells <num>   Maximum shell number (default: 5)\n"
              << "  -proton         Calculate for protons (default)\n"
              << "  -neutron        Calculate for neutrons\n"
              << "  -v              Verbose output\n"
              << "  -h              Show this help message\n"
              << "\nExample:\n"
              << "  " << program << " -A 120 -Z 50 -G 0.25 -neutron\n";
}

int main(int argc, char* argv[]) {
    // Default parameters for Sn-120 (tin isotope)
    int A = 120;        // Mass number
    int Z = 50;         // Proton number
    int N = -1;         // Neutron number (auto-calculated if -1)
    double G = 0.25;    // Pairing strength (MeV)
    double hw = -1;     // Oscillator parameter (auto if < 0)
    int max_shells = 5; // Maximum major shell
    bool is_proton = false;  // Calculate for neutrons by default
    bool verbose = false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "-A" && i + 1 < argc) {
            A = std::stoi(argv[++i]);
        } else if (arg == "-Z" && i + 1 < argc) {
            Z = std::stoi(argv[++i]);
        } else if (arg == "-N" && i + 1 < argc) {
            N = std::stoi(argv[++i]);
        } else if (arg == "-G" && i + 1 < argc) {
            G = std::stod(argv[++i]);
        } else if (arg == "-hw" && i + 1 < argc) {
            hw = std::stod(argv[++i]);
        } else if (arg == "-shells" && i + 1 < argc) {
            max_shells = std::stoi(argv[++i]);
        } else if (arg == "-proton") {
            is_proton = true;
        } else if (arg == "-neutron") {
            is_proton = false;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            std::cerr << "Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Calculate neutron number if not specified
    if (N < 0) {
        N = A - Z;
    }

    // Calculate oscillator parameter if not specified
    if (hw < 0) {
        hw = 41.0 / std::pow(A, 1.0/3.0);  // Standard formula
    }

    // Determine particle number
    int particle_number = is_proton ? Z : N;

    // Ensure even particle number
    if (particle_number % 2 != 0) {
        std::cerr << "Warning: Particle number " << particle_number 
                  << " is odd. Adjusting to " << (particle_number - 1) 
                  << " for BCS calculation.\n";
        particle_number -= 1;
    }

    std::cout << "======================================\n"
              << " Nuclear Pair Approximation (BCS)\n"
              << "======================================\n\n"
              << "Nucleus: A=" << A << ", Z=" << Z << ", N=" << N << "\n"
              << "Calculating for: " << (is_proton ? "Protons" : "Neutrons") << "\n"
              << "Particle number: " << particle_number << "\n"
              << "Pairing strength G: " << G << " MeV\n"
              << "Oscillator hw: " << std::fixed << std::setprecision(3) << hw << " MeV\n"
              << "Max shells: " << max_shells << "\n\n";

    // Generate single-particle states
    auto states = Utils::generateHarmonicOscillatorStates(hw, max_shells);

    if (verbose) {
        std::cout << "Single-particle states:\n"
                  << std::setw(5) << "n" 
                  << std::setw(5) << "l" 
                  << std::setw(7) << "j"
                  << std::setw(10) << "Energy"
                  << std::setw(8) << "Deg"
                  << "\n";
        std::cout << std::string(35, '-') << "\n";
        for (const auto& state : states) {
            std::cout << std::setw(5) << state.quantum_n
                      << std::setw(5) << state.quantum_l
                      << std::setw(7) << std::setprecision(1) << state.quantum_j
                      << std::setw(10) << std::setprecision(3) << state.energy
                      << std::setw(8) << std::setprecision(0) << state.degeneracy
                      << "\n";
        }
        std::cout << "\n";
    }

    // Verify we have enough states
    int total_capacity = 0;
    for (const auto& state : states) {
        total_capacity += static_cast<int>(state.degeneracy);
    }

    if (total_capacity < particle_number) {
        std::cerr << "Error: Not enough single-particle states. "
                  << "Capacity: " << total_capacity 
                  << ", Required: " << particle_number << "\n"
                  << "Try increasing the number of shells with -shells option.\n";
        return 1;
    }

    // Create pair approximation solver
    try {
        PairApproximation solver(states, particle_number, G);

        // Solve the BCS equations
        double initial_gap = 1.0;  // Initial guess in MeV
        BCSResult result = solver.solve(initial_gap);

        // Output results
        std::cout << "BCS Solution:\n"
                  << std::string(35, '-') << "\n"
                  << std::fixed << std::setprecision(6)
                  << "Pairing gap (Delta): " << result.gap << " MeV\n"
                  << "Chemical potential (lambda): " << result.chemical_potential << " MeV\n"
                  << "Pairing energy: " << result.pairing_energy << " MeV\n"
                  << "Converged: " << (result.converged ? "Yes" : "No") << "\n"
                  << "Iterations: " << result.iterations << "\n\n";

        if (verbose && result.converged) {
            std::cout << "Occupation probabilities (v^2):\n"
                      << std::setw(5) << "n"
                      << std::setw(5) << "l"
                      << std::setw(7) << "j"
                      << std::setw(10) << "Energy"
                      << std::setw(10) << "v^2"
                      << std::setw(10) << "u"
                      << std::setw(10) << "v"
                      << "\n";
            std::cout << std::string(57, '-') << "\n";
            
            for (size_t i = 0; i < states.size(); ++i) {
                std::cout << std::setw(5) << states[i].quantum_n
                          << std::setw(5) << states[i].quantum_l
                          << std::setw(7) << std::setprecision(1) << states[i].quantum_j
                          << std::setw(10) << std::setprecision(3) << states[i].energy
                          << std::setw(10) << std::setprecision(5) << result.occupation[i]
                          << std::setw(10) << result.u_factors[i]
                          << std::setw(10) << result.v_factors[i]
                          << "\n";
            }
        }

        // Summary
        double fermi = Utils::calculateFermiEnergy(states, particle_number);
        std::cout << "\nSummary:\n"
                  << "Fermi energy (no pairing): " << std::setprecision(3) << fermi << " MeV\n";
        
        if (result.gap > 1e-6) {
            std::cout << "Pairing is present with gap = " << result.gap << " MeV\n";
        } else {
            std::cout << "No pairing (normal state)\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
