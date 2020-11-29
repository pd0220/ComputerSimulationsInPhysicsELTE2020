// used headers and libraries
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// ----------------------------------------------------------------------------------

// constant
// number of Monte Carlo sweeps
int const N_MC = 100000;

// ----------------------------------------------------------------------------------

// lambda to calculate square
auto sq = [](double const &x) { return x * x; };

// ----------------------------------------------------------------------------------

// calculate propability of changing site
auto Rate = [](double const &deltaS) {
    if (deltaS < 0)
        return 1.;
    else
        return std::exp(-deltaS);
};

// ----------------------------------------------------------------------------------

// periodic boundary conditions in time
auto PeriodicBoundary = [](int const &tau, int const &Nt) {
    // taking periodic boundary conditions into account for time --> **MAGIC**
    int tauBefore = (tau + Nt - 1) % Nt;
    int tauAfter = (tau + 1) % Nt;

    // return time slices
    return std::pair<int, int>(tauBefore, tauAfter);
};

// ----------------------------------------------------------------------------------

// calculate mean
auto Mean = [](std::vector<double> const &vec) {
    double tmp = 0;
    for (double i : vec)
        tmp += i;
    return tmp / static_cast<int>(vec.size());
};

// ----------------------------------------------------------------------------------

// harmonic oscillator Euclidean action contribution from one coordinate
auto HarmOsc_S = [](double const &mass, double const &frequency, double const &pathBefore, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = sq(pathAfter - pathBefore);
    // potential part
    double pot = sq(frequency * pathNow);

    // return action contribution
    return 0.5 * mass * (kin + pot);
};

// ----------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // initial boundary for variations of paths
    double hInit = 0.5;
    // number of discretised time steps
    int Nt = 30;
    // mass
    double mass = 1.;
    // frequency
    double frequency = 1.;

    // random number generation with Mersenne Twister
    std::random_device rd{};
    std::mt19937 gen(rd());
    // uniform real distribution with fixed symmetric boundaries --> path variation
    std::uniform_real_distribution distrReal(-hInit, hInit);
    // uniform integer distributin from 0 to (Nt-1) --> site visiting
    std::uniform_int_distribution distrInt(0, Nt - 1);
    // uniform real distribution from 0 to 1 --> change acceptance for desired distribution ~ detailed balance, etc.
    std::uniform_real_distribution distr01(0., 1.);
    // random number generator lambdas
    auto randReal = [&distrReal, &gen]() { return (double)distrReal(gen); };
    auto randInt = [&distrInt, &gen]() { return (int)distrInt(gen); };
    auto rand_01 = [&distr01, &gen]() { return (double)distr01(gen); };

    // initial path (hot start)
    std::vector<double> pathInit(Nt);
    // fill initial path
    std::generate(pathInit.begin(), pathInit.end(), []() { return 0.; });

    //  ///////////\\\\\\\\\\\ 
    // ||| SIMULATION START |||
    //  \\\\\\\\\\\///////////

    // declaring vector containers
    // which sites to visit in one MC sweep
    std::vector<int> visitSites(Nt);
    // acceptance for updates
    std::vector<double> acceptanceVec(Nt);
    // path to update iteratively
    std::vector<double> pathN = pathInit;
    // start loop for sweeps
    for (int iSweep = 0; iSweep < N_MC; iSweep++)
    {
        // determine which sites to visit in given MC sweep
        std::generate(visitSites.begin(), visitSites.end(), randInt);
        // generate random numbers for acceptance of site updates
        std::generate(acceptanceVec.begin(), acceptanceVec.end(), rand_01);

        // loop for sites
        for (int iSite = 0; iSite < Nt; iSite++)
        {
            // which site to update in the time lattice
            int tau = visitSites[iSite];
            // time slices ~ before and after
            std::pair<int, int> tauBA = PeriodicBoundary(tau, Nt);
            int tauBefore = tauBA.first, tauAfter = tauBA.second;
            // possible new coordinate in path at the chosen site
            double tmpSite = pathN[tau] + randReal();

            // calculate difference in Euclidean action (only one member of the summation) ~ S_new - S_old
            double sOld = HarmOsc_S(mass, frequency, pathN[tauBefore], pathN[tau], pathN[tauAfter]);
            double sNew = HarmOsc_S(mass, frequency, pathN[tauBefore], tmpSite, pathN[tauAfter]);
            double deltaS = sNew - sOld;

            // accept or reject change in site
            if (Rate(deltaS) > acceptanceVec[iSite])
                pathN[tau] = tmpSite;
        }

        // calculate mean
        double mean = Mean(pathN);

        for (double i : pathN)
            std::cout << i << " ";
        std::cout << std::endl;
    }
}
