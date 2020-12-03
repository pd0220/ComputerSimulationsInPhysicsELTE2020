// used headers and libraries
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

// ----------------------------------------------------------------------------------

// constant
// number of Monte Carlo sweeps...
// ... for autocorrelation times estimation
int const N_MC_autoCorr = 1000;
// ... for simulation
int const N_MC_simulation = (int)(1e6);

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

// calculate nth moment of position
auto MomentNth = [](std::vector<double> const &vec, int const &order, int const &Nt) {
    double tmp = 0.;
    for (double i : vec)
        tmp += std::pow(i, order);
    return tmp / (double)Nt;
};

// ----------------------------------------------------------------------------------

// calculate two-point correlation function with given time difference
auto TwoPointCorr = [](std::vector<double> const &vec, int const &Nt, int const &deltaTau) {
    double tmp = 0.;
    // loop for data points
    for (int iOuter = 0; iOuter < Nt; iOuter++)
    {
        // index for time translated datapoint (with periodic boundary conditions)
        int iInner = (iOuter + deltaTau) % Nt;
        tmp += vec[iOuter] * vec[iInner];
    }
    // return result
    return tmp / (double)(Nt - 1);
};

// ----------------------------------------------------------------------------------

// harmonic oscillator Euclidean action contribution from one coordinate
auto HarmOsc_S = [](double const &latticeSpacing, double const &mass, double const &omega, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = sq(omega * pathNow);

    // return action contribution
    return latticeSpacing * 0.5 * mass * (kin + pot);
};

// ----------------------------------------------------------------------------------

// anharmonic oscillator Euclidean action contribution from one coordinate
auto AnHarmOsc_S = [](double const &latticeSpacing, double const &mass, double const &mu, double const &lambda, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = -0.5 * mu * sq(pathNow) + lambda * sq(sq(pathNow)) + 1;

    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// Morse potential Euclidean action contribution from one coordinate
auto Morse_S = [](double const &latticeSpacing, double const &mass, double const &lambda, double const &x, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double diff = pathNow - x;
    double pot = sq(lambda) * (std::exp(-2 * diff) - 2 * std::exp(-diff));

    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// Pöschl-Teller potential Euclidean action contribution from one coordinate
auto PT_S = [](double const &latticeSpacing, double const &mass, double const &lambda, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = -lambda * (lambda + 1) / 2 / sq(std::cosh(pathNow));

    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// finite well potential Euclidean action contribution from one coordinate
auto Well_S = [](double const &latticeSpacing, double const &mass, double const &lambda, double const &a, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = 0.;
    if (std::abs(pathNow) <= a)
        pot = -lambda;
   
    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// Lennard-Jones potential Euclidean action contribution from one coordinate
auto LJ_S = [](double const &latticeSpacing, double const &mass, double const &epsilon, double const &sigma, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = 4 * epsilon * (std::pow(sigma / pathNow, 12) - std::pow(sigma / pathNow, 6));
   
    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// Coulomb potential Euclidean action contribution from one coordinate
auto Coulomb_S = [](double const &latticeSpacing, double const &mass, double const &alpha, double const &pathNow, double const &pathAfter) {
    // kinetic part
    double kin = 0.5 * mass * sq((pathAfter - pathNow) / latticeSpacing);
    // potential part
    double pot = - alpha / pathNow;
   
    // return action contribution
    return latticeSpacing * (kin + pot);
};

// ----------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // initial boundary for variations of paths
    double hInit = 0.5;
    // number of discretised time steps
    int Nt = 120;
    // harmonic oscillator parameters
    double mass = 1.;
    double omega = 1.;
    // lattice spacing
    double latticeSpacing = 1.;

    // random number generation with Mersenne Twister
    std::random_device rd{};
    std::mt19937 gen(rd());
    // uniform real distribution with fixed symmetric boundaries --> path variation & change acceptance for desired distribution ~ detailed balance, etc.
    std::uniform_real_distribution distrReal(0., 1.);
    // uniform integer distributin from 0 to (Nt - 1) --> site visiting
    std::uniform_int_distribution distrInt(0, Nt - 1);
    // random number generator lambdas
    auto randReal = [&distrReal, &gen]() { return (double)distrReal(gen); };
    auto randInt = [&distrInt, &gen]() { return (int)distrInt(gen); };

    // initial path (cold start)
    std::vector<double> pathInit(Nt);
    // fill initial path
    std::generate(pathInit.begin(), pathInit.end(), []() { return 0.; });

    //  ///////////\\\\\\\\\\\ 
    // ||| AUTOCORRELATIONS |||
    //  \\\\\\\\\\\///////////

    // declaring vector containers
    // which sites to visit in one MC sweep
    std::vector<int> visitSites(Nt);
    // acceptance for updates
    std::vector<double> acceptanceVec(Nt);
    // path to update iteratively
    std::vector<double> pathN = pathInit;
    // start loop for sweeps
    /*
    for (int iSweep = 0; iSweep < N_MC_autoCorr; iSweep++)
    {
        // determine which sites to visit in given MC sweep
        std::generate(visitSites.begin(), visitSites.end(), randInt);
        // generate random numbers for acceptance of site updates
        std::generate(acceptanceVec.begin(), acceptanceVec.end(), randReal);

        // loop for sites
        for (int iSite = 0; iSite < Nt; iSite++)
        {
            // which site to update in the time lattice
            int tau = visitSites[iSite];
            // time slices ~ before and after
            std::pair<int, int> tauBA = PeriodicBoundary(tau, Nt);
            //int tauBefore = tauBA.first;
            int tauAfter = tauBA.second;
            // possible new coordinate in path at the chosen site
            double tmpSite = pathN[tau] + hInit * (randReal() - 0.5);

            // calculate difference in Euclidean action (only one member of the summation) ~ S_new - S_old
            double sOld = HarmOsc_S(latticeSpacing, mass, omega, pathN[tau], pathN[tauAfter]);
            double sNew = HarmOsc_S(latticeSpacing, mass, omega, tmpSite, pathN[tauAfter]);
            double deltaS = sNew - sOld;

            // accept or reject change in site
            if (Rate(deltaS) > acceptanceVec[iSite] || Rate(deltaS) == 1)
                pathN[tau] = tmpSite;
        }
    }
    */

    // ideal step size determined from previous part via python analysis ~ hard coded here...
    // step size for fastest termalisation
    double h = 5;
    // estimated autocorrelation time for the corresponding step size (we will mulitply it with some bigger number... just in case)
    int const tauExp = 50;

    //  ///////////\\\\\\\\\\\ 
    // ||| SIMULATION START |||
    //  \\\\\\\\\\\///////////

    // setting lattice parameters
    Nt = 120;
    latticeSpacing = 1;

    // harmonic oscillator parameters
    mass = 1.;
    omega = 1.;

    // anharmonic oscillator parameters
    double const mu = 5.;
    double const lambda = 1.;

    // path to update iteratively
    pathN = pathInit;
    // separation number
    int sepTrigger = 1;

    // MEASUREMENT
    std::vector<double> meansMeasured;
    std::vector<double> varianceMeasured;
    std::vector<double> autoCorrMeasured;

    // generate initial path (cold start)
    pathInit.clear();
    pathInit.resize(Nt);
    std::generate(pathInit.begin(), pathInit.end(), []() { return 4.5; });

    pathN = pathInit;

    // start loop for sweeps
    for (int iSweep = 0; iSweep < N_MC_simulation; iSweep++)
    {
        // random number generation for site visiting
        std::uniform_int_distribution distrIntAdapt(0, Nt - 1);
        // random number generator lambdas
        auto randIntAdapt = [&distrIntAdapt, &gen]() { return (int)distrIntAdapt(gen); };

        // determine which sites to visit in given MC sweep
        visitSites.clear();
        visitSites.resize(Nt);
        std::generate(visitSites.begin(), visitSites.end(), randIntAdapt);
        // generate random numbers for acceptance of site updates
        acceptanceVec.clear();
        acceptanceVec.resize(Nt);
        std::generate(acceptanceVec.begin(), acceptanceVec.end(), randReal);

        // loop for sites
        for (int iSite = 0; iSite < Nt; iSite++)
        {
            // which site to update in the time lattice
            int tau = visitSites[iSite];
            // time slices ~ before and after
            std::pair<int, int> tauBA = PeriodicBoundary(tau, Nt);
            //int tauBefore = tauBA.first;
            int tauAfter = tauBA.second;
            // possible new coordinate in path at the chosen site
            double tmpSite = pathN[tau] + h * (randReal() - 0.5);

            // calculate difference in Euclidean action (only one member of the summation) ~ S_new - S_old
            // harmonic
            //double sOld = HarmOsc_S(latticeSpacing, mass, omega, pathN[tau], pathN[tauAfter]);
            //double sNew = HarmOsc_S(latticeSpacing, mass, omega, tmpSite, pathN[tauAfter]);
            // anharmonic
            //double sOld = AnHarmOsc_S(latticeSpacing, mass, mu, lambda, pathN[tau], pathN[tauAfter]);
            //double sNew = AnHarmOsc_S(latticeSpacing, mass, mu, lambda, tmpSite, pathN[tauAfter]);
            // Morse
            //double sOld = Morse_S(latticeSpacing, mass, lambda, 5, pathN[tau], pathN[tauAfter]);
            //double sNew = Morse_S(latticeSpacing, mass, lambda, 5, tmpSite, pathN[tauAfter]);
            // Pöschl-Teller
            //double sOld = PT_S(latticeSpacing, mass, lambda, pathN[tau], pathN[tauAfter]);
            //double sNew = PT_S(latticeSpacing, mass, lambda, tmpSite, pathN[tauAfter]);
            // finite potential well
            double sOld = Well_S(latticeSpacing, mass, lambda, 5, pathN[tau], pathN[tauAfter]);
            double sNew = Well_S(latticeSpacing, mass, lambda, 5, tmpSite, pathN[tauAfter]); 
            double deltaS = sNew - sOld;

            // accept or reject change in site
            if (Rate(deltaS) > acceptanceVec[iSite] || Rate(deltaS) == 1)
                pathN[tau] = tmpSite;
        }

        // MEASUREMENTS TRIGGERED
        // to save or not to save...
        sepTrigger++;
        if (sepTrigger % tauExp == 0)
        {
            // calculate mean
            meansMeasured.push_back(MomentNth(pathN, 1, Nt));
            // calculate second moment
            varianceMeasured.push_back(MomentNth(pathN, 2, Nt));
            // path
            for (double i : pathN)
                std::cout << i << " ";
            std::cout << std::endl;
        }
    }
}