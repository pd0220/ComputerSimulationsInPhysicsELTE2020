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
int const N_MC_simulation = 1000000;

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
auto Mean = [](std::vector<double> const &vec, int const &Nt) {
    double tmp = 0.;
    for (double i : vec)
        tmp += i;
    return tmp / Nt;
};

// ----------------------------------------------------------------------------------

// calculate variance
auto Variance = [](std::vector<double> const &vec, int const &Nt) {
    double tmp = 0.;
    for (double i : vec)
        tmp += sq(i);
    return tmp / (Nt - 1);
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
    return tmp / (Nt - 1);
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
// argv[1] --> h
int main(int, char **argv)
{
    // initial boundary for variations of paths
    double hInit = std::atof(argv[1]);
    // number of discretised time steps
    int Nt = 120;
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
    std::generate(pathInit.begin(), pathInit.end(), []() { return 10.; });

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
    for (int iSweep = 0; iSweep < N_MC_autoCorr; iSweep++)
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
        double mean = Mean(pathN, Nt);
        // calculate variance
        double var = Variance(pathN, Nt);

        // calculate two-point correlation function to estimate autocorrelation times
        double tmpCorr = 0.;
        for (int tau = 0; tau < Nt; tau++)
        {
            tmpCorr += TwoPointCorr(pathN, Nt, tau);
        }
        //std::cout << tmpCorr / Nt << std::endl;
    }

    //  ///////////\\\\\\\\\\\ 
    // ||| SIMULATION START |||
    //  \\\\\\\\\\\///////////

    // ideal step size determined from previous part via python analysis ~ hard coded here...
    // step size for fastest termalisation
    int const hStep = 25;
    // determined autocorrelation time for the corresponding step size (we will mulitply it with 100... just in case)
    int const tauExp = 200;

    // MEASUREMENT
    std::vector<double> meansMeasured;
    std::vector<double> varianceMeasured;
    std::vector<double> autoCorrMeasured;

    // path to update iteratively
    pathN = pathInit;
    // separation number
    int sepTrigger = 1;
    // start loop for sweeps
    for (int iSweep = 0; iSweep < N_MC_simulation; iSweep++)
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

        // MEASUREMENTS
        // to save or not to save...
        sepTrigger++;
        if (sepTrigger % tauExp == 0)
        {
            // calculate mean
            meansMeasured.push_back(Mean(pathN, Nt));
            // calculate variance
            varianceMeasured.push_back(Variance(pathN, Nt));

            // calculate two-point correlation function to estimate autocorrelation times
            double tmpCorr = 0.;
            for (int tau = 0; tau < Nt; tau++)
            {
                tmpCorr += TwoPointCorr(pathN, Nt, tau);
            }

            // averaged autocorrelation function over time separation
            autoCorrMeasured.push_back(tmpCorr / Nt);

            // jackknife analysis in outer python notebook
        }
    }

    // write results to screen
    for (int i = 0; i < static_cast<int>(meansMeasured.size()); i++)
        std::cout << meansMeasured[i] << " " << varianceMeasured[i] << " " << autoCorrMeasured[i] << std::endl;
}