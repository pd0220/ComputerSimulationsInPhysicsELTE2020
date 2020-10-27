// used headers and libraries
#include <fstream>
#include <vector>
#include "DiffEqSolver.hh"

// -----------------------------------------------------------------------------------------------------------------

// constants and parameters
//
// error parameter
double const eps = 1e-5;
// pi
double const pi = 3.14159265358979323846;
// gravitatinal acceleration [m / s^2]
double const g = 9.80665;
// 6 * pi * eta * R / m = zeta ~ Stokes drag parameter [1 / s]
double zeta = 1;
// 0.5 * C * rho * A / m = zeta' ~ Newton drag parameter [1 / m]
// (can be used as density parameter for barometric simulations)
double zetaPrime = 1.2250;
// density of air at zero height above sea level [kg / m^3]
double const rho0 = 1.2250;
// molar mass of air [kg / mol]
double const M = 0.0289644;
// universal gas constant [N m / (mol K)]
double const R = 8.3144598;
// standard temperature at zero height above sea level [K]
double const T = 288.15;
// exponent coefficient in barometric formula [1 / m]
//double const lambda = M * g / R / T;
double const lambda = 1e-2;
// Newton's law of gravity pre-factor (gravitational constan * mass of Earth)
double const preFactorNewton = 6.67408e-11 * 5.9722e24;
// radius of Earth
double const radiusEarth = 6.371e6;

// -----------------------------------------------------------------------------------------------------------------

// equations of motion for the simulations (RHSs or derivative vectors)
//
// simple gravitational acceleration
auto GravityBasic = [&g = g](double, vector4<double> vec) -> vector4<double> {
    return {vec.v1, vec.v2, 0, -g};
};

// gravitational acceleration with Stokes drag
auto GravityStokes = [&g = g, &zeta = zeta](double, vector4<double> vec) -> vector4<double> {
    return {vec.v1, vec.v2, -zeta * vec.v1, -g - zeta * vec.v2};
};

// gravitational acceleration with Newton drag
auto GravityNewton = [&g = g, &zetaPrime = zetaPrime](double, vector4<double> vec) -> vector4<double> {
    double velMagnitude = std::sqrt(vec.v1 * vec.v1 + vec.v2 * vec.v2);
    return {vec.v1, vec.v2, -zetaPrime * velMagnitude * vec.v1, -g - zetaPrime * velMagnitude * vec.v2};
};

// gravitational acceleration with Newton drag including barometric formula for density of air (set *zetaPrime = 1* apart from density ~ multiply only with *barometric*)
auto GravityNewtonBarometric = [&lambda = lambda, &rho0 = rho0](double, vector4<double> vec) -> vector4<double> {
    double velMagnitude = std::sqrt(vec.v1 * vec.v1 + vec.v2 * vec.v2);
    double barometric = rho0 * std::exp(-lambda * vec.x2);
    return {vec.v1, vec.v2, -barometric * velMagnitude * vec.v1, -g - barometric * velMagnitude * vec.v2};
};

// Newton's law of gravity
auto GravityLaw = [&preFactorNewton = preFactorNewton, &radiusEarth = radiusEarth](double, vector4<double> vec) -> vector4<double> {
    double denominator = std::pow(vec.x1 * vec.x1 + vec.x2 * vec.x2, 1.5);
    return {vec.v1, vec.v2, -preFactorNewton / denominator * vec.x1, -preFactorNewton / denominator * vec.x2};
};

// -----------------------------------------------------------------------------------------------------------------

// further functions
//
// filename to save data
std::string fileName = "None.txt";

// callback functions ~ ho to write current data to file
// 1.
auto Callback1 = [&fileName = fileName](double t, vector4<double> vec, double h) {
    std::ofstream file;
    file.open(fileName, std::fstream::app);
    // print everything
    file << t << " " << vec.x1 << " " << vec.x2 << " " << vec.v1 << " " << vec.v2 << " " << h << std::endl;
    file.close();
};

// 2.
auto Callback2 = [&fileName = fileName](double, vector4<double> vec, double) {
    std::ofstream file;
    file.open(fileName, std::fstream::app);
    // only horizontal and vertical coordinates
    file << vec.x1 << " " << vec.x2 << std::endl;
    file.close();
};

// 3.
auto Callback3 = [&fileName = fileName](double, vector4<double>, double) {
    // do nothing
};

// 4.
auto Callback4 = [](double t, vector4<double> vec, double) {
    // write to terminal only
    std::cout << t << " " << vec.x1 << " " << vec.x2 << std::endl;
};

// break function ~ simulate only positive heigths
auto BreakLoop = [](vector4<double> vec) {
    // check if we are "above ground"
    if (vec.x2 < 0)
        return true;
    else
        return false;
};

// -----------------------------------------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // initial time
    double t0 = 0;
    // final time
    double t1 = 1000000;
    // step size
    double h = 0.01;

    // velocity magnitude
    double vMagnitude = 1;
    // file name to save data
    fileName = "Angles/Newton.txt";
    // number of datapoints
    int N = 200;
    // simulations
    for (int i = 0; i < N; i++)
    {
        // set zeta [0, 20]
        zetaPrime = double(i) / double(N) * 20;
        // run simulation with given initial angle and zeta
        std::ofstream file;
        // loop for angle values
        for (int j = 0; j < N; j++)
        {
            // set angle
            double phi = j * pi / 2 / N;
            // determine initial velocities for x and z directions
            double vx = vMagnitude * std::cos(phi);
            double vz = vMagnitude * std::sin(phi);
            // initial values
            vector4<double> init = {0., 0., vx, vz};
            // results
            vector4<double> yRes = CashKarpSolver(init, t0, t1, h, GravityNewton, Callback3, eps, BreakLoop);
            // write endpoint horizontal distance to file
            file.open(fileName, std::ofstream::app);
            file << yRes.x1 << " ";
            file.close();
        }
        file << std::endl;
    }

    //RK4Solver(init, t0, t1, h, GravityStokes, Callback1, BreakLoop);
}
