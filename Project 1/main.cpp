// used headers and libraries
#include <fstream>
#include <vector>
#include "DiffEqSolver.hh"

// -----------------------------------------------------------------------------------------------------------------

// constants and parameters
//
// error parameter
double const eps = 1e-5;
// gravitatinal acceleration [m / s^2]
double const g = 9.80665;
// 6 * pi * eta * R / m = zeta ~ Stokes drag parameter [1 / s]
double const zeta = 1;
// 0.5 * C * rho * A / m = zeta' ~ Newton drag parameter [1 / m]
// (can be used as density parameter for barometric simulations)
double const zetaPrime = 1.2250;
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

// -----------------------------------------------------------------------------------------------------------------

// further functions
//
// filename to save data
std::string fileName = "Barometric/Newton1.txt";

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
    // initial values
    vector4<double> initContainer = {0., 0., 150, 150.};

    // file names
    std::string fileNameNewton = "BarometricTest/NewtonDists.txt";
    std::string fileNameNewtonBar = "BarometricTest/NewtonBarDists.txt";

    // initial time
    double t0 = 0;
    // final time
    double t1 = 100;
    // step size
    double h = 0.01;

    // simulations
    for (int i = 1; i < 19; i++)
    {
        vector4<double> init{0., 0., i * 50., i * 50.};
        std::ofstream file;
        file.open(fileNameNewton, std::ofstream::app);
        vector4<double> yRes = CashKarpSolver(init, t0, t1, h, GravityNewton, Callback3, eps, BreakLoop);
        file << yRes.x1 << " " << yRes.x2 << std::endl;
        file.close();
    }
    for (int i = 1; i < 19; i++)
    {
        vector4<double> init{0., 0., i * 50., i * 50.};
        std::ofstream file;
        file.open(fileNameNewtonBar, std::ofstream::app);
        vector4<double> yRes = CashKarpSolver(init, t0, t1, h, GravityNewtonBarometric, Callback3, eps, BreakLoop);
        file << yRes.x1 << " " << yRes.x2 << std::endl;
        file.close();
    }
    
    //RK4Solver(init, t0, t1, h, GravityStokes, Callback1, BreakLoop);
}
