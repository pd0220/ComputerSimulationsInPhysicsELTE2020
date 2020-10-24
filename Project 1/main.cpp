// used headers and libraries
#include <fstream>
#include "DiffEqSolver.hh"

// -----------------------------------------------------------------------------------------------------------------

// constants and parameters
//
// error parameter
double const eps = 1e-5;
// gravitatinal acceleration [m / s^2]
double const g = 9.80665;
// 6 * pi * eta * R / m = zeta ~ Stokes drag parameter [1 / s]
double const zeta = 10;

// -----------------------------------------------------------------------------------------------------------------

// equations of motion
//
// simple gravitational acceleration
auto GravityBasic = [&g = g](double, vector4<double> vec) -> vector4<double> {
    return {vec.v1, vec.v2, 0, -g};
};

// gravitational acceleration with Stokes drag
auto GravityStokes = [&g = g, &zeta = zeta](double, vector4<double> vec) -> vector4<double> {
    return {vec.v1, vec.v2, -zeta * vec.v1, -g - zeta * vec.v2};
};

// -----------------------------------------------------------------------------------------------------------------

// further functions
//
// filename to save data
std::string const fileName = "Stokes_RK3.txt";

// callback function ~ write current data to file
auto Callback = [&fileName = fileName](double t, vector4<double> vec, double h) {
    std::ofstream file;
    file.open(fileName, std::fstream::app);
    file << t << " " << vec.x1 << " " << vec.x2 << " " << vec.v1 << " " << vec.v2 << " " << h << std::endl;
    file.close();
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
    vector4<double> init{0, 0, 10, 10};
    // initial time
    double t0 = 0;
    // final time
    double t1 = 10;
    // step size
    double h = 0.001;

    // integration methods
    //CashKarpSolver(init, t0, t1, h, GravityStokes, Callback, eps, BreakLoop);
    RK4Solver(init, t0, t1, h, GravityStokes, Callback, BreakLoop);
}
