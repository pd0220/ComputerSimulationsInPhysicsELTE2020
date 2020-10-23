// used headers and libraries
#include "DiffEqSolver.hh"

// -----------------------------------------------------------------------------------------------------------------

// main function
int main(int, char **)
{
    // delta error
    double eps = 1e-3;
    // parameter(s)
    double g = 9.81;
    // test equaion
    auto test = [&g](double, vector4<double> vec) -> vector4<double> {
        return {vec.v1, vec.v2, 0, -g};
    };

    // callback function
    auto Callback = [](double t, vector4<double> vec, double h) {
        std::cout << t << " " << vec.x1 << " " << vec.x2 << " " << vec.v1 << " " << vec.v2 << " " << h << std::endl;
    };

    // initial values
    vector4<double> init{0, 0, 10, 10};
    // initial time
    double t0 = 0;
    // final time
    double t1 = 10;
    // step size
    double h = 0.01;

    // integration
    CashKarpSolver(init, t0, t1, h, test, Callback, eps);
    //RK4Solver(init, t0, t1, h, test, Callback);
}
