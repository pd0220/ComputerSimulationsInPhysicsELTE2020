// numerical differential equation solvers
// adaptive Cash-Karp and 4th order Runge-Kutta methods implemented

// used headers and libraries
#include <fstream>
#include <iostream>
#include "vector4D.hh"
// -----------------------------------------------------------------------------------------------------------------

// adaptive Cash-Karp method implementation for numerical integration
// via one RK5 and one embedded 4th order step
template <typename State, typename T, typename RHS, typename Callback>
auto CashKarpSolver(State const &y0, T const &t0, T const &t1, T h, RHS f, Callback cb, T const &eps)
{
    // setting initial values
    // initial time
    T t = t0;
    // initial state
    State y = y0;
    // adaptive steps until reaching integration bound
    while (t < t1)
    {
        // set minimal error parameter
        T Delta0 = length(eps * h * y);

        // helper states to calculate steps
        State k1, k2, k3, k4, k5, k6;
        // future error
        T Delta = (T)0.;

        // adaptive steps ~ do until error is small enough
        do
        {
            // fifth order step
            k1 = h * f(t, y);
            k2 = h * f(t + h * (T)(1. / 5.), y + (T)(1. / 5.) * k1);
            k3 = h * f(t + h * (T)(3. / 10.), y + (T)(3. / 40.) * k1 + (T)(9. / 40.) * k2);
            k4 = h * f(t + h * (T)(3. / 5.), y + (T)(3. / 10.) * k1 - (T)(9. / 10.) * k2 + (T)(6. / 5.) * k3);
            k5 = h * f(t + h, y - (T)(11. / 54.) * k1 + (T)(5. / 2.) * k2 - (T)(70. / 27.) * k3 + (T)(35. / 27.) * k4);
            k6 = h * f(t + h * (T)(7. / 8.), y + (T)(1631. / 55296.) * k1 + (T)(175. / 512.) * k2 + (T)(575. / 13824.) * k3 + (T)(44275. / 110592.) * k4 + (T)(253. / 4096.) * k5);

            // error estimation via fourth order step
            State Delta_state = ((T)(37. / 378.) - (T)(2825. / 27648.)) * k1 +
                                ((T)(250. / 621.) - (T)(18575. / 48384.)) * k3 +
                                ((T)(125. / 594.) - (T)(13525. / 55296.)) * k4 +
                                ((T)(-277. / 14336.)) * k5 +
                                ((T)(512. / 1771.) - (T)(1. / 4.)) * k6;

            // reduce error state via some nomarlizing function
            Delta = length(Delta_state);
            // update step size
            h = h * 0.3 * std::pow(std::abs(Delta0 / Delta), 0.2);

        } while (Delta0 < Delta);

        // update step size
        h = h * 0.3 * std::pow(std::abs(Delta0 / Delta), 0.25);

        // taking the step
        y = y + (T)(37. / 378.) * k1 + (T)(250. / 621.) * k3 + (T)(125. / 594.) * k4 + (T)(512. / 1771.) * k6;
        t = t + h;

        // callback function
        cb(t, y, h);
    }
    return y;
}

// 4th order Runge-Kutta method
template <typename State, typename T, typename RHS, typename Callback>
auto RK4Solver(State y0, T t0, T t1, T h, RHS f, Callback cb)
{

    // setting initial values
    // initial time
    T t = t0;
    // initial state
    State y = y0;
    // steps until reaching integration bound
    while (t < t1)
    {
        // handle overstep
        if (t + h > t1)
        {
            h = t1 - t;
        }
        // fourth order step
        State k1 = f(t, y);
        State k2 = f(t + h * (T)0.5, y + (h * (T)0.5) * k1);
        State k3 = f(t + h * (T)0.5, y + (h * (T)0.5) * k2);
        State k4 = f(t + h, y + h * k3);

        // taking the step
        y = y + (k1 + k4 + (T)2. * (k2 + k3)) * (h / (T)6.);
        t = t + h;

        // callback function
        cb(t, y, h);
    }
    return y;
}