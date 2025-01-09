#include <cmath>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include "CelestialPhysics.hpp"


CelestialPhysics *CelestialPhysics::singleton = nullptr;

void CelestialPhysics::_bind_methods()
{
    godot::ClassDB::bind_method(godot::D_METHOD("get_time_scale"), &CelestialPhysics::get_time_scale);
    godot::ClassDB::bind_method(godot::D_METHOD("set_time_scale", "new_time_scale"), &CelestialPhysics::set_time_scale);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "time_scale", godot::PROPERTY_HINT_NONE), "set_time_scale", "get_time_scale");

    godot::ClassDB::bind_static_method("CelestialPhysics", godot::D_METHOD("true_anomaly_to_eccentric_anomaly", "true_anomaly", "eccentricity"), &CelestialPhysics::true_anomaly_to_eccentric_anomaly);
    godot::ClassDB::bind_static_method("CelestialPhysics", godot::D_METHOD("eccentric_anomaly_to_true_anomaly", "eccentric_anomaly", "eccentricity"), &CelestialPhysics::eccentric_anomaly_to_true_anomaly);
}


void CelestialPhysics::set_time_scale(double new_time_scale)
{
    time_scale = new_time_scale;
}


CelestialPhysics::CelestialPhysics()
{
    ERR_FAIL_COND(singleton != nullptr);
    singleton = this;
}

CelestialPhysics::~CelestialPhysics()
{
    ERR_FAIL_COND(singleton != this);
    singleton = nullptr;
}


// Conics math.
double CelestialPhysics::true_anomaly_to_eccentric_anomaly(double true_anomaly, double eccentricity)
{
    double beta = eccentricity / (1.0 + std::sqrt(1.0 - eccentricity * eccentricity));
    return true_anomaly - 2.0 * std::atan(beta * std::sin(true_anomaly) / (1.0 + beta * std::cos(true_anomaly)));
}

double CelestialPhysics::eccentric_anomaly_to_true_anomaly(double eccentric_anomaly, double eccentricity)
{
    double beta = eccentricity / (1.0 + std::sqrt(1.0 - eccentricity * eccentricity));
    return eccentric_anomaly + 2.0 * std::atan(beta * std::sin(eccentric_anomaly) / (1.0 - beta * std::cos(eccentric_anomaly)));
}


double CelestialPhysics::mean_anomaly_to_true_anomaly(double mean_anomaly, double eccentricity,
        double true_anomaly_hint, double tolerance, unsigned int max_iter)
{
    //double eccentric_anomaly = mean_anomaly_to_eccentric_anomaly(mean_anomaly, eccentricity);
    //return eccentric_anomaly_to_true_anomaly(eccentric_anomaly, eccentricity);

    // Using Newton method to convert mean anomaly to true anomaly.
    int orbit_type = (eccentricity >= 1.0) + (eccentricity > 1.0); // 0 - elliptical, 1 - parabolical, 2 - hyperbolical.
    switch (orbit_type)
    {
        default:
        case 0:
        {
            mean_anomaly = std::remainder(mean_anomaly, 2.0 * PI);
            double eccentric_anomaly = 2.0 * std::atan(std::tan(0.5 * true_anomaly_hint) * std::sqrt((1.0 - eccentricity) / (1.0 + eccentricity)));

            double region_min = mean_anomaly - eccentricity;
            double region_max = mean_anomaly + eccentricity;

            for (unsigned int iter = 0; iter < max_iter; ++iter) // Newton iteration for Kepler equation;
            {
                eccentric_anomaly = clamp(eccentric_anomaly, region_min, region_max);

                double residual = eccentric_anomaly - eccentricity * std::sin(eccentric_anomaly) - mean_anomaly;
                double derivative = 1.0 - eccentricity * std::cos(eccentric_anomaly);

                double delta = -residual / derivative;
                eccentric_anomaly += delta;
                if (std::abs(delta) < tolerance) { break; }

                if (iter + 1 == max_iter)
                { godot::UtilityFunctions::print("Mean anomaly to true anomaly conversion failed: the solver did not converge."); }
            }

            return 2.0 * std::atan(std::tan(0.5 * eccentric_anomaly) * std::sqrt((1.0 + eccentricity) / (1.0 - eccentricity)));
        }
        case 1:
        {
            double z = std::cbrt(3.0 * mean_anomaly + std::sqrt(1 + 9.0 * mean_anomaly * mean_anomaly));
            return 2.0 * std::atan(z - 1.0 / z);
        }
        case 2:
        {
            double eccentric_anomaly = 2.0 * std::atanh(std::tan(0.5 * true_anomaly_hint) * std::sqrt((eccentricity - 1.0) / (eccentricity + 1.0)));
            for (unsigned int iter = 0; iter < max_iter; ++iter) // Newton iteration for Kepler equation;
            {
                double residual = eccentricity * std::sinh(eccentric_anomaly) - eccentric_anomaly - mean_anomaly;
                double derivative = eccentricity * std::cosh(eccentric_anomaly) - 1.0;

                double delta = -residual / derivative;
                eccentric_anomaly += delta;
                if (std::abs(delta) < tolerance) { break; }

                if (iter + 1 == max_iter)
                { godot::UtilityFunctions::print("Mean anomaly to true anomaly conversion failed: the solver did not converge."); }
            }

            return 2.0 * std::atan(std::tanh(0.5 * eccentric_anomaly) * std::sqrt((eccentricity + 1.0) / (eccentricity - 1.0)));
        }
    }
}

double CelestialPhysics::true_anomaly_to_mean_anomaly(double true_anomaly, double eccentricity)
{
    //double eccentric_anomaly = true_anomaly_to_eccentric_anomaly(true_anomaly, eccentricity);
    //return eccentric_anomaly_to_mean_anomaly(eccentric_anomaly, eccentricity);

    int orbit_type = (eccentricity >= 1.0) + (eccentricity > 1.0); // 0 - elliptic, 1 - parabolic, 2 - hyperbolic.
    switch (orbit_type)
    {
        default:
        case 0:
        {
            double eccentric_anomaly = 2.0 * std::atan(std::tan(0.5 * true_anomaly) * std::sqrt((1.0 - eccentricity) / (1.0 + eccentricity)));
            return eccentric_anomaly - eccentricity * std::sin(eccentric_anomaly);
        }
        case 1:
        {
            double hta_tan = std::tan(0.5 * true_anomaly);
            return 0.5 * hta_tan * (1.0 + hta_tan * hta_tan / 3.0);
        }
        case 2:
        {
            double eccentric_anomaly = 2.0 * std::atanh(std::tan(0.5 * true_anomaly) * std::sqrt((eccentricity - 1.0) / (eccentricity + 1.0)));
            return eccentricity * std::sinh(eccentric_anomaly) - eccentric_anomaly;
        }
    }
}

