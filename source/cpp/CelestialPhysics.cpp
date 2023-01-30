#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include "CelestialPhysics.hpp"


CelestialPhysics *CelestialPhysics::singleton = nullptr;

void CelestialPhysics::_bind_methods()
{
    godot::ClassDB::bind_method(godot::D_METHOD("get_time_scale"), &CelestialPhysics::get_time_scale);
    godot::ClassDB::bind_method(godot::D_METHOD("set_time_scale", "new_time_scale"), &CelestialPhysics::set_time_scale);
	ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "time_scale", godot::PROPERTY_HINT_NONE), "set_time_scale", "get_time_scale");

    godot::ClassDB::bind_method(godot::D_METHOD("true_anomaly_to_eccentric_anomaly", "true_anomaly", "eccentricity"), &CelestialPhysics::true_anomaly_to_eccentric_anomaly);
    godot::ClassDB::bind_method(godot::D_METHOD("eccentric_anomaly_to_true_anomaly", "eccentric_anomaly", "eccentricity"), &CelestialPhysics::eccentric_anomaly_to_true_anomaly);
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


double CelestialPhysics::true_anomaly_to_eccentric_anomaly(double true_anomaly, double eccentricity) const
{
    double beta = eccentricity / (1.0 + std::sqrt(1.0 - eccentricity * eccentricity));
    return true_anomaly - 2.0 * std::atan(beta * std::sin(true_anomaly) / (1.0 + beta * std::cos(true_anomaly)));
}

double CelestialPhysics::eccentric_anomaly_to_true_anomaly(double eccentric_anomaly, double eccentricity) const
{
    double beta = eccentricity / (1.0 + std::sqrt(1.0 - eccentricity * eccentricity));
    return eccentric_anomaly + 2.0 * std::atan(beta * std::sin(eccentric_anomaly) / (1.0 - beta * std::cos(eccentric_anomaly)));
}
