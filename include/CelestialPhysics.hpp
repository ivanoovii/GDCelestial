#pragma once
#ifndef GDCELESTIAL_CELESTIAL_PHYSICS_HPP
#define GDCELESTIAL_CELESTIAL_PHYSICS_HPP


#include <godot_cpp/classes/object.hpp>
#include <godot_cpp/core/class_db.hpp>


class CelestialPhysics : public godot::Object
{
    GDCLASS(CelestialPhysics, godot::Object);

    static CelestialPhysics *singleton;

protected:
    static void _bind_methods();

    double get_time_scale() { return time_scale; }
    void set_time_scale(double new_time_scale);

public:
    static CelestialPhysics *get_singleton() { return singleton; }

    // Constants.
    // Default mass unit: 10^24 kg.
    // Default length unit: 10^6 m (1 pixel = 1000 km).
    // Default time unit: 1 s.
    double gravitational_constant = 6.6743015e-5;
    double time_scale = 3600.0; // Default: one real-time hour per second.

    CelestialPhysics();
    ~CelestialPhysics();

    // Conics math.
    static double true_anomaly_to_eccentric_anomaly(double true_anomaly, double eccentricity);
    static double eccentric_anomaly_to_true_anomaly(double eccentric_anomaly, double eccentricity);

    //static double mean_anomaly_to_eccentric_anomaly(double mean_anomaly, double eccentricity);
    //static double eccentric_anomaly_to_mean_anomaly(double eccentric_anomaly, double eccentricity);

    static double mean_anomaly_to_true_anomaly(double mean_anomaly, double eccentricity,
            double true_anomaly_hint = 0.0, double tolerance = 1e-9, unsigned int max_iters = 64);
    static double true_anomaly_to_mean_anomaly(double true_anomaly, double eccentricity);
};

#endif
