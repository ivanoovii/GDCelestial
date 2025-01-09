#pragma once
#ifndef GDCELESTIAL_CELESTIAL_BODY_HPP
#define GDCELESTIAL_CELESTIAL_BODY_HPP

#include <cmath>

#include <godot_cpp/core/binder_common.hpp>
#include <godot_cpp/core/gdvirtual.gen.inc>

#include <godot_cpp/classes/node2d.hpp>
#include <godot_cpp/variant/vector2.hpp>
#include <godot_cpp/core/class_db.hpp>

#include "CelestialPhysics.hpp"

// Base class for both 2D and 3D celestial bodies.
// Contains shared math and parameters.
class CelestialBody
{
    /*
     * Dependancy graph for member variables.
     *
     * mu           -> specific angular momentum, specific mechanical energy
     * periapsis    -> specific angular momentum, specific mechanical energy, distance
     * eccentricity -> specific angular momentum, specific mechanical energy, distance
     * true anomlay -> distance
     *
     * specific angular momentum -> velocity
     * distance                  -> velocity, position
     * longitude of perigee      -> velocity, position
     *
     * mass -> childrens' mu
     *
     */

protected:
    const static godot::Vector2 reference_direction;

    bool enabled = true;

    double mass = 0.0;

    // Common orbital parameters.
    double mu = 0.0;           // Standard gravitational parameter of the parent celestial body.
    double periapsis = 0.0;    // Distance to the lowest point of the orbit.
    double eccentricity = 0.0; // Eccentricity of the orbit;
                               // 0 - circular orbit, [0; 1) - elliptical, 1 - parabolic, (1; +\infty) - hyperbolic.
    double longitude_of_perigee = 0.0;
    bool clockwise = true;

    double true_anomaly = 0.0;
    double mean_anomaly = 0.0;
    double distance = 0.0;

    // Integrals.
    double specific_angular_momentum = 0.0;
    double specific_mechanical_energy = 0.0;
    double mean_motion = 0.0;

    // Local cartesian coordinates.
    godot::Vector2 local_position = godot::Vector2(0.0, 0.0);
    godot::Vector2 local_velocity = godot::Vector2(0.0, 0.0);


    // On keplerian parameters changed callback.
    virtual void on_keplerian_parameters_changed() = 0;

    // Integrals calculation from keplerian parameters.
    double get_specific_angular_momentum_from_keplerian() const { return (clockwise ? 1.0 : -1.0) * std::sqrt(mu * get_semi_latus_rectum()); }
    double get_specific_mechanical_energy_from_keplerian() const { return -mu / (2.0 * get_semi_major_axis()); }
    double get_mean_motion_from_keplerian() const
    {
        //double abs_semi_major_axis = std::abs(get_semi_major_axis());
        //return std::sqrt(mu / abs_semi_major_axis) / abs_semi_major_axis;

        double multiplier = (eccentricity == 1.0) ? 1.0 : std::abs(1.0 - eccentricity * eccentricity);
        multiplier = std::sqrt(multiplier * multiplier * multiplier);
        return multiplier * mu * mu / (specific_angular_momentum * specific_angular_momentum * specific_angular_momentum);
    }

    // Conversion between local (on orbital plane) and global cartesian coordinates.
    virtual void local_cartesian_to_cartesian() = 0;
    virtual void cartesian_to_local_cartesian() = 0;

public:
    CelestialBody();
    ~CelestialBody();


    bool get_enabled() const { return enabled; }
    void set_enabled(bool new_enabled);


    double get_mass() const { return mass; }
    void set_mass(double new_mass);


    // Keplerian parameters.
    double get_mu() const { return mu; }
    void set_mu(double new_mu);

    double get_periapsis() const { return periapsis; }
    void set_periapsis(double new_periapsis);

    double get_eccentricity() const { return eccentricity; }
    void set_eccentricity(double new_eccentricity);

    double get_longitude_of_perigee() const { return longitude_of_perigee; }
    void set_longitude_of_perigee(double new_longitude_of_perigee);

    bool get_clockwise() const { return clockwise; }
    void set_clockwise(bool new_enabled);

    godot::Dictionary get_keplerian_parameters() const;


    // Supplementary orbital parameters.
    double get_semi_latus_rectum() const { return periapsis * (1.0 + eccentricity); }
    void set_semi_latus_rectum(double new_semi_latus_rectum) { set_periapsis(new_semi_latus_rectum / (1.0 + eccentricity)); };

    double get_semi_major_axis() const { return periapsis / (1.0 - eccentricity); }

    double get_apoapsis() const
    { return eccentricity < 1.0 ? periapsis * (1.0 + eccentricity) / (1.0 - eccentricity) : std::numeric_limits<double>::infinity(); }


    // Anomalies.
    double get_true_anomaly() const { return true_anomaly; }
    void set_true_anomaly(double new_true_anomaly);

    double get_mean_anomaly() const { return mean_anomaly; }
    void set_mean_anomaly(double new_mean_anomaly);

    double get_true_anomaly_at_distance(double _distance) const
    {
        // TODO: optimize: get_apoapsis and get_semi_latus_rectum can be combined.
        return ((_distance >= periapsis) && (_distance <= get_apoapsis())) ?
            std::acos((get_semi_latus_rectum() / _distance - 1.0) / eccentricity) :
            std::numeric_limits<double>::quiet_NaN();
    }


    // Get gravitational parameter used by children.
    virtual double get_children_mu() const
    { return CelestialPhysics::get_singleton()->gravitational_constant * mass; }

    virtual void update_children_mu() const = 0;


    // Spheres of influence for pathced conics or system generation.
    double get_influence_radius() const         { return periapsis * std::pow(get_children_mu() / mu, 0.4); }
    double get_influence_radius_squared() const { return periapsis * periapsis * std::pow(get_children_mu() / mu, 0.8); }
    double get_Hill_radius() const              { return periapsis * std::pow(get_children_mu() / (3.0 * mu), 1.0 / 3.0); }
    double get_Hill_radius_squared() const      { return periapsis * periapsis * std::pow(get_children_mu() / (3.0 * mu), 2.0 / 3.0); }


    // Conversion between coordinates.
    template<bool update_integrals = true>
    void keplerian_to_cartesian();

    void cartesian_to_keplerian();

    // Physics process step (orbtal motion).
    void celestial_physics_process(double delta);
};


// Class for 2D celestial body.
// TODO: use AnimatableBody2D as the base class when it is ready.
class CelestialBody2D : public godot::Node2D, public CelestialBody
{
    GDCLASS(CelestialBody2D, godot::Node2D);

protected:
    static void _bind_methods();

    godot::Vector2 linear_velocity = godot::Vector2(0.0, 0.0);
    double angular_velocity = 0.0;

    void on_keplerian_parameters_changed() override;

    void local_cartesian_to_cartesian() override;
    void cartesian_to_local_cartesian() override;

public:
    CelestialBody2D();
    ~CelestialBody2D();

    void update_children_mu() const override;

    godot::Vector2 get_linear_velocity() const { return linear_velocity; }
    void set_linear_velocity(godot::Vector2 new_linear_velocity);

    double get_angular_velocity() const { return angular_velocity; }
    void set_angular_velocity(double new_angular_velocity);

    void reparent_up();
    void reparent_down(CelestialBody2D* new_parent);
    void patch_conics();

    void _ready() override;
    void _physics_process(double delta) override;
};

#endif
