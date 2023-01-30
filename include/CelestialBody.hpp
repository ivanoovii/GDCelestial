#pragma once
#ifndef GDCELESTIAL_CELESTIAL_BODY_HPP
#define GDCELESTIAL_CELESTIAL_BODY_HPP

#include <cmath>

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
    double mu = 0.0;           // Standart gravitational parameter of parent celestial body.
    double periapsis = 0.0;    // Distance to the lowest point of the orbit.
    double eccentricity = 0.0; // Eccentricity of an orbit;
                               // 0 - circular orbit, [0; 1) - elliptical, 1 - parabolic, (1; +\infty) - hyperbolic.
    double longitude_of_perigee = 0.0;
    bool clockwise = true;

    double true_anomaly = 0.0;
    double distance = 0.0;

    // Integrals.
    double specific_angular_momentum = 0.0;
    double specific_mechanical_energy = 0.0;

    // Local cartesian coordinates.
    godot::Vector2 local_position = godot::Vector2(0.0, 0.0);
    godot::Vector2 local_velocity = godot::Vector2(0.0, 0.0);


    virtual void on_keplerian_parameters_changed() = 0;

    double get_specific_angular_momentum_from_keplerian() const { return (clockwise ? 1.0 : -1.0) * std::sqrt(mu * get_semi_latus_rectum()); }
    double get_specific_mechanical_energy_from_keplerian() const { return -mu / (2.0 * get_semi_major_axis()); }

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


    double get_semi_latus_rectum() const { return periapsis * (1.0 + eccentricity); }
    void set_semi_latus_rectum(double new_semi_latus_rectum) { set_periapsis(new_semi_latus_rectum / (1.0 + eccentricity)); };

    double get_semi_major_axis() const { return periapsis / (1.0 - eccentricity); }


    double get_true_anomaly() const { return true_anomaly; }
    void set_true_anomaly(double new_true_anomaly);


    // Get gravitational parameter used by children.
    virtual double get_children_mu() const
    { return CelestialPhysics::get_singleton()->gravitational_constant * mass; }

    virtual void update_children_mu() const = 0;


    // Spheres of influence for pathced conics or system generation.
    double get_sphere_of_influence() const { return periapsis * std::pow(get_children_mu() / mu, 0.4); }
    double get_Hill_sphere() const { return periapsis * std::pow(get_children_mu() / (3.0 * mu), 1.0 / 3.0); }


    // Conversion between coordinates.
    template<bool update_integrals = true>
    void keplerian_to_cartesian();

    void cartesian_to_keplerian();

    // Physics process step (orbtal motion).
    void celestial_physics_process(double delta);
};


// Class for 2D celestial body.
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

    void _ready() override;
    void _physics_process(double delta) override;
};

#endif
