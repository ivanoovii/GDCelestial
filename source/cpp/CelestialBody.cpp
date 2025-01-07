#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

#include <godot_cpp/classes/sprite2d.hpp>

#include "CelestialBody.hpp"


//////// CelestialBody ////////
// private:

// protected:
const godot::Vector2 CelestialBody::reference_direction = godot::Vector2(1.0, 0.0);

// public:
CelestialBody::CelestialBody() { }
CelestialBody::~CelestialBody() { }

void CelestialBody::set_enabled(bool new_enabled)
{
    enabled = new_enabled;
}


void CelestialBody::set_mass(double new_mass)
{
    mass = std::max(new_mass, 0.0);
}


void CelestialBody::set_mu(double new_mu)
{
    mu = std::max(new_mu, 0.0);
    keplerian_to_cartesian();
    on_keplerian_parameters_changed();
}

void CelestialBody::set_periapsis(double new_periapsis)
{
    periapsis = std::max(new_periapsis, 0.0);
    keplerian_to_cartesian();
    on_keplerian_parameters_changed();
}

void CelestialBody::set_eccentricity(double new_eccentricity)
{
    eccentricity = std::max(new_eccentricity, 0.0);

    // Reset anomalies, as they are invalidated when eccentricity changes.
    true_anomaly = 0.0;
    mean_anomaly = 0.0;

    keplerian_to_cartesian();
    on_keplerian_parameters_changed();
}

void CelestialBody::set_longitude_of_perigee(double new_longitude_of_perigee)
{
    longitude_of_perigee = std::remainder(new_longitude_of_perigee, 2.0 * PI);
    keplerian_to_cartesian<false>(); // Do not update the integrals, they are the same.
    on_keplerian_parameters_changed();
}

void CelestialBody::set_clockwise(bool new_clockwise)
{
    // Change revolvation direction by changing the sign of the specific angular momentum.
    if (clockwise != new_clockwise)
    { specific_angular_momentum = -specific_angular_momentum; }

    clockwise = new_clockwise;
    keplerian_to_cartesian<false>(); // Just need to update the sign of the specific angular momentum, no need for a full update.
    on_keplerian_parameters_changed();
}

godot::Dictionary CelestialBody::get_keplerian_parameters() const
{
    godot::Dictionary result;

    result["mu"] = mu;
    result["periapsis"] = periapsis;
    result["eccentricity"] = eccentricity;
    result["longitude_of_perigee"] = longitude_of_perigee;

    return result;
}


void CelestialBody::set_true_anomaly(double new_true_anomaly)
{
    if (eccentricity < 1.0) { new_true_anomaly = std::remainder(new_true_anomaly, 2.0 * PI); }
    else
    {
        double max_true_anomaly_absolute_value = std::acos(-1.0 / eccentricity) -
            static_cast<double>(std::numeric_limits<float>::epsilon());
        // `float` machine epsilon for stability.
        new_true_anomaly = clamp(new_true_anomaly, -max_true_anomaly_absolute_value, max_true_anomaly_absolute_value);
    }
    true_anomaly = new_true_anomaly;

    mean_anomaly = CelestialPhysics::true_anomaly_to_mean_anomaly(new_true_anomaly, eccentricity);
    keplerian_to_cartesian<false>(); // Do not update the integrals, they are the same.
}


void CelestialBody::set_mean_anomaly(double new_mean_anomaly)
{
    if (eccentricity < 1.0) { new_mean_anomaly = std::remainder(new_mean_anomaly, 2.0 * PI); }
    mean_anomaly = new_mean_anomaly;

    true_anomaly = CelestialPhysics::mean_anomaly_to_true_anomaly(new_mean_anomaly, eccentricity, true_anomaly);
    keplerian_to_cartesian<false>(); // Do not update the integrals, they are the same.
}


template<bool update_integrals>
void CelestialBody::keplerian_to_cartesian()
{
    if constexpr(update_integrals)
    {
        // Updating the integrals.
        specific_angular_momentum  = get_specific_angular_momentum_from_keplerian();
        specific_mechanical_energy = get_specific_mechanical_energy_from_keplerian();
        mean_motion = get_mean_motion_from_keplerian();
    }

    // Updating the distance.
    distance = get_semi_latus_rectum() / (1.0 + eccentricity * std::cos(true_anomaly));

    // Updating local cartesian coordinates.
    godot::Vector2 direction = reference_direction.rotated(longitude_of_perigee + true_anomaly);
    local_position = distance * direction;
    local_velocity = specific_angular_momentum * (
            std::sin(true_anomaly) * eccentricity * direction / get_semi_latus_rectum() +
            direction.rotated(PI / 2.0) / distance
    );

    // Updating global cartesian coordinates.
    // (must be implemented in derived classes)
    local_cartesian_to_cartesian();
}

void CelestialBody::cartesian_to_keplerian()
{
    // Updating global cartesian coordinates.
    // (must be implemented in derived classes)
    cartesian_to_local_cartesian();

    distance = local_position.length();
    double velocity_quad = local_velocity.length_squared();

    // Eccentricity and longitude of perigee.
    godot::Vector2 eccentricity_vector = ((velocity_quad - mu / distance) * local_position - local_position.dot(local_velocity) * local_velocity) / mu;
    eccentricity = eccentricity_vector.length();
    longitude_of_perigee = reference_direction.angle_to(eccentricity_vector);
    longitude_of_perigee = std::remainder(longitude_of_perigee, 2.0 * PI);

    // Integrals.
    specific_angular_momentum  = local_position.cross(local_velocity);
    specific_mechanical_energy = velocity_quad / 2.0 - mu / distance;
    mean_motion = get_mean_motion_from_keplerian(); // TODO direct calculation from cartesian.

    clockwise = (specific_angular_momentum >= 0.0);

    // Periapsis height.
    periapsis = specific_angular_momentum * specific_angular_momentum / (mu * (1.0 + eccentricity));

    // True and mean anomaly.
    true_anomaly = reference_direction.angle_to(local_position) - longitude_of_perigee;
    true_anomaly = std::remainder(true_anomaly, 2.0 * PI);
    mean_anomaly = CelestialPhysics::true_anomaly_to_mean_anomaly(true_anomaly, eccentricity);

    on_keplerian_parameters_changed();
}



void CelestialBody::celestial_physics_process(double delta)
{
    if (enabled && !godot::Engine::get_singleton()->is_editor_hint())
    {
        double physics_dt = delta * CelestialPhysics::get_singleton()->time_scale;
        set_mean_anomaly(mean_anomaly + physics_dt * mean_motion);
    }
}



//////// CelestialBody2D ////////
// private:

// protected:
void CelestialBody2D::_bind_methods()
{
    godot::ClassDB::bind_method<godot::MethodDefinition, bool (CelestialBody2D::*)() const>(godot::D_METHOD("get_enabled"), &CelestialBody2D::get_enabled);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(bool)>(godot::D_METHOD("set_enabled", "new_enabled"), &CelestialBody2D::set_enabled);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::BOOL, "enabled"), "set_enabled", "get_enabled");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_mass"), &CelestialBody2D::get_mass);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_mass", "new_mass"), &CelestialBody2D::set_mass);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "mass"), "set_mass", "get_mass");


    ADD_GROUP("Keplerian parameters", "keplerian_");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_mu"), &CelestialBody2D::get_mu);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_mu", "new_mu"), &CelestialBody2D::set_mu);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_mu", godot::PROPERTY_HINT_RANGE, "0,99999,0.00001,or_greater,hide_slider"), "set_mu", "get_mu");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_periapsis"), &CelestialBody2D::get_periapsis);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_periapsis", "new_periapsis"), &CelestialBody2D::set_periapsis);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_periapsis", godot::PROPERTY_HINT_RANGE, "0,99999,0.001,or_greater,hide_slider,suffix:px"), "set_periapsis", "get_periapsis");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_eccentricity"), &CelestialBody2D::get_eccentricity);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_eccentricity", "new_eccentricity"), &CelestialBody2D::set_eccentricity);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_eccentricity", godot::PROPERTY_HINT_RANGE, "0,99999,0.001,or_greater,hide_slider"), "set_eccentricity", "get_eccentricity");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_longitude_of_perigee"), &CelestialBody2D::get_longitude_of_perigee);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_longitude_of_perigee", "new_longitude_of_perigee"),
            &CelestialBody2D::set_longitude_of_perigee);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_longitude_of_perigee", godot::PROPERTY_HINT_RANGE,
                "-180,180,0.1,or_less,or_greater,radians"), "set_longitude_of_perigee", "get_longitude_of_perigee");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_true_anomaly"), &CelestialBody2D::get_true_anomaly);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_true_anomaly", "new_true_anomaly"), &CelestialBody2D::set_true_anomaly);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_true_anomaly", godot::PROPERTY_HINT_RANGE, "-180,180,0.1,or_less,or_greater,radians"), "set_true_anomaly", "get_true_anomaly");

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_mean_anomaly"), &CelestialBody2D::get_mean_anomaly);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(double)>(godot::D_METHOD("set_mean_anomaly", "new_mean_anomaly"), &CelestialBody2D::set_mean_anomaly);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "keplerian_mean_anomaly", godot::PROPERTY_HINT_RANGE, "-180,180,0.1,or_less,or_greater,radians"), "set_mean_anomaly", "get_mean_anomaly");

    godot::ClassDB::bind_method<godot::MethodDefinition, bool (CelestialBody2D::*)() const>(godot::D_METHOD("get_clockwise"), &CelestialBody2D::get_clockwise);
    godot::ClassDB::bind_method<godot::MethodDefinition, void (CelestialBody2D::*)(bool)>(godot::D_METHOD("set_clockwise", "new_clockwise"), &CelestialBody2D::set_clockwise);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::BOOL, "keplerian_clockwise"), "set_clockwise", "get_clockwise");


    ADD_SIGNAL(godot::MethodInfo("keplerian_parameters_changed", godot::PropertyInfo(godot::Variant::DICTIONARY, "new_parameters")));


    ADD_GROUP("Cartesian parameters", "cartesian_");

    godot::ClassDB::bind_method(godot::D_METHOD("get_linear_velocity"), &CelestialBody2D::get_linear_velocity);
    godot::ClassDB::bind_method(godot::D_METHOD("set_linear_velocity", "new_linear_velocity"), &CelestialBody2D::set_linear_velocity);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::VECTOR2, "cartesian_linear_velocity", godot::PROPERTY_HINT_NONE, "suffix:px/s"), "set_linear_velocity", "get_linear_velocity");

    godot::ClassDB::bind_method(godot::D_METHOD("get_angular_velocity"), &CelestialBody2D::get_angular_velocity);
    godot::ClassDB::bind_method(godot::D_METHOD("set_angular_velocity", "new_angular_velocity"), &CelestialBody2D::set_angular_velocity);
    ADD_PROPERTY(godot::PropertyInfo(godot::Variant::FLOAT, "cartesian_angular_velocity", godot::PROPERTY_HINT_NONE, U"radians,suffix:\u00B0/s"), "set_angular_velocity", "get_angular_velocity");


    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_semi_latus_rectum"), &CelestialBody2D::get_semi_latus_rectum);
    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_semi_major_axis"), &CelestialBody2D::get_semi_major_axis);

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_children_mu"), &CelestialBody2D::get_children_mu);

    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_sphere_of_influence"), &CelestialBody2D::get_sphere_of_influence);
    godot::ClassDB::bind_method<godot::MethodDefinition, double (CelestialBody2D::*)() const>(godot::D_METHOD("get_Hill_sphere"), &CelestialBody2D::get_Hill_sphere);
}


void CelestialBody2D::on_keplerian_parameters_changed()
{
    emit_signal("keplerian_parameters_changed", get_keplerian_parameters());
}

void CelestialBody2D::local_cartesian_to_cartesian()
{
    set_position(local_position);
    linear_velocity = local_velocity;
}

void CelestialBody2D::cartesian_to_local_cartesian()
{
    local_position = get_position();
    local_velocity = linear_velocity;
}


// public:
CelestialBody2D::CelestialBody2D()
{ }

CelestialBody2D::~CelestialBody2D()
{ }


void CelestialBody2D::update_children_mu() const
{
    double children_mu = get_children_mu();

    godot::TypedArray<Node> children = get_children();
    for (int64_t index = 0; index < children.size(); ++index)
    {
        CelestialBody2D *child_as_celestial_body_2d = godot::Node::cast_to<CelestialBody2D>(children[index]);

        if (child_as_celestial_body_2d)
        {
            child_as_celestial_body_2d->set_mu(children_mu);
        }
    }
}


void CelestialBody2D::set_linear_velocity(godot::Vector2 new_linear_velocity)
{
    linear_velocity = new_linear_velocity;
    if (is_inside_tree())
    {
        cartesian_to_keplerian();
    }
}

void CelestialBody2D::set_angular_velocity(double new_angular_velocity)
{
    angular_velocity = new_angular_velocity;
}


void CelestialBody2D::_ready()
{
    on_keplerian_parameters_changed();
    update_children_mu();
    add_to_group("celestial_bodies");
}

void CelestialBody2D::_physics_process(double delta)
{
    celestial_physics_process(delta);
}
