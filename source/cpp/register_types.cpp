#include "register_types.h"

#include <gdextension_interface.h>

#include <godot_cpp/core/class_db.hpp>
#include <godot_cpp/core/defs.hpp>
#include <godot_cpp/classes/engine.hpp>
#include <godot_cpp/godot.hpp>

#include <iostream>
#include <fstream>

#include "CelestialBody.hpp"
#include "CelestialPhysics.hpp"


static CelestialPhysics *_celestial_physics;

void gdcelestial_initialize(godot::ModuleInitializationLevel p_level)
{
    if (p_level == godot::MODULE_INITIALIZATION_LEVEL_SCENE)
    {
        godot::ClassDB::register_class<CelestialBody2D>();
        godot::ClassDB::register_class<CelestialPhysics>();

        _celestial_physics = memnew(CelestialPhysics);
        godot::Engine::get_singleton()->register_singleton("CelestialPhysics", CelestialPhysics::get_singleton());
    }
}

void gdcelestial_terminate(godot::ModuleInitializationLevel p_level)
{
    if (p_level == godot::MODULE_INITIALIZATION_LEVEL_SCENE)
    {
        godot::Engine::get_singleton()->unregister_singleton("CelestialPhysics");
        memdelete(_celestial_physics);
    }
}

extern "C"
{
    GDExtensionBool GDE_EXPORT gdextension_init(GDExtensionInterfaceGetProcAddress p_get_proc_address, GDExtensionClassLibraryPtr p_library, GDExtensionInitialization *r_initialization)
    //GDExtensionBool GDE_EXPORT gdextension_init(const GDExtensionInterface *p_interface, GDExtensionClassLibraryPtr p_library, GDExtensionInitialization *r_initialization)
    {
        godot::GDExtensionBinding::InitObject init_obj(p_get_proc_address, p_library, r_initialization);
        //godot::GDExtensionBinding::InitObject init_obj(p_interface, p_library, r_initialization);

        init_obj.register_initializer(gdcelestial_initialize);
        init_obj.register_terminator(gdcelestial_terminate);
        init_obj.set_minimum_library_initialization_level(godot::MODULE_INITIALIZATION_LEVEL_SCENE);

        return init_obj.init();
    }
}
