/* register_types.cpp */

#include "register_types.h"
#include <core/class_db.h>
#include "godot_gsl.h"

void register_godot_gsl_types() {
    ClassDB::register_class<GodotGSL>();
}

void unregister_godot_gsl_types() {
    //nothing to do here
}
