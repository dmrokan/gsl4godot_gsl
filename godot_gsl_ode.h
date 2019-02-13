#ifndef GODOT_GSL_ODE_H
#define GODOT_GSL_ODE_H

#include "godot_gsl_function.h"
#include "ode-initval2/gsl_odeiv2.h"

#define PARAM_COUNT_STEP 5
#define ODE_DELTA_T 1e-3

int __function(double t, const double y[], double f[], void *params);
int __jacobian(double t, const double y[], double f[], void *params);

class GodotGSLODE
{
public:
    GodotGSLODE() { }
    GodotGSLODE(size_t dim) { init(dim); };
    ~GodotGSLODE();
    void init(size_t dim);
    void set_function(GodotGSLFunction *fnc);
    void set_jacobian(GodotGSLFunction *fnc);
    void run(double t0, double t1, double _dt);
    int step(double dt);
    void func_execute() { function->execute(); }
    void jac_execute() { jacobian->execute(); }
    void set_node_path(Object *obj, const String subpath, const int index);
    void set_node_path_as_input(Object *obj, const String subpath, const String vn, const int index);
    void set_node_path_as_output(Object *obj, const String subpath, const String vn, const int index);
    void update_node_properties();
    void set_initial_conditions(const Array &x0_arr, const double t0);
    GodotGSLMatrix *x = NULL;
    GodotGSLMatrix *xdot = NULL;
    size_t dimension = 0;

private:
    GodotGSLFunction *function = NULL;
    GodotGSLFunction *jacobian = NULL;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver *driver;
    double t;
    double dt;
    double property_refresh_rate;
    unsigned int t_property = 0;
};

#endif
