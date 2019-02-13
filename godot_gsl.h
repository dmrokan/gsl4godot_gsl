/* godot_gsl.h */
#ifndef GODOT_GSL_H
#define GODOT_GSL_H

#include "core/reference.h"
#include "godot_gsl_matrix.h"
#include "godot_gsl_function.h"
#include "godot_gsl_ode.h"

class GodotGSL : public Reference {
    GDCLASS(GodotGSL, Reference);

protected:
    static void _bind_methods();

public:
    GodotGSL();
    ~GodotGSL();
    void matrix_new(String vn, const size_t row_count, const size_t col_count);
    void matrix_new_from_array(String vn, const Array a);
    void add(String vn1, String vn2);
    Array to_array(String vn);
    void set_zero(String vn);
    void prod(String vn1, String vn2, String to);
    void vector_new(String vn, const size_t size, bool is_row_vector);
    void vector_new_from_array(String vn, const Array a, bool is_row_vector);
    void fx(const String fn, const String a, String to);
    void function(const String fn, const Array args);
    void instruction(const String in, const Array args);
    void callf(const String fn);
    void ode(const String on, const size_t dim);
    void ode_set_fx(const String on, const String fn);
    void ode_run(const String on, const Array x0, const Array time_interval, const double dt);
    void ode_set_node_path(const String on, Variant obj_var, const String subpath, const int index);
    void ode_set_init_cond(const String on, const Array x0, const double t0);
    void ode_run_delta(const String on, const double delta);
    void matrix_set_identity(const String vn);
    void ode_set_node_path_as_input(const String on, Variant obj_var, const String subpath, const String vn, const int index);

private:
    void _add_variable(String vn, GodotGSLMatrix* mtx);
    void _add_function(String fn, GodotGSLFunction* fnc);
    Map<String, GodotGSLMatrix*> variables;
    Map<String, GodotGSLFunction*> functions;
    Map<String, GodotGSLODE*> odes;
    GodotGSLFunction* current;
};

#endif // GODOT_GSL_H
