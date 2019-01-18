/* godot_gsl.h */
#ifndef GODOT_GSL_H
#define GODOT_GSL_H

#include "core/reference.h"
#include "core/os/memory.h"
#include "core/map.h"
#include <gsl_vector_float.h>
#include <gsl_matrix_float.h>
#include <gsl_blas.h>

#define ERR_DIMENSION_MISMATCH -1
#define ERR_NULL_VALUE -2

#define MATRIX_VECTOR_PROD 1
#define MATRIX_MATRIX_PROD 2
#define VECTOR_VECTOR_PROD 3

typedef enum val_type {
    GSL_SCALER,
    GSL_VECTOR,
    GSL_MATRIX,
} val_type;

#define FLOAT_BLOCK_SIZE (sizeof(gsl_block_float))
#define FLOAT_MATRIX_SIZE (sizeof(gsl_matrix_float))
#define FLOAT_MATRIX_BLOCK_ADDR_SKIP (3 * sizeof(size_t) + sizeof(float*) \
                                      + sizeof(gsl_block*) + sizeof(int))
#define FLOAT_MATRIX_DATA_ADDR_SKIP (FLOAT_MATRIX_BLOCK_ADDR_SKIP + sizeof(gsl_block))
#define FLOAT_VECTOR_SIZE (sizeof(gsl_matrix_float))
#define FLOAT_VECTOR_DATAPTR_SKIP (3 * sizeof(size_t))
#define FLOAT_VECTOR_BLOCKPTR_SKIP (FLOAT_VECTOR_DATAPTR_SKIP + sizeof(float*))
#define FLOAT_VECTOR_DATA_SKIP (FLOAT_VECTOR_BLOCKPTR_SKIP + sizeof(float*) + sizeof(int))

/*
 * Module class
 */

class GodotGSLMatrix;

class GodotGSL : public Reference {
    GDCLASS(GodotGSL, Reference);

protected:
    static void _bind_methods();

public:
    GodotGSL();
    void matrix_new(String vn, const size_t row_count, const size_t col_count);
    void matrix_new_from_array(String vn, const Array a);
    void add(String vn1, String vn2);
    Array to_array(String vn);
    void set_zero(String vn);
    void prod(String vn1, String vn2, String to);
    void vector_new(String vn, const size_t size, bool is_row_vector);
    void vector_new_from_array(String vn, const Array a, bool is_row_vector);

private:
    void _add_variable(String vn, GodotGSLMatrix* mtx);
    Map<String, GodotGSLMatrix*> variables;
};

/*
**************** Matrix class
*/

#define DTYPE gsl_matrix_float
#define VTYPE gsl_vector_float
#define STYPE float
class GodotGSLMatrix {
public:
    GodotGSLMatrix(const size_t row_count, const size_t col_count) { init(row_count, col_count); }
    GodotGSLMatrix(const Array a) { init(a); }
    GodotGSLMatrix() { }
    ~GodotGSLMatrix();
    void init(const size_t row_count, const size_t col_count);
    void init(const Array a);
    void set_zero();
    void set_identity();
    int add(const GodotGSLMatrix &a);
    int sub(const GodotGSLMatrix &a);
    int mul_elements(const GodotGSLMatrix &a);
    int scale(const STYPE a);
    STYPE get(const size_t row, const size_t col);
    DTYPE *get_ptr() const;
    VTYPE *get_vec_ptr() const;
    void set(const size_t row, const size_t col, const STYPE a);
    GodotGSLMatrix *prod(const GodotGSLMatrix &a);

    size_t size[2];

private:
    DTYPE *gsl_mtx;
    VTYPE *gsl_vec;
    val_type type = GSL_MATRIX;
    bool is_row_vector = false;
    DTYPE *_alloc(const size_t row_count, const size_t col_count);
    void _vector_alloc(size_t s);
};
static GodotGSLMatrix *NULLMTX;

#undef DTYPE

#endif // GODOT_GSL_H
