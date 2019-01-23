/* godot_gsl_matrix.h */
#ifndef GODOT_GSL_MATRIX_H
#define GODOT_GSL_MATRIX_H

#include "core/reference.h"
#include "core/os/memory.h"
#include "core/map.h"
#include <gsl_vector_double.h>
#include <gsl_matrix_double.h>
#include <gsl_blas.h>
#include <math.h>
#include "godot_gsl_macros.h"

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

typedef enum cond_type {
    EQUAL_SIZE,
    NONZERO_SIZE,
} cond_type;

/* Math function names */
#define FN_EQ "="
#define FN_ADD "+"
#define FN_SUB "-"
#define FN_MUL "*"
#define FN_SIN "sin"

/* Float size definitions */
#define FLOAT_BLOCK_SIZE (sizeof(gsl_block_float))
#define FLOAT_MATRIX_SIZE (sizeof(gsl_matrix_float))
#define FLOAT_MATRIX_BLOCK_ADDR_SKIP (3 * sizeof(size_t) + sizeof(float*) \
                                      + sizeof(gsl_block*) + sizeof(int))
#define FLOAT_MATRIX_DATA_ADDR_SKIP (FLOAT_MATRIX_BLOCK_ADDR_SKIP + sizeof(gsl_block))
#define FLOAT_VECTOR_SIZE (sizeof(gsl_matrix_float))
#define FLOAT_VECTOR_DATAPTR_SKIP (3 * sizeof(size_t))
#define FLOAT_VECTOR_BLOCKPTR_SKIP (FLOAT_VECTOR_DATAPTR_SKIP + sizeof(float*))
#define FLOAT_VECTOR_DATA_SKIP (FLOAT_VECTOR_BLOCKPTR_SKIP + sizeof(float*) + sizeof(int))

/* Double size definitions */
#define DOUBLE_BLOCK_SIZE (sizeof(gsl_block))
#define DOUBLE_MATRIX_SIZE (sizeof(gsl_matrix))
#define DOUBLE_MATRIX_BLOCK_ADDR_SKIP (3 * sizeof(size_t) + sizeof(double*) \
                                      + sizeof(gsl_block*) + sizeof(int))
#define DOUBLE_MATRIX_DATA_ADDR_SKIP (DOUBLE_MATRIX_BLOCK_ADDR_SKIP + sizeof(gsl_block))
#define DOUBLE_VECTOR_SIZE (sizeof(gsl_matrix))
#define DOUBLE_VECTOR_DATAPTR_SKIP (3 * sizeof(size_t))
#define DOUBLE_VECTOR_BLOCKPTR_SKIP (DOUBLE_VECTOR_DATAPTR_SKIP + sizeof(double*))
#define DOUBLE_VECTOR_DATA_SKIP (DOUBLE_VECTOR_BLOCKPTR_SKIP + sizeof(double*) + sizeof(int))

/* Base data type definitons */
#define GSL_BASE_DOUBLE

#ifdef GSL_BASE_FLOAT
#define DTYPE gsl_matrix_float
#define VTYPE gsl_vector_float
#define BTYPE gsl_block_float
#define STYPE float
#define STRIDE 1
#define BLOCK_SIZE FLOAT_BLOCK_SIZE
#define MATRIX_SIZE FLOAT_MATRIX_SIZE
#define MULTIPLICITY 1
#define ATOMIC_SIZE (sizeof(float))
#define MATRIX_BLOCK_ADDR_SKIP FLOAT_MATRIX_BLOCK_ADDR_SKIP
#define MATRIX_DATA_ADDR_SKIP FLOAT_MATRIX_DATA_ADDR_SKIP
#define MATRIX_SIZE FLOAT_MATRIX_SIZE
#define BLOCK_SIZE FLOAT_BLOCK_SIZE
/* Function macros */
#define gsl_matrix_set gsl_matrix_float_set
#define gsl_matrix_set_zero gsl_matrix_float_set_zero
#define gsl_matrix_set_identity gsl_matrix_float_set_identity
#define gsl_matrix_add gsl_matrix_float_add
#define gsl_matrix_get gsl_matrix_float_get
#define gsl_matrix_sub gsl_matrix_float_sub
#define gsl_matrix_mul_elements gsl_matrix_float_mul_elements
#define gsl_matrix_scale gsl_matrix_float_scale
#define gsl_blas_gemv gsl_blas_sgemv
#define gsl_blas_gemm gsl_blas_sgemm
#endif

#ifdef GSL_BASE_DOUBLE
#define DTYPE gsl_matrix
#define VTYPE gsl_vector
#define BTYPE gsl_block
#define STYPE double
#define STRIDE 1
#define BLOCK_SIZE DOUBLE_BLOCK_SIZE
#define MATRIX_SIZE DOUBLE_MATRIX_SIZE
#define MULTIPLICITY 1
#define ATOMIC_SIZE (sizeof(double))
#define MATRIX_BLOCK_ADDR_SKIP DOUBLE_MATRIX_BLOCK_ADDR_SKIP
#define MATRIX_DATA_ADDR_SKIP DOUBLE_MATRIX_DATA_ADDR_SKIP
#define MATRIX_SIZE DOUBLE_MATRIX_SIZE
#define BLOCK_SIZE DOUBLE_BLOCK_SIZE
/* Function macros */
#define gsl_matrix_set gsl_matrix_set
#define gsl_matrix_set_zero gsl_matrix_set_zero
#define gsl_matrix_set_identity gsl_matrix_set_identity
#define gsl_matrix_add gsl_matrix_add
#define gsl_matrix_get gsl_matrix_get
#define gsl_matrix_sub gsl_matrix_sub
#define gsl_matrix_mul_elements gsl_matrix_mul_elements
#define gsl_matrix_scale gsl_matrix_scale
#define gsl_blas_gemv gsl_blas_dgemv
#define gsl_blas_gemm gsl_blas_dgemm
#endif

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
    void fx(const String fn, const GodotGSLMatrix *a, GodotGSLMatrix *to);
    void fx(const String fn, GodotGSLMatrix *to);
    void fx(const String fn, const GodotGSLMatrix *a);
    void fx(const String fn);
    void copy(GodotGSLMatrix* to);

    size_t size[2];

private:
    DTYPE *gsl_mtx;
    VTYPE *gsl_vec;
    STYPE *gsl_scalar;
    val_type type = GSL_MATRIX;
    bool is_row_vector = false;
    DTYPE *_alloc(const size_t row_count, const size_t col_count);
    void _vector_alloc(size_t s);
    STYPE (*math_func1)(STYPE a);
    STYPE (*math_func2)(STYPE a, STYPE b);
    bool _condition(const cond_type cond, const GodotGSLMatrix *a);
    bool _condition(const cond_type cond);
    void _fx_elements1(GodotGSLMatrix *out);
    void _fx_elements2(GodotGSLMatrix *out, const GodotGSLMatrix *a);
};
static GodotGSLMatrix NULLMTX;

#endif // GODOT_GSL_MATRIX_H
