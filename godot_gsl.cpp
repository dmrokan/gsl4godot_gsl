/* godot_gsl.cpp */

#include "godot_gsl.h"

void GodotGSL::_bind_methods() {
    ClassDB::bind_method(D_METHOD("matrix_new", "vn", "row_count", "col_count"), &GodotGSL::matrix_new);
    ClassDB::bind_method(D_METHOD("matrix_new_from_array", "vn", "a"), &GodotGSL::matrix_new_from_array);
    ClassDB::bind_method(D_METHOD("add", "vn1", "vn2"), &GodotGSL::add);
    ClassDB::bind_method(D_METHOD("to_array", "vn"), &GodotGSL::to_array);
    ClassDB::bind_method(D_METHOD("set_zero", "vn"), &GodotGSL::set_zero);
    ClassDB::bind_method(D_METHOD("prod", "vn1", "vn2", "to"), &GodotGSL::prod);
    ClassDB::bind_method(D_METHOD("vector_new", "vn", "size", "is_row_vector"), &GodotGSL::vector_new, NULL, NULL, false);
    ClassDB::bind_method(D_METHOD("vector_new_from_array", "vn", "a", "is_row_vector"), &GodotGSL::vector_new_from_array, NULL, NULL, false);
}

GodotGSL::GodotGSL() {
    return;
}

void GodotGSL::_add_variable(String vn, GodotGSLMatrix* mtx)
{
    if (variables.has(vn))
    {
        memdelete(variables[vn]);
        variables[vn] = mtx;
    }

    variables.insert(vn, mtx);
}

void GodotGSL::matrix_new(String vn, const size_t row_count, const size_t col_count)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);
    mtx->init(row_count, col_count);

    _add_variable(vn, mtx);
}

void GodotGSL::matrix_new_from_array(String vn, const Array a)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);
    mtx->init(a);

    _add_variable(vn, mtx);
}

void GodotGSL::add(String vn1, String vn2)
{
    GodotGSLMatrix *mtx1 = variables[vn1];
    GodotGSLMatrix *mtx2 = variables[vn2];

    mtx1->add(*mtx2);
}

void GodotGSL::set_zero(String vn)
{
    GodotGSLMatrix *mtx = variables[vn];

    mtx->set_zero();
}

Array GodotGSL::to_array(String vn)
{
    GodotGSLMatrix *mtx = variables[vn];

    size_t row_count = mtx->size[0];
    size_t col_count = mtx->size[1];

    Array rows;
    for (size_t k = 0; k < row_count; k++)
    {
        Array row;
        for (size_t l = 0; l < col_count; l++)
        {
            row.append(mtx->get(k, l));
        }

        rows.append(row);
    }

    return rows;
}

void GodotGSL::prod(String vn1, String vn2, String to)
{
    GodotGSLMatrix *mtx1 = variables[vn1];
    GodotGSLMatrix *mtx2 = variables[vn2];

    GodotGSLMatrix *result = mtx2->prod(*mtx1);
    _add_variable(to, result);
}


void GodotGSL::vector_new(String vn, const size_t size, bool is_row_vector)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);

    size_t row_count = size;
    size_t col_count = 1;

    if (is_row_vector)
    {
        row_count = 1;
        col_count = size;
    }

    mtx->init(row_count, col_count);

    _add_variable(vn, mtx);
}

void GodotGSL::vector_new_from_array(String vn, const Array a, bool is_row_vector)
{
    GodotGSLMatrix *mtx = memnew(GodotGSLMatrix);

    Array b;

    if (is_row_vector)
    {
        b.append(a);
    }
    else
    {
        size_t size = a.size();
        for (size_t k = 0; k < size; k++)
        {
            STYPE value = a.get(k);
            Array c;
            c.append(value);
            b.append(c);
        }
    }

    mtx->init(b);

    _add_variable(vn, mtx);
}

/*
*********** Matrix class
*/

#define DTYPE gsl_matrix_float
#define VTYPE gsl_vector_float
#define BTYPE gsl_block_float
#define STYPE float
#define STRIDE 1

#if DTYPE == gsl_matrix_float
#define BLOCK_SIZE FLOAT_BLOCK_SIZE
#define MATRIX_SIZE FLOAT_MATRIX_SIZE
#define MULTIPLICITY 1
#define ATOMIC_SIZE (sizeof(float))
#define MATRIX_DATAPTR_SKIP FLOAT_MATRIX_DATAPTR_SKIP
#define MATRIX_BLOCKPTR_SKIP FLOAT_MATRIX_BLOCKPTR_SKIP
#define MATRIX_DATAPTR2_SKIP FLOAT_MATRIX_DATAPTR2_SKIP
#endif

void GodotGSLMatrix::init(const size_t row_count, const size_t col_count)
{
    gsl_mtx = _alloc(row_count, col_count);

    // gsl_mtx = gsl_matrix_float_alloc(row_count, col_count);

    size[0] = row_count;
    size[1] = col_count;

    if ((size[0] == 1 || size[1] == 1) && !(size[0] == 1 && size[1] == 1))
    {
        type = GSL_VECTOR;

        if (size[0] == 1)
        {
            is_row_vector = true;
        }

        _vector_alloc((size[0] == 1 ? size[1] : size[0]));
    }
}

void GodotGSLMatrix::init(const Array a)
{
    size_t row_count = a.size();
    if (row_count == 0)
    {
        ERR_FAIL_COND("GodotGSLMatrix::GodotGSLMatrix: row_count == 0");
    }

    Array b = (Array) a.get(0);
    size_t col_count = b.size();
    if (col_count == 0)
    {
        ERR_FAIL_COND("GodotGSLMatrix::GodotGSLMatrix: col_count == 0");
    }

    gsl_mtx = _alloc(row_count, col_count);

    for (int k = 0; k < row_count; k++)
    {
        for (int l = 0; l < col_count; l++)
        {
            Array row = (Array) a.get(k);
            if (row.size() <= l)
            {
                ERR_FAIL_COND("GodotGSLMatrix::GodotGSLMatrix: row->size() <= l");
            }

            STYPE val = row.get(l);

            gsl_matrix_float_set(gsl_mtx, k, l, val);
        }
    }

    size[0] = row_count;
    size[1] = col_count;

    if ((size[0] == 1 || size[1] == 1) && !(size[0] == 1 && size[1] == 1))
    {
        type = GSL_VECTOR;

        if (size[0] == 1)
        {
            is_row_vector = true;
        }

        _vector_alloc((size[0] == 1 ? size[1] : size[0]));
    }
}

GodotGSLMatrix::~GodotGSLMatrix()
{
    memdelete_arr((uint8_t*) gsl_mtx);
    memdelete(gsl_vec);
}


DTYPE *GodotGSLMatrix::_alloc(const size_t row_count, const size_t col_count)
{
    size_t data_size = row_count * col_count;
    size_t total_bytes = MATRIX_SIZE;
    total_bytes += BLOCK_SIZE;
    total_bytes += MULTIPLICITY * data_size * ATOMIC_SIZE;
    uint8_t *mtx_raw = memnew_arr(uint8_t, total_bytes);

    BTYPE block;
    DTYPE mtx;
    memcpy(mtx_raw, &mtx, sizeof(mtx));
    memcpy(mtx_raw + FLOAT_MATRIX_BLOCK_ADDR_SKIP, &block, sizeof(block));

    BTYPE *block_ptr = (BTYPE*) (mtx_raw + FLOAT_MATRIX_BLOCK_ADDR_SKIP);
    DTYPE *mtx_ptr = (DTYPE*) mtx_raw;

    mtx_ptr->data = (STYPE*) (mtx_raw + FLOAT_MATRIX_DATA_ADDR_SKIP);
    mtx_ptr->block = block_ptr;
    block_ptr->data = ((DTYPE*) mtx_raw)->data;
    block_ptr->size = data_size;

    mtx_ptr->size1 = row_count;
    mtx_ptr->size2 = col_count;
    mtx_ptr->tda = col_count;
    mtx_ptr->owner = 1;

    /*
    DTYPE *mtx = gsl_matrix_float_alloc(row_count, col_count);
    */

    return mtx_ptr;
}

void GodotGSLMatrix::_vector_alloc(size_t s)
{
    gsl_vec = (VTYPE*) memnew(VTYPE);
    gsl_vec->size = s;
    gsl_vec->stride = STRIDE;
    gsl_vec->data = gsl_mtx->data;
    gsl_vec->block = gsl_mtx->block;
    gsl_vec->owner = gsl_mtx->owner;
}

void GodotGSLMatrix::set_zero()
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND("GodotGSLMatrix::set_zero: gsl_mtx == NULL");
    }

    gsl_matrix_float_set_zero(gsl_mtx);
}

void GodotGSLMatrix::set_identity()
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND("GodotGSLMatrix::set_identity: gsl_mtx == NULL");
    }

    gsl_matrix_float_set_identity(gsl_mtx);
}

int GodotGSLMatrix::add(const GodotGSLMatrix &a)
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::add: gsl_mtx == NULL", ERR_NULL_VALUE);
    }

    if (size[0] != a.size[0] || size[1] != a.size[1])
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::add: size[0] != a.size[0] || size[1] != a.size[1]", ERR_DIMENSION_MISMATCH);
    }

    int ret = gsl_matrix_float_add(gsl_mtx, a.get_ptr());

    return ret;
}

STYPE GodotGSLMatrix::get(const size_t row, const size_t col)
{
    if (row >= size[0] || col >= size[1])
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::get: row >= size[0] || col >= size[1]", 0.0);
    }

    STYPE ret = gsl_matrix_float_get(gsl_mtx, row, col);

    return ret;
}

void GodotGSLMatrix::set(const size_t row, const size_t col, const STYPE a)
{
    if (row >= size[0] || col >= size[1])
    {
        ERR_FAIL_COND("GodotGSLMatrix::set: row >= size[0] || col >= size[1]");
    }

    gsl_matrix_float_set(gsl_mtx, row, col, a);
}

DTYPE *GodotGSLMatrix::get_ptr() const
{
    return gsl_mtx;
}

VTYPE *GodotGSLMatrix::get_vec_ptr() const
{
    return gsl_vec;
}

int GodotGSLMatrix::sub(const GodotGSLMatrix &a)
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::sub: gsl_mtx == NULL", ERR_NULL_VALUE);
    }

    if (size[0] != a.size[0] || size[1] != a.size[1])
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::sub: size[0] != a.size[0] || size[1] != a.size[1]", ERR_DIMENSION_MISMATCH);
    }

    int ret = gsl_matrix_float_add(gsl_mtx, a.get_ptr());

    return ret;
}

int GodotGSLMatrix::mul_elements(const GodotGSLMatrix &a)
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::mul: gsl_mtx == NULL", ERR_NULL_VALUE);
    }

    if (size[0] != a.size[0] || size[1] != a.size[1])
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::mul: size[0] != a.size[0] || size[1] != a.size[1]", ERR_DIMENSION_MISMATCH);
    }

    int ret = gsl_matrix_float_mul_elements(gsl_mtx, a.get_ptr());

    return ret;
}

int GodotGSLMatrix::scale(const STYPE a)
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::mul: gsl_mtx == NULL", ERR_NULL_VALUE);
    }

    int ret = gsl_matrix_float_scale(gsl_mtx, a);

    return ret;
}


GodotGSLMatrix *GodotGSLMatrix::prod(const GodotGSLMatrix &a)
{
    int prod_type = -1;

    if ((type == GSL_MATRIX && a.type == GSL_MATRIX)
        || (type == GSL_VECTOR && a.type == GSL_VECTOR))
    {
        prod_type = MATRIX_MATRIX_PROD;
    }
    else
    {
        prod_type = MATRIX_VECTOR_PROD;
    }

    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::prod: gsl_mtx == NULL", NULLMTX);
    }

    if (prod_type == MATRIX_VECTOR_PROD)
    {
        if (size[0] != a.size[1])
        {
            ERR_FAIL_COND_V("GodotGSLMatrix::prod: size[0] != a.size[1]", NULLMTX);
        }
    }
    else
    {
        if (size[0] != a.size[1] || size[1] != a.size[0])
        {
            ERR_FAIL_COND_V("GodotGSLMatrix::prod: size[0] != a.size[1] || size[1] != a.size[0]", NULLMTX);
        }
    }

    GodotGSLMatrix *b = memnew(GodotGSLMatrix);

    if (prod_type == MATRIX_VECTOR_PROD)
    {
        if (type == GSL_VECTOR)
        {
            if (!is_row_vector)
            {
                b->init(size[1], a.size[0]);
                gsl_blas_sgemv(CblasNoTrans, 1.0, a.get_ptr(), gsl_vec, 0.0, b->get_vec_ptr());
            }
            else
            {
                b->init(size[0], a.size[1]);
                gsl_blas_sgemv(CblasTrans, 1.0, a.get_ptr(), gsl_vec, 0.0, b->get_vec_ptr());
            }
        }
        else if (a.type == GSL_VECTOR)
        {
            if (!a.is_row_vector)
            {
                b->init(a.size[1], size[0]);
                gsl_blas_sgemv(CblasNoTrans, 1.0, gsl_mtx, a.get_vec_ptr(), 0.0, b->get_vec_ptr());
            }
            else
            {
                b->init(a.size[0], size[1]);
                gsl_blas_sgemv(CblasTrans, 1.0, gsl_mtx, a.get_vec_ptr(), 0.0, b->get_vec_ptr());
            }
        }
    }
    else
    {
        b->init(a.size[0], size[1]);
        gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1.0, a.get_ptr(), gsl_mtx, 0.0, b->get_ptr());
    }

    return b;
}

#undef DTYPE
#undef BTYPE
#undef STYPE
