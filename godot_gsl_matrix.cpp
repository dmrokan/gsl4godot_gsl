#include "godot_gsl_matrix.h"

void GodotGSLMatrix::init(const size_t row_count, const size_t col_count)
{
    gsl_mtx = _alloc(row_count, col_count);

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
    else if (size[0] == 1 && size[1] == 1)
    {
        type = GSL_SCALER;
        gsl_scalar = gsl_mtx->data;
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

            gsl_matrix_set(gsl_mtx, k, l, val);
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
    else if (size[0] == 1 && size[1] == 1)
    {
        type = GSL_SCALER;
        gsl_scalar = gsl_mtx->data;
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
    memcpy(mtx_raw + MATRIX_BLOCK_ADDR_SKIP, &block, sizeof(block));

    BTYPE *block_ptr = (BTYPE*) (mtx_raw + MATRIX_BLOCK_ADDR_SKIP);
    DTYPE *mtx_ptr = (DTYPE*) mtx_raw;

    mtx_ptr->data = (STYPE*) (mtx_raw + MATRIX_DATA_ADDR_SKIP);
    mtx_ptr->block = block_ptr;
    block_ptr->data = ((DTYPE*) mtx_raw)->data;
    block_ptr->size = data_size;

    mtx_ptr->size1 = row_count;
    mtx_ptr->size2 = col_count;
    mtx_ptr->tda = col_count;
    mtx_ptr->owner = 1;

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

    gsl_matrix_set_zero(gsl_mtx);
}

void GodotGSLMatrix::set_identity()
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND("GodotGSLMatrix::set_identity: gsl_mtx == NULL");
    }

    gsl_matrix_set_identity(gsl_mtx);
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

    int ret = 0;

    if (a.type == GSL_SCALER && type == GSL_SCALER)
    {
        *gsl_scalar += *(a.gsl_scalar);
    }
    else
    {
        ret = gsl_matrix_add(gsl_mtx, a.get_ptr());
    }

    return ret;
}

STYPE GodotGSLMatrix::get(const size_t row, const size_t col)
{
    if (row >= size[0] || col >= size[1])
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::get: row >= size[0] || col >= size[1]", 0.0);
    }

    STYPE ret = gsl_matrix_get(gsl_mtx, row, col);

    return ret;
}

void GodotGSLMatrix::set(const size_t row, const size_t col, const STYPE a)
{
    if (row >= size[0] || col >= size[1])
    {
        ERR_FAIL_COND("GodotGSLMatrix::set: row >= size[0] || col >= size[1]");
    }

    gsl_matrix_set(gsl_mtx, row, col, a);
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

    int ret = 0;

    if (a.type == GSL_SCALER && type == GSL_SCALER)
    {
        *gsl_scalar -= *(a.gsl_scalar);
    }
    else
    {
        ret = gsl_matrix_sub(gsl_mtx, a.get_ptr());
    }

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

    int ret = 0;

    if (a.type == GSL_SCALER && type == GSL_SCALER)
    {
        *gsl_scalar *= *(a.gsl_scalar);
    }
    else
    {
        ret = gsl_matrix_mul_elements(gsl_mtx, a.get_ptr());
    }

    return ret;
}

int GodotGSLMatrix::scale(const STYPE a)
{
    if (gsl_mtx == NULL)
    {
        ERR_FAIL_COND_V("GodotGSLMatrix::mul: gsl_mtx == NULL", ERR_NULL_VALUE);
    }

    int ret = 0;

    if (type == GSL_SCALER)
    {
        *gsl_scalar *= a;
    }
    else
    {
        ret = gsl_matrix_scale(gsl_mtx, a);
    }

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
        ERR_FAIL_COND_V("GodotGSLMatrix::prod: gsl_mtx == NULL", &NULLMTX);
    }

    if (prod_type == MATRIX_VECTOR_PROD)
    {
        if (size[0] != a.size[1])
        {
            ERR_FAIL_COND_V("GodotGSLMatrix::prod: size[0] != a.size[1]", &NULLMTX);
        }
    }
    else
    {
        if (size[0] != a.size[1] || size[1] != a.size[0])
        {
            ERR_FAIL_COND_V("GodotGSLMatrix::prod: size[0] != a.size[1] || size[1] != a.size[0]", &NULLMTX);
        }
    }

    GodotGSLMatrix *b;

    if (prod_type == MATRIX_VECTOR_PROD)
    {
        b = memnew(GodotGSLMatrix);

        if (type == GSL_VECTOR)
        {
            if (!is_row_vector)
            {
                b->init(size[1], a.size[0]);
                gsl_blas_gemv(CblasNoTrans, 1.0, a.get_ptr(), gsl_vec, 0.0, b->get_vec_ptr());
            }
            else
            {
                b->init(size[0], a.size[1]);
                gsl_blas_gemv(CblasTrans, 1.0, a.get_ptr(), gsl_vec, 0.0, b->get_vec_ptr());
            }
        }
        else if (a.type == GSL_VECTOR)
        {
            if (!a.is_row_vector)
            {
                b->init(a.size[1], size[0]);
                gsl_blas_gemv(CblasNoTrans, 1.0, gsl_mtx, a.get_vec_ptr(), 0.0, b->get_vec_ptr());
            }
            else
            {
                b->init(a.size[0], size[1]);
                gsl_blas_gemv(CblasTrans, 1.0, gsl_mtx, a.get_vec_ptr(), 0.0, b->get_vec_ptr());
            }
        }
    }
    else if (prod_type == MATRIX_MATRIX_PROD || prod_type == VECTOR_VECTOR_PROD)
    {
        b = memnew(GodotGSLMatrix);

        b->init(a.size[0], size[1]);
        gsl_blas_gemm(CblasNoTrans, CblasNoTrans, 1.0, a.get_ptr(), gsl_mtx, 0.0, b->get_ptr());
    }
    else
    {
        WARN_PRINT("GodotGSLMatrix::prod: For scalar products use matrix_scale function.");
    }

    return b;
}

void GodotGSLMatrix::fx(const String fn, const GodotGSLMatrix *a, GodotGSLMatrix *to)
{
    printf("aaaaaaaaaaaaaaaaaaaaaaaaaaa %d %d %d\n", a->size[0], to->size[1]);
    if (fn == FN_SIN)
    {
        if (_condition(EQUAL_SIZE, to))
        {
            math_func1 = sin;
            _fx_elements1(to);
            return;
        }
    }
    else if (fn == FN_ADD)
    {
        if (_condition(EQUAL_SIZE, a))
        {
            if (this != to)
            {
                this->copy(to);
            }
            else if (!_condition(EQUAL_SIZE, to))
            {
                return;
            }

            to->add(*a);
            return;
        }
    }
}

void GodotGSLMatrix::fx(const String fn, const GodotGSLMatrix *a)
{
    fx(fn, a, this);
}

void GodotGSLMatrix::fx(const String fn, GodotGSLMatrix *to)
{
    fx(fn, NULL, to);
}

void GodotGSLMatrix::fx(const String fn)
{
    fx(fn, NULL, this);
}

void GodotGSLMatrix::_fx_elements1(GodotGSLMatrix *out)
{
    if (math_func1 == NULL)
    {
        ERR_FAIL_COND("GodotGSLMatrix::_fx_elements1: math_func1 == NULL");
    }

    for (size_t k = 0; k < size[0]; k++)
    {
        for (size_t l = 0; l < size[1]; l++)
        {
            STYPE val = get(k, l);
            val = math_func1(val);
            out->set(k, l, val);
        }
    }

    math_func1 = NULL;
}

bool GodotGSLMatrix::_condition(const cond_type cond)
{
    return _condition(cond, NULL);
}

bool GodotGSLMatrix::_condition(const cond_type cond, const GodotGSLMatrix *a)
{
    if (cond == EQUAL_SIZE)
    {
        if (size[0] == a->size[0] && size[1] == a->size[1])
        {
            return true;
        }
    }
    else if (cond == NONZERO_SIZE)
    {
        if (size[0] > 0 && size[1] > 0)
        {
            return true;
        }
    }

    return false;
}

void GodotGSLMatrix::copy(GodotGSLMatrix* to)
{
    if (to == NULL)
    {
        GGSL_MESSAGE("GodotGSLMatrix::copy: to == NULL");
        return;
    }

    gsl_matrix_memcpy(to->gsl_mtx, this->gsl_mtx);
}
