/* blas/source_gemm_c.h
 * 
 * Copyright (C) 2001 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

{
    size_t i, j, k;
    size_t n1, n2;
    size_t ldf, ldg;
    int conjF, conjG, TransF, TransG;
    const BASE *F, *G;

    const BASE alpha_real = REAL0(alpha), alpha_imag = IMAG0(alpha);
    const BASE beta_real = REAL0(beta), beta_imag = IMAG0(beta);

    if ((alpha_real == 0.0 && alpha_imag == 0.0)
	&& (beta_real == 1.0 && beta_imag == 0.0))
	return;

    if (Order == CblasRowMajor)
      {
        n1 = M;
        n2 = N;
        F = A;
        ldf = lda;
        conjF = (TransA == CblasConjTrans) ? -1 : 1;
        TransF = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
        G = B;
        ldg = ldb;
        conjG = (TransB == CblasConjTrans) ? -1 : 1;
        TransG = (TransB == CblasNoTrans) ? CblasNoTrans : CblasTrans;
      }
    else
      {
        n1 = N;
        n2 = M;
        F = B;
        ldf = ldb;
        conjF = (TransB == CblasConjTrans) ? -1 : 1;
        TransF = (TransB == CblasNoTrans) ? CblasNoTrans : CblasTrans;
        G = A;
        ldg = lda;
        conjG = (TransA == CblasConjTrans) ? -1 : 1;
        TransG = (TransA == CblasNoTrans) ? CblasNoTrans : CblasTrans;
      }

    /* form  y := beta*y */
    if (beta_real == 0.0 && beta_imag == 0.0) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          REAL(C, ldc * i + j) = 0.0;
          IMAG(C, ldc * i + j) = 0.0;
        }
      }
    } else if (!(beta_real == 1.0 && beta_imag == 0.0)) {
      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          const BASE Cij_real = REAL(C, ldc * i + j);
          const BASE Cij_imag = IMAG(C, ldc * i + j);
          REAL(C, ldc * i + j) = beta_real * Cij_real - beta_imag * Cij_imag;
          IMAG(C, ldc * i + j) = beta_real * Cij_imag + beta_imag * Cij_real;
        }
      }
    }

    if (alpha_real == 0.0 && alpha_imag == 0.0)
	return;

    if (TransF == CblasNoTrans && TransG == CblasNoTrans) {

	/* form  C := alpha*A*B + C */

      for (k = 0; k < K; k++) {
        for (i = 0; i < n1; i++) {
          const BASE Fik_real =  REAL(F, ldf * i + k);
          const BASE Fik_imag =  conjF * IMAG(F, ldf * i + k);
          const BASE temp_real = alpha_real * Fik_real - alpha_imag * Fik_imag;
          const BASE temp_imag = alpha_real * Fik_imag + alpha_imag * Fik_real;
          if (!(temp_real == 0.0 && temp_imag == 0.0)) {
            for (j = 0; j < n2; j++) {
              const BASE Gkj_real = REAL(G,ldg * k + j);
              const BASE Gkj_imag = conjG * IMAG(G,ldg * k + j);
              REAL(C,ldc * i + j) += temp_real * Gkj_real - temp_imag * Gkj_imag;
              IMAG(C,ldc * i + j) += temp_real * Gkj_imag + temp_imag * Gkj_real;
            }
          }
        }
      }
      
    } else if (TransF == CblasNoTrans && TransG == CblasTrans) {

	/* form  C := alpha*A*B' + C */

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          BASE temp_real = 0.0;
          BASE temp_imag = 0.0;
          for (k = 0; k < K; k++) {
            const BASE Fik_real = REAL(F,ldf * i + k);
            const BASE Fik_imag = conjF * IMAG(F,ldf * i + k);
            const BASE Gjk_real = REAL(G,ldg * j + k);
            const BASE Gjk_imag = conjG * IMAG(G,ldg * j + k);
            temp_real +=  Fik_real * Gjk_real - Fik_imag * Gjk_imag ;
            temp_imag +=  Fik_real * Gjk_imag + Fik_imag * Gjk_real ;
          }
          REAL(C,ldc * i + j) += alpha_real * temp_real - alpha_imag * temp_imag;
          IMAG(C,ldc * i + j) += alpha_real * temp_imag + alpha_imag * temp_real;
        }
      }

    } else if (TransF == CblasTrans && TransG == CblasNoTrans) {

      for (k = 0; k < K; k++) {
        for (i = 0; i < n1; i++) {
          const BASE Fki_real =  REAL(F, ldf * k + i);
          const BASE Fki_imag =  conjF * IMAG(F, ldf * k + i);
          const BASE temp_real = alpha_real * Fki_real - alpha_imag * Fki_imag;
          const BASE temp_imag = alpha_real * Fki_imag + alpha_imag * Fki_real;
          if (!(temp_real == 0.0 && temp_imag == 0.0)) {
            for (j = 0; j < n2; j++) {
              const BASE Gkj_real = REAL(G,ldg * k + j);
              const BASE Gkj_imag = conjG * IMAG(G,ldg * k + j);
              REAL(C,ldc * i + j) += temp_real * Gkj_real - temp_imag * Gkj_imag;
              IMAG(C,ldc * i + j) += temp_real * Gkj_imag + temp_imag * Gkj_real;
            }
          }
        }
      }

    } else if (TransF == CblasTrans && TransG == CblasTrans) {

      for (i = 0; i < n1; i++) {
        for (j = 0; j < n2; j++) {
          BASE temp_real = 0.0;
          BASE temp_imag = 0.0;
          for (k = 0; k < K; k++) {
            const BASE Fki_real = REAL(F,ldf * k + i);
            const BASE Fki_imag = conjF * IMAG(F,ldf * k + i);
            const BASE Gjk_real = REAL(G,ldg * j + k);
            const BASE Gjk_imag = conjG * IMAG(G,ldg * j + k);

            temp_real +=  Fki_real * Gjk_real - Fki_imag * Gjk_imag ;
            temp_imag +=  Fki_real * Gjk_imag + Fki_imag * Gjk_real ;
          }
          REAL(C,ldc * i + j) += alpha_real * temp_real - alpha_imag * temp_imag;
          IMAG(C,ldc * i + j) += alpha_real * temp_imag + alpha_imag * temp_real;
        }
      }

    } else {
	BLAS_ERROR("unrecognized operation");
    }
}