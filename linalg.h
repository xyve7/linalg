#ifndef LINALG_H
#define LINALG_H

#include <stddef.h>

typedef struct {
	double **data;
	size_t row, col;
} Mat;

Mat mat_new(size_t row, size_t col);
Mat mat_from(size_t row, size_t col, ...);
Mat mat_fill(size_t row, size_t col, double value);
Mat mat_zero(size_t row, size_t col);
Mat mat_identity(size_t row, size_t col);

Mat mat_dup(Mat *self);

double mat_at(Mat *self, size_t row, size_t col);
Mat mat_row(Mat *self, size_t row);
Mat mat_col(Mat *self, size_t col);

Mat mat_add(Mat *A, Mat *B);
Mat mat_sub(Mat *A, Mat *B);
Mat mat_mul(Mat *A, Mat *B);
Mat mat_scale(Mat *A, double scalar);

Mat mat_trans(Mat *self);
Mat mat_minor(Mat *self, size_t row, size_t col);

double mat_cofactor_at(Mat *self, size_t row, size_t col);
Mat mat_cofactor(Mat *self);
Mat mat_cramer(Mat *coefficients, Mat *constants);

double mat_det(Mat *self);
Mat mat_inverse(Mat *self);
void mat_free(Mat *self);

void mat_print(Mat *self);

#endif 
