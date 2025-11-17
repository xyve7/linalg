#include <linalg.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void die(const char *restrict format, ...) {
    va_list args;
    va_start(args, format);

    fprintf(stderr, "linalg: ");
    vfprintf(stderr, format, args);
    fputc('\n', stderr);

    exit(EXIT_FAILURE);

    va_end(args);
}

Mat mat_new(size_t row, size_t col) {
    Mat mat;
    mat.row = row;
    mat.col = col;

    mat.data = calloc(row, sizeof(double *));
    for (size_t i = 0; i < row; i++) {
        mat.data[i] = calloc(col, sizeof(double));
    }
    return mat;
}
Mat mat_from(size_t row, size_t col, ...) {
    Mat mat = mat_new(row, col);

    va_list values;
    va_start(values, col);
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
            mat.data[i][j] = va_arg(values, double);
        }
    }
    va_end(values);
    return mat;
}
Mat mat_fill(size_t row, size_t col, double value) {
    Mat mat = mat_new(row, col);
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
            mat.data[i][j] = value;
        }
    }
    return mat;
}
Mat mat_zero(size_t row, size_t col) {
    return mat_fill(row, col, 0.0);
}
Mat mat_identity(size_t row, size_t col) {
    if (row != col) {
        die("mat_identity: row != col");
    }
    Mat mat = mat_zero(row, col);
    for (size_t i = 0; i < row; i++) {
        mat.data[i][i] = 1;
    }
    return mat;
}

Mat mat_dup(Mat *self) {
    Mat mat;
    mat.row = self->row;
    mat.col = self->col;

    mat.data = calloc(mat.row, sizeof(double *));
    for (size_t i = 0; i < mat.row; i++) {
        mat.data[i] = calloc(mat.col, sizeof(double));
        memcpy(mat.data[i], self->data[i], self->col * sizeof(double));
    }
    return mat;
}

double mat_at(Mat *self, size_t row, size_t col) {
    if (row < 1 || col < 1) {
        die("mat_at: mat is indexed from 1, not 0");
    }
    if (row > self->row || col > self->col) {
        die("mat_at: index out of bounds");
    }
    return self->data[row - 1][col - 1];
}
Mat mat_row(Mat *self, size_t row) {
    if (row < 1) {
        die("mat_row: mat is indexed from 1, not 0");
    }
    Mat mat = mat_new(1, self->col);
    memcpy(mat.data[0], self->data[row - 1], self->col * sizeof(double));
    return mat;
}
Mat mat_col(Mat *self, size_t col) {
    if (col < 1) {
        die("mat_row: mat is indexed from 1, not 0");
    }
    Mat mat = mat_new(self->row, 1);
    for (size_t i = 0; i < mat.row; i++) {
        mat.data[i][0] = self->data[i][col - 1];
    }
    return mat;
}

Mat mat_add(Mat *A, Mat *B) {
    if (A->row != B->row || A->col != B->col) {
        die("mat_add: mat A and B must be the same size");
    }
    Mat mat = mat_new(A->row, A->col);
    for (size_t i = 0; i < A->row; i++) {
        for (size_t j = 0; j < A->col; j++) {
            mat.data[i][j] = A->data[i][j] + B->data[i][j];
        }
    }
    return mat;
}
Mat mat_sub(Mat *A, Mat *B) {
    if (A->row != B->row || A->col != B->col) {
        die("mat_add: mat A and B must be the same size");
    }
    Mat mat = mat_new(A->row, A->col);
    for (size_t i = 0; i < A->row; i++) {
        for (size_t j = 0; j < A->col; j++) {
            mat.data[i][j] = A->data[i][j] - B->data[i][j];
        }
    }
    return mat;
}
Mat mat_mul(Mat *A, Mat *B) {
    if (A->col != B->row) {
        die("mat_mul: mat A and B must be complimentary");
    }
    size_t l = A->col;
    Mat mat = mat_new(A->row, B->col);
    for (size_t i = 0; i < mat.row; i++) {
        for (size_t j = 0; j < mat.col; j++) {
            mat.data[i][j] = 0;
            for (size_t k = 0; k < l; k++) {
                mat.data[i][j] += A->data[i][k] * B->data[k][j];
            }
        }
    }
    return mat;
}
Mat mat_scale(Mat *A, double scalar) {
    Mat mat = mat_new(A->row, A->col);
    for (size_t i = 0; i < A->row; i++) {
        for (size_t j = 0; j < A->col; j++) {
            mat.data[i][j] = scalar * A->data[i][j];
        }
    }
    return mat;
}

Mat mat_trans(Mat *self) {
    Mat mat = mat_new(self->col, self->row);
    for (size_t i = 0; i < self->row; i++) {
        for (size_t j = 0; j < self->col; j++) {
            mat.data[j][i] = self->data[i][j];
        }
    }
    return mat;
}
Mat mat_minor(Mat *self, size_t row, size_t col) {
    Mat mat = mat_new(self->row - 1, self->col - 1);
    size_t k = 0;
    size_t l = 0;
    for (size_t i = 0; i < self->row; i++) {
        if (i == row - 1) {
            continue;
        }
        for (size_t j = 0; j < self->col; j++) {
            if (j == col - 1) {
                continue;
            }
            mat.data[k][l] = self->data[i][j];
            l++;
        }
        k++;
        l = 0;
    }
    return mat;
}

double mat_cofactor_at(Mat *self, size_t row, size_t col) {
    Mat minor = mat_minor(self, row, col);
	int sign = ((row + col) % 2 == 0) ? 1 : -1;
    double cofactor = (double)sign * mat_det(&minor);
    mat_free(&minor);
    return cofactor;
}
Mat mat_cofactor(Mat *self) {
    if (self->col == 2 && self->row == 2) {
        Mat mat = mat_from(2, 2, self->data[1][1], -self->data[1][0], -self->data[0][1], self->data[0][0]);
        return mat;
    }
    Mat mat = mat_new(self->row, self->col);
    for (size_t i = 0; i < self->row; i++) {
        for (size_t j = 0; j < self->col; j++) {
            mat.data[i][j] = mat_cofactor_at(self, i + 1, j + 1);
        }
    }
    return mat;
}
Mat mat_cramer(Mat *coefficients, Mat *constants) {
	if (coefficients->row != coefficients->col) {
		die("mat_cramer: mat must be square");
	}
	if (coefficients->col != constants->row) {
		die("mat_cramer: constants are not equal to coefficients");
	}

	double coefficient_det = mat_det(coefficients); 
	Mat result = mat_new(constants->row, 1);
	for (size_t i = 0; i < constants->row; i++) {
		Mat A = mat_dup(coefficients);

		for (size_t j = 0; j < constants->row; j++) {
			A.data[j][i] = constants->data[j][0];
		}

		double A_det = mat_det(&A);
		result.data[i][0] = A_det / coefficient_det;

		mat_free(&A);
	}
	return result;
}
double mat_det(Mat *self) {
    if (self->col != self->row) {
        die("mat_det: mat must be square");
    }
    if (self->col == 2 && self->row == 2) {
        double a = mat_at(self, 1, 1);
        double b = mat_at(self, 1, 2);
        double c = mat_at(self, 2, 1);
        double d = mat_at(self, 2, 2);

        return a * d - b * c;
    }
    double det = 0;
    size_t i = 1;
    for (size_t j = 1; j <= self->col; j++) {
        double alpha = mat_at(self, i, j);
        det += alpha * mat_cofactor_at(self, i, j);
    }
    return det;
}
Mat mat_inverse(Mat *self) {
    double det = mat_det(self);
    if (det == 0.0) {
        die("mat_inverse: mat not invertable");
    }
    Mat C = mat_cofactor(self);
    Mat Ct = mat_trans(&C);

    Mat inv = mat_scale(&Ct, 1.0 / det);

    mat_free(&C);
    mat_free(&Ct);

    return inv;
}
void mat_free(Mat *self) {
    for (size_t i = 0; i < self->row; i++) {
        free(self->data[i]);
    }
    free(self->data);
    self->col = 0;
    self->row = 0;
}
void mat_print(Mat *self) {
    for (size_t i = 0; i < self->row; i++) {
        printf("[");
        for (size_t j = 0; j < self->col; j++) {
            printf("%3.2f", self->data[i][j]);
            if (j + 1 != self->col) {
                putc(' ', stdout);
            }
        }
        printf("]\n");
    }
}
