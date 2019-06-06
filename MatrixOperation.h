#pragma once
#ifndef MAREIXOPERATION
#define MAREIXOPERATION


/* Matrix math. */

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"ParaDefine.h"


/* 为矩阵分配初始空间 */
Matrix alloc_matrix(int rows, int cols) {
	Matrix m;
	int i;
	int j;
	m.rows = rows;
	m.cols = cols;
	m.data = (float**)malloc(sizeof(float*) * m.rows);

	for (i = 0; i < m.rows; ++i)
	{
		m.data[i] = (float*)malloc(sizeof(float) * m.cols);
		assert(m.data[i]);
		for (j = 0; j < m.cols; ++j) {
			m.data[i][j] = 0.0;
		}
	}
	return m;
}
/* 释放空间 */
void free_matrix(Matrix m) {
	int i;
	assert(m.data != NULL);
	for (i = 0; i < m.rows; ++i) {
		free(m.data[i]);
	}
	free(m.data);
}
/* 初始化矩阵 */
void set_matrix(Matrix m, ...) {
	va_list ap;
	int i, j;
	va_start(ap, m);

	for (i = 0; i < m.rows; ++i) {
		for (j = 0; j < m.cols; ++j) {
			m.data[i][j] = va_arg(ap, float);
		}
	}

	va_end(ap);
}
/* 转换为单元矩阵 */
void set_identity_matrix(Matrix m) {
	int i;
	int j;
	assert(m.rows == m.cols);
	for (i = 0; i < m.rows; ++i) {
		for (j = 0; j < m.cols; ++j) {
			if (i == j) {
				m.data[i][j] = 1.0;
			}
			else {
				m.data[i][j] = 0.0;
			}
		}
	}
}
/* 复制矩阵 */
void copy_matrix(Matrix source, Matrix destination) {
	int i;
	int j;
	assert(source.rows == destination.rows);
	assert(source.cols == destination.cols);
	for (i = 0; i < source.rows; ++i) {
		for (j = 0; j < source.cols; ++j) {
			destination.data[i][j] = source.data[i][j];
		}
	}
}
/* 打印矩阵 */
void print_matrix(Matrix m) {
	int i;
	int j;
	for (i = 0; i < m.rows; ++i) {
		for (j = 0; j < m.cols; ++j) {
			if (j > 0) {
				printf(" ");
			}
			printf("%f ", m.data[i][j]);
		}
		printf("\n");
	}
}
/* 矩阵相加 */
Matrix add_matrix(Matrix a, Matrix b) {
	Matrix c;
	int i;
	int j;
	assert(a.rows == b.rows);
	assert(a.cols == b.cols);
	c = alloc_matrix(a.rows, a.cols);
	for (i = 0; i < a.rows; ++i) {
		for (j = 0; j < a.cols; ++j) {
			c.data[i][j] = a.data[i][j] + b.data[i][j];
		}
	}
	return c;
}
/* 矩阵相减 */
void subtract_matrix(Matrix a, Matrix b, Matrix c) {
	int i;
	int j;
	assert(a.rows == b.rows);
	assert(a.rows == c.rows);
	assert(a.cols == b.cols);
	assert(a.cols == c.cols);
	for (i = 0; i < a.rows; ++i) {
		for (j = 0; j < a.cols; ++j) {
			c.data[i][j] = a.data[i][j] - b.data[i][j];
		}
	}
}
/* 用单元矩阵减去该矩阵 */
void subtract_from_identity_matrix(Matrix a) {
	int i;
	int j;
	assert(a.rows == a.cols);
	for (i = 0; i < a.rows; ++i) {
		for (j = 0; j < a.cols; ++j) {
			if (i == j) {
				a.data[i][j] = 1.0 - a.data[i][j];
			}
			else {
				a.data[i][j] = 0.0 - a.data[i][j];
			}
		}
	}
}
/* 矩阵相乘 */
void multiply_matrix(Matrix a, Matrix b, Matrix c) {
	int i;
	int j;
	int k;
	assert(a.cols == b.rows);
	assert(a.rows == c.rows);
	assert(b.cols == c.cols);
	for (i = 0; i < c.rows; ++i) {
		for (j = 0; j < c.cols; ++j) {
			/* Calculate element c.data[i][j] via a dot product of one row of a
			with one column of b */
			c.data[i][j] = 0.0;
			for (k = 0; k < a.cols; ++k) {
				c.data[i][j] += a.data[i][k] * b.data[k][j];
			}
		}
	}
}

/* 乘以一个矩阵的转置矩阵. */
void multiply_by_transpose_matrix(Matrix a, Matrix b, Matrix c) {
	int i;
	int j;
	int k;
	assert(a.cols == b.cols);
	assert(a.rows == c.rows);
	assert(b.rows == c.cols);
	for (i = 0; i < c.rows; ++i) {
		for (j = 0; j < c.cols; ++j) {
			/* Calculate element c.data[i][j] via a dot product of one row of a
			with one row of b */
			c.data[i][j] = 0.0;
			for (k = 0; k < a.cols; ++k) {
				c.data[i][j] += a.data[i][k] * b.data[j][k];
			}
		}
	}
}
/* 矩阵转置 */
void transpose_matrix(Matrix input, Matrix output) {
	int i;
	int j;
	//int k;
	assert(input.rows == output.cols);
	assert(input.cols == output.rows);
	for (i = 0; i < input.rows; ++i) {
		for (j = 0; j < input.cols; ++j) {
			output.data[j][i] = input.data[i][j];
		}
	}
}
/* 两矩阵是否相等 */
int equal_matrix(Matrix a, Matrix b, float tolerance) {
	int i;
	int j;
	//int k;
	assert(a.rows == b.rows);
	assert(a.cols == b.cols);
	for (i = 0; i < a.rows; ++i) {
		for (j = 0; j < a.cols; ++j) {
			if (abs(a.data[i][j] - b.data[i][j]) > tolerance) {
				return 0;
			}
		}
	}
	return 1;
}
/* 矩阵乘以一个系数 */
void scale_matrix(Matrix m, float scalar) {
	int i;
	int j;

	assert(scalar != 0.0);
	for (i = 0; i < m.rows; ++i) {
		for (j = 0; j < m.cols; ++j) {
			m.data[i][j] *= scalar;
		}
	}
}

/* 交换矩阵的两行 */
void swap_rows(Matrix m, int r1, int r2) {
	float *tmp;
	assert(r1 != r2);
	tmp = m.data[r1];
	m.data[r1] = m.data[r2];
	m.data[r2] = tmp;
}
/* 矩阵某行乘以一个系数 */
void scale_row(Matrix m, int r, float scalar) {
	int i;
	assert(scalar != 0.0);
	for (i = 0; i < m.cols; ++i) {
		m.data[r][i] *= scalar;
	}
}

/* Add scalar * row r2 to row r1. */
void shear_row(Matrix m, int r1, int r2, float scalar) {
	int i;
	assert(r1 != r2);
	for (i = 0; i < m.cols; ++i) {
		m.data[r1][i] += scalar * m.data[r2][i];
	}
}

/* 矩阵的求逆(借鉴他人) */
/* Uses Gauss-Jordan elimination.

The elimination procedure works by applying elementary row
operations to our input matrix until the input matrix is reduced to
the identity matrix.
Simultaneously, we apply the same elementary row operations to a
separate identity matrix to produce the inverse matrix.
If this makes no sense, read wikipedia on Gauss-Jordan elimination.

This is not the fastest way to invert matrices, so this is quite
possibly the bottleneck. */
int destructive_invert_matrix(Matrix input, Matrix output) {
	int i;
	int j;
	int r;
	float scalar;
	float shear_needed;
	assert(input.rows == input.cols);
	assert(input.rows == output.rows);
	assert(input.rows == output.cols);

	set_identity_matrix(output);

	/* Convert input to the identity matrix via elementary row operations.
	The ith pass through this loop turns the element at i,i to a 1
	and turns all other elements in column i to a 0. */

	for (i = 0; i < input.rows; ++i) {

		if (input.data[i][i] == 0.0) {
			/* We must swap rows to get a nonzero diagonal element. */

			for (r = i + 1; r < input.rows; ++r) {
				if (input.data[r][i] != 0.0) {
					break;
				}
			}
			if (r == input.rows) {
				/* Every remaining element in this column is zero, so this
				matrix cannot be inverted. */
				return 0;
			}
			swap_rows(input, i, r);
			swap_rows(output, i, r);
		}

		/* Scale this row to ensure a 1 along the diagonal.
		We might need to worry about overflow from a huge scalar here. */
		scalar = 1.0 / input.data[i][i];
		scale_row(input, i, scalar);
		scale_row(output, i, scalar);

		/* Zero out the other elements in this column. */
		for (j = 0; j < input.rows; ++j) {
			if (i == j) {
				continue;
			}
			shear_needed = -input.data[j][i];
			shear_row(input, j, i, shear_needed);
			shear_row(output, j, i, shear_needed);
		}
	}

	return 1;
}

/* 3D矩阵乘以一个1D矩阵，对应位置想乘 */
void scale_mul_matrix(Matrix m, Matrix scalar, bool flag) {
	int i;
	int j;
	// mul
	if (flag) {
		for (i = 0; i < m.rows; ++i)
			for (j = 0; j < m.cols; ++j)
				m.data[i][j] *= scalar.data[0][j];
	}
	else //devide
	{
		for (i = 0; i < m.rows; ++i)
			for (j = 0; j < m.cols; ++j)
				if (abs(scalar.data[2][j]) < 0.000001)
					m.data[i][j] = 0.0;
				else
					m.data[i][j] *= 1.0 / scalar.data[2][j];
	}
	// assert(scalar != 0.0);
	
}

// dian cheng .* 
Matrix dot_multipy(Matrix a, Matrix b) {
	Matrix m;
	assert(a.rows == b.rows);
	assert(a.cols == b.cols);
	int r = a.rows; int c = a.cols;
	m = alloc_matrix(r, c);
	for (int i = 0; i < r; i++){
		for (int j = 0; j < c; j++){
			m.data[i][j] = a.data[i][j] * b.data[i][j];
		}
	}
	return m;
}

/* 矩阵添加一个维度 */
void add_oneD(Matrix source, Matrix destion) {
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < source.cols; j++)
		{
			destion.data[i][j] = source.data[i][j];
		}

	}
	// add 1D
	for (int j = 0; j < source.cols; j++)
	{
		destion.data[3][j] = 1.0;
	}
	return;
}

// resize matrix
void resize_Matrix(DepthImage in, int height, int width, Matrix out) {
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			out.data[0][height*i + j] = in.data[j][i];
		}
	}
	return;
}

/* init depth image Matrix */
DepthImage alloc_depth(int rows, int cols) {
	DepthImage m;
	int i;
	int j;
	m.rows = rows;
	m.cols = cols;
	m.data = (float**)malloc(sizeof(float*) * m.rows);

	for (i = 0; i < m.rows; ++i)
	{
		m.data[i] = (float*)malloc(sizeof(float) * m.cols);
		assert(m.data[i]);
		for (j = 0; j < m.cols; ++j) {
			m.data[i][j] = -1.0;
		}
	}
	return m;
}
/* free the memory */
void free_depth(DepthImage m) {
	int i;
	assert(m.data != NULL);
	for (i = 0; i < m.rows; ++i) {
		free(m.data[i]);
	}
	free(m.data);
}

RGBImage alloc_RGB(int rows, int cols) {
	RGBImage m;
	int i;
	int j;
	m.rows = rows;
	m.cols = cols;
	m.channel = 3;
	m.data = (int**)malloc(sizeof(int*) * m.rows);

	for (i = 0; i < m.rows; ++i)
	{
		m.data[i] = (int*)malloc(sizeof(int) * m.cols * 3);
		assert(m.data[i]);
		memset(m.data[i], 0, sizeof(int)*m.cols * 3);
		/*for (j = 0; j < m.cols; ++j) {
			m.data[i][j] = 0;
		}*/
	}
	return m;
}
/* free the memory */
void free_RGB(RGBImage m) {
	int i;
	assert(m.data != NULL);
	for (i = 0; i < m.rows; ++i) {
		free(m.data[i]);
	}
	free(m.data);
}

#endif // !MAREIXOPERATION
