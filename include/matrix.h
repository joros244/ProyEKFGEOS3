/**
 * Computes the transpose of a matrix.
 *
 * @param mat Pointer to the input matrix
 * @param res Pointer to the result
 * @param f   Number of rows in mat
 * @param c   Number of columns in mat
 */
void transpose(double **mat, double **res, int f, int c);

/**
 * Computes the multiplication of two matrices.
 *
 * @param mat1 Pointer to the first input matrix
 * @param f1   Number of rows in mat1
 * @param c1   Number of columns in mat1
 * @param mat2 Pointer to the second input matrix
 * @param f2   Number of rows in mat2
 * @param c2   Number of columns in mat2
 * @param res  Pointer to the result
 */
void mult(double **mat1, int f1, int c1, double **mat2, int f2, int c2,
          double **res);

/**
 * Computes the multiplication of three matrices.
 *
 * @param mat1 Pointer to the first input matrix
 * @param f1   Number of rows in mat1
 * @param c1   Number of columns in mat1
 * @param mat2 Pointer to the second input matrix
 * @param f2   Number of rows in mat2
 * @param c2   Number of columns in mat2
 * @param mat3 Pointer to the third input matrix
 * @param f3   Number of rows in mat3
 * @param c3   Number of columns in mat3
 * @param res  Pointer to the result
 */
void mult3(double **mat1, int f1, int c1, double **mat2, int f2, int c2,
           double **mat3, int f3, int c3, double **res);

/**
 * Calculates the inverse of matrix.
 *
 * @param M Input matrix
 * @param r Rows of the matrix
 * @param c Columns of the matrix
 * @param inverted The resulting inverted matrix
 */
void inverMat(double **M, int r, int c, double **inverted);
