/**
 * Computes the Euclidean norm of a vector.
 * 
 * @param v Pointer to the input vector
 * @param n Size of the vector
 * @return The Euclidean norm of the vector
 */
double norm(double *v, int n);

/**
 * Computes the dot product between two vectors.
 * 
 * @param v1  Pointer to the first input vector
 * @param n1  Size of the first vector
 * @param v2  Pointer to the second input vector
 * @param n2  Size of the second vector
 * @return    The dot product of the two vectors
 */
double dot(double *v1, int n1, double *v2, int n2);

/**
 * Computes the cross product between two vectors.
 * 
 * @param v1  Pointer to the first input vector
 * @param v2  Pointer to the second input vector
 * @param v3  Pointer to the result
 */
void cross(double *v1, double *v2, double *v3);
