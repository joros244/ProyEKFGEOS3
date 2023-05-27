#include "../include/LTC.h"
#include "../include/R_y.h"
#include "../include/R_z.h"
#include "../include/matrix.h"

void LTC(double lon, double lat, double **M){

double **mat1= new double*[3];
double **mat2= new double*[3];

for(int i=0;i<3;i++){
	mat1[i]=new double[3];
	mat2[i]=new double[3];
}

R_y(-1.0*lat,mat1);
R_z(lon,mat2);
mult(mat1,3,3,mat2,3,3,M);

double Aux;
for(int i=0; i<3;i++){
    Aux=M[0][i]; 
    M[0][i]=M[1][i]; 
    M[1][i]=M[2][i]; 
    M[2][i]= Aux;


}
for(int i=0;i<3;i++){
	delete[] mat1[i];
	delete[] mat2[i];
}

delete[] mat1;
delete[] mat2;

}
