#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>

using namespace std;
/**
 * 修正コレスキー分解
 * 行列はvector< vector< > > で表現
 * 行列Dは対角行列なので、1次元のvectorで十分。
 * 動作確認済み。
 */
void ModifiedCholeskyDecomposition(vector< vector<double> > &A, vector< vector<double> > &L, vector<double> &D){
	int i, j ,k;
	double sum;
	int N = A.size();
	D.at(0) = A.at(0).at(0);
	L.at(0).at(0) = 1.0;
	for(i=1;i<N;i++){
		//i > jの場合
		for(j=0;j<i;j++){
			sum = 0.0;
			for(k=0;k<j;k++){
				sum += L.at(i).at(k) * D.at(k) * L.at(j).at(k);
			}
			L.at(i).at(j) = (A.at(i).at(j) - sum) / D.at(j);
		}
		//i == jの場合
		L.at(i).at(i) = 1.0;
		sum = 0.0;
		for(k=0;k<i;k++){
			sum += L.at(i).at(k) * L.at(i).at(k) * D.at(k);
		}
		D.at(k) = A.at(i).at(i) - sum;
	}
}

/**
 * L^-1 * y = z <=> L * z = y をガウスの消去法を用いて解く。
 * ここでは、Lは対角成分が全て1の下三角行列とする。
 * 動作確認済み。
 */
vector<double> solveGaussElimination(vector< vector<double> > &L, vector<double> &y){
	vector<double> z = vector<double>(y.size());
	z.at(0) = y.at(0);
	for(int i=1;i<y.size();i++){
		double sum = 0.0;
		for(int k=0;k<i;k++){
			sum += L.at(i).at(k) * z.at(k);
		}
		z.at(i) = y.at(i) - sum;
	}
	return z;
}

double pxGivenC(vector<double> &x, vector<double> &mu, vector< vector<double> > &sigma){
	int N = x.size();
	//分散共分散行列をコレスキー分解する。
	vector< vector<double> > L = vector< vector<double> >(N, vector<double>(N, 0.0));
	vector<double> D = vector<double>(N, 0.0);
	ModifiedCholeskyDecomposition(sigma, L, D);
	//xと平均ベクトルmuとの差を取る。
	vector<double> y = vector<double>(N, 0.0);
	for(int i=0;i<N;i++){
		y.at(i) = x.at(i) - mu.at(i);
	}
	//exp()の中身を計算
	vector<double> z = solveGaussElimination(L, y);
	double sum = 0.0;
	for(int i=0;i<N;i++){
		sum += z.at(i) * z.at(i) / D.at(i);
	}
	//分散共分散行列の行列式を計算
	double detSigma = D.at(0);
	for(int k=1;k<D.size();k++){
		detSigma *= D.at(k);
	}
	double zenhan = 1.0 / ( pow( sqrt(2 * M_PI), N ) * sqrt(detSigma) );
	double kohan = exp(- sum / 2.0);
	return zenhan * kohan;
}

int main(int argc, char *argv[]){
	vector< vector<double> > A(3), L;
	vector<double> D;
//	A = vector< vector<double> >(3, vector<double>(3, 0)); // 2次元vectorの全ての要素を0で初期化。
	L = vector< vector<double> >(3, vector<double>(3, 0));
	D = vector<double>(3, 0);
	A.at(0).push_back(1.0);
	A.at(0).push_back(4.0);
	A.at(0).push_back(5.0);
	A.at(1).push_back(4.0);
	A.at(1).push_back(2.0);
	A.at(1).push_back(4.0);
	A.at(2).push_back(5.0);
	A.at(2).push_back(4.0);
	A.at(2).push_back(3.0);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			printf("%4.2f ", A.at(i).at(j));
		}
		printf("\n");
	}
	printf("\n");
	ModifiedCholeskyDecomposition(A, L, D);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			printf("%4.2f ", L.at(i).at(j));
		}
		printf("\n");
	}
	printf("%4.2f\n", D.at(0));
	printf("  %4.2f\n", D.at(1));
	printf("    %4.2f\n\n", D.at(2));
	vector<double> y, z;
	y.push_back(1.0);
	y.push_back(2.0);
	y.push_back(3.0);
	z = solveGaussElimination(L, y);
	for(int l=0;l<z.size();l++){
		printf("%f\n", z.at(l));
	}
	vector<double> x(2), mu(2);
	vector< vector<double> > sigma = vector< vector<double> >(2, vector<double>(2, 0));
	mu.at(0) = 0.0;
	mu.at(1) = 0.0;
	sigma.at(0).at(0) = 0.9;
	sigma.at(0).at(1) = 0.1;
	sigma.at(1).at(0) = 0.1;
	sigma.at(1).at(1) = 0.9;
	FILE *output;
	if( (output = fopen("p_given_c.txt", "w")) == NULL ){
		perror("open output file");
		return -1;
	}
	for(int i=-30;i<=30;i++){
		for(int j=-30;j<=30;j++){
			x.at(0) = (double)i / 10.0;
			x.at(1) = (double)j / 10.0;
			fprintf(output, "%f %f %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
		}
	}
	x.at(0) = 0.0;
	x.at(1) = 0.0;
	printf("(%4.1f, %4.1f) -> %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
	x.at(0) = 0.9;
	x.at(1) = 0.1;
	printf("(%4.1f, %4.1f) -> %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
	x.at(0) = 0.1;
	x.at(1) = 0.9;
	printf("(%4.1f, %4.1f) -> %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
	x.at(0) = -0.9;
	x.at(1) = -0.1;
	printf("(%4.1f, %4.1f) -> %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
	x.at(0) = -0.1;
	x.at(1) = -0.9;
	printf("(%4.1f, %4.1f) -> %f\n", x.at(0), x.at(1), pxGivenC(x, mu, sigma));
	return 0;
}