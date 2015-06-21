/**
 * 多次元データをk-meansクラスタリングで分類する。
 * まずは2次元データを扱う。
 * 多次元ガウス分布(共分散=0)から出力された正解のデータをクラスタリングすることを考える。
 * k = 2 (2つのクラスタに分類)とする。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>

#define OUTPUT1 "result/cluster1_k_means.dat"
#define OUTPUT2 "result/cluster2_k_means.dat"
#define VECTOR_DIMENSION 2
#define K 2

using namespace std;
/**
 * 一様分布(出力確率が全て1)で0〜1の乱数を出力
 */
double uniformRand(){
	return ( (double)rand() + 1.0 ) / ( (double)RAND_MAX + 2.0 );
}
/**
 * 平均mu、標準偏差sigmaの正規分布に従って乱数を出力
 * cos型、sin型両方用意したが、どちらを使っても構わない。
 */
double gaussRandCos(double mu, double sigma){
	double z1 = sqrt( -2.0 * log(uniformRand()) ) * cos( 2.0 * M_PI * uniformRand() );
	return sigma * z1 + mu;
}

double gaussRandSin(double mu, double sigma){
	double z2 = sqrt( -2.0 * log(uniformRand()) ) * sin( 2.0 * M_PI * uniformRand() );
	return sigma * z2 + mu;
}

vector<double> generateGaussVector(double mu, double sigma, int size){
	vector<double> v;
	for(int i=0;i<size;i++){
		v.push_back(gaussRandCos(mu, sigma));
	}
	return v;
}
int main(int argc, char *argv[]){
	/**
	 * データ列を生成
	 */
	double mu1, mu2, sigma1, sigma2;
	cout << "mu1 = ";
	cin >> mu1;
	cout << "mu2 = ";
	cin >> mu2;
	cout << "sigma1 = ";
	cin >> sigma1;
	cout << "sigma2 = ";
	cin >> sigma2;
	printf("---------\nanswer\n---------\n");
	printf("<cluster1>\nmu = (%4.2f %4.2f )", mu1, mu1);
	printf("sigma\n%4.2f 0.0\n0.0 %4.2f\n\n", sigma1, sigma1);
	printf("<cluster2>\nmu = (%4.2f %4.2f )", mu2, mu2);
	printf("sigma\n%4.2f 0.0\n0.0 %4.2f\n", sigma2, sigma2);

	FILE *output;
	if( (output = fopen(OUTPUT, "w")) == NULL ){
		perror("open output file");
		return -1;
	}
	srand((unsigned int)time(NULL));
	vector< vector<double> > data;
	int c1, c2;
	c1 = c2 = 0;
	for(int i=0;i<MAX_POINT;i++){
		int tmp = rand() % 10; //uniformRand()の時とあまり変わらなかったため、こちらを採用している。
		if(tmp < C1_RATE){
			data.push_back(generateGaussVector(mu1, sigma1, VECTOR_DIMENSION));
			c1++;
		}else{
			data.push_back(generateGaussVector(mu2, sigma2, VECTOR_DIMENSION));
			c2++;
		}
	}
	for(int i=0;i<data.size();i++){
		fprintf(output, "%f %f\n", data[i].at(0), data[i].at(1));
	}
	printf("---------\ncluster1 : %4.2f%%, cluster2 : %4.2f%%\n", (double)c1/MAX_POINT, (double)c2/MAX_POINT);
	return 0;
}