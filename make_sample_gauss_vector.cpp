/**
 * ガウス分布に従った多次元データを作成する。
 * 各次元の相関は無いとする。(つまり共分散は0)
 * ユーザが指定するのは以下のパラメータ
 * int K : ガウス分布の個数
 * int vectorDimension : データのベクトルの次元数
 * int dataSize : データ数
 * int dataRate[K] : 各ガウス分布の生成率(%)
 * double mu[K][vectorDimension] : 各ガウス分布の平均ベクトル
 * double sigma[K][vectorDimension] : 各ガウス分布における、各次元の分散(共分散は0)
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>

#define OUTPUT "result/gauss.dat"

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

vector<double> generateGaussVector(double mu[], double sigma[], int size){
	vector<double> v;
	for(int i=0;i<size;i++){
		v.push_back(gaussRandCos(mu[i], sigma[i]));
	}
	return v;
}
int main(int argc, char *argv[]){
	/**
	 * パラメータ読み込み
	 */
	int K;
	int vectorDimension;
	int dataSize;
	cout << "gauss num K = ";
	cin >> K;
	cout << "vector dimension = ";
	cin >> vectorDimension;
	cout << "data size = ";
	cin >> dataSize;
	double mu[K][vectorDimension];
	double sigma[K][vectorDimension];
	int dataRate[K];
	for(int i=0;i<K;i++){
		printf("<cluster%d>\n", i);
		for(int j=0;j<vectorDimension;j++){
			cout << "mu[" << j << "] = ";
			cin >> mu[i][j];
		}
		for(int j=0;j<vectorDimension;j++){
			cout << "sigma[" << j << "] = ";
			cin >> sigma[i][j];
		}
		cout << "data rate(%%) = ";
		cin >> dataRate[i];
	}
	int sum = 0;
	for(int i=0;i<K;i++) sum += dataRate[i];
	assert(sum == 100);
	/**
	 * データ列の生成
	 */
	FILE *output;
	if( (output = fopen(OUTPUT, "w")) == NULL ){
		perror("open output file");
		return -1;
	}
	srand((unsigned int)time(NULL));
	vector< vector<double> > data;
	int count[K];
	sum = 0;
	for(int i=0;i<K;i++){
		int num = dataSize * dataRate[i] / 100;
		for(int l=0;l<num;l++){
			data.push_back(generateGaussVector(mu[i], sigma[i], vectorDimension));
		}
		count[i] = num;
		sum += num;
	}
	if(sum < dataSize){
		int remain = dataSize - sum;
		for(int i=0;i<remain;i++){
			int target = rand() % K;
			data.push_back(generateGaussVector(mu[target], sigma[target], vectorDimension));
		}
	}else if(sum > dataSize){
		int extra = sum - dataSize;
		for(int i=0;i<extra;i++){
			int target = rand() % data.size();
			data.erase(data.begin() + target);
		}
	}
	assert(data.size() == dataSize);

	for(int i=0;i<data.size();i++){
		for(int j=0;j<vectorDimension;j++){
			fprintf(output, "%f ", data.at(i).at(j));
		}
		fprintf(output, "\n");
	}
	return 0;
}