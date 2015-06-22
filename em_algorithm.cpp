/**
 * EMアルゴリズムを実装
 * サンプルデータとして、正規分布に独立に従う多次元ベクトルをファイルから読み込む。
 * これをEMアルゴリズムを用いて、適切にクラスタリングする。
 * データはvector<double>で表現。
 * データ列はvector< vector<double> >で表現。
 * 生成したデータ列を標準入力で指定した数の正規分布の混合で表す。
 * その正規分布の平均と標準偏差を推定する。
 * 推定すべきパラメータthetaは、平均と標準偏差とクラスタ確率の組(mu, sigma, P(c))
 * 初期値は以下のように決定する。
 * mu : sigmaの各行を採用。sigmaの大きさ(データの次元)以上のクラスタにはランダムに2つ選んだsigmaの行を足したものを採用。
 * sigma : 全データを一つのガウス分布として捉えた時の分散共分散行列を採用。
 * P(c) : 全てのクラスタに等確率で割り振る。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>

#define EPSILON 0.000001

using namespace std;

/**
 * 2つのvector<double>の距離の2乗
 */
double distance2(vector<double> &x, vector<double> &y){
	assert(x.size() == y.size());
	double sum = 0.0;
	for(int i=0;i<x.size();i++){
		double tmp = x.at(i) - y.at(i);
		sum += tmp * tmp;
	}
	return sum;
}
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
/**
 * クラスタcからデータxが出力される確率P(x|c)
 * 正規分布であると仮定
 */
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
	srand((unsigned int)time(NULL));
	/**
	 * データ準備
	 */
	char inputName[128];
	FILE *input;
	int gmmNum, vectorDimension;
	cout << "input file name : ";
	cin >> inputName;
	cout << "gauss num = ";
	cin >> gmmNum;
	cout << "vector dimension = ";
	cin >> vectorDimension;
	if( (input = fopen(inputName, "r")) == NULL ){
		perror("open input file");
		return -1;
	}
	//ファイルから読み取る。
	vector< vector<double> > data;
	double buf;
	int dim = 0;
	vector<double> v(vectorDimension);
	while( fscanf(input, "%lf", &buf) != EOF ){	//ファイルが終わるまで読み込む。
		if(dim == vectorDimension){
			data.push_back(v);
			dim = 0;
			v = vector<double>(vectorDimension);
		}
		v.at(dim) = buf;
		dim++;
	}
	if(dim == vectorDimension){ //最後の行だけが上記whileループで追加されないので。
		data.push_back(v);
	}
	printf("data size = %d\n", (int)data.size());

	/**
	 * EMアルゴリズムでパラメータを推定
	 */
	//初期値の準備のため、全体を一つのガウス分布として捉える。
	vector<double> uniqMu = vector<double>(vectorDimension, 0.0);
	vector< vector<double> > uniqSigma;
	uniqSigma = vector< vector<double> >(vectorDimension, vector<double>(vectorDimension, 0.0));
	for(int i=0;i<data.size();i++){
		for(int j=0;j<vectorDimension;j++){
			uniqMu.at(j) += data.at(i).at(j);
		}
	}
	for(int j=0;j<vectorDimension;j++){
		uniqMu.at(j) /= data.size();
	}
	for(int i=0;i<data.size();i++){
		for(int k=0;k<vectorDimension;k++){
			for(int l=0;l<vectorDimension;l++){
				double yk = data.at(i).at(k) - uniqMu.at(k);
				double yl = data.at(i).at(l) - uniqMu.at(l);
				uniqSigma.at(k).at(l) += yk * yl;
			}
		}
	}
	for(int k=0;k<vectorDimension;k++){
		for(int l=0;l<vectorDimension;l++){
			uniqSigma.at(k).at(l) /= data.size();
		}
	}
	//パラメータの初期設定
	vector< vector<double> > mu;
	mu = vector< vector<double> >(gmmNum, vector<double>(vectorDimension, 0.0));
	vector< vector< vector<double> > > sigma;
	sigma = vector< vector< vector<double> > >(gmmNum, vector< vector<double> >(vectorDimension, vector<double>(vectorDimension, 0.0)));
	double pc[gmmNum];
	for(int i=0;i<gmmNum;i++){
		if(i > vectorDimension - 1){
			int target1 = rand() % vectorDimension;
			int target2 = rand() % vectorDimension;
			for(int j=0;j<vectorDimension;j++){
				mu.at(i).at(j) = uniqSigma.at(target1).at(j) + uniqSigma.at(target2).at(j);
			}
		}else{
			mu.at(i) = uniqSigma.at(i);
		}
		sigma.at(i) = uniqSigma;
		pc[i] = 1.0 / (double)gmmNum;
	}
	double Q = 0.0;
	int turn = 0;
	//計算開始
	printf("\ncaluculating...\n");
	while(1){
		//vectorであるpcGivenXのサイズはデータ数分。
		vector< vector<double> > pcGivenX = vector< vector<double> >(gmmNum, vector<double>(data.size(), 0.0));
		for(int c=0;c<gmmNum;c++){
			pcGivenX[c] = vector<double>(data.size(), 0.0);
		}
		//Eステップ
		for(int i=0;i<data.size();i++){
			double nowPxGivenC[gmmNum];
			double sum = 0.0;
			for(int c=0;c<gmmNum;c++){
				nowPxGivenC[c] = pxGivenC(data.at(i), mu[c], sigma[c]);
				sum += nowPxGivenC[c];
			}
			//P(c|x;theta')を計算
			for(int c=0;c<gmmNum;c++){
				pcGivenX[c].at(i) = nowPxGivenC[c] / sum;
			}
		}
		//Mステップ
		double Nc[gmmNum]; //sum( P(c|x;theta') )
		double SumAll = 0.0;
		for(int c=0;c<gmmNum;c++){
			Nc[c] = 0.0;
			for(int i=0;i<data.size();i++){
				Nc[c] += pcGivenX[c].at(i);
			}
			SumAll += Nc[c];
		}
		for(int c=0;c<gmmNum;c++){
			vector<double> nextMu = vector<double>(vectorDimension, 0.0);
			vector< vector<double> > nextSigma;
			nextSigma = vector< vector<double> >(vectorDimension, vector<double>(vectorDimension, 0.0));
			//平均ベクトルの更新
			for(int i=0;i<data.size();i++){
				for(int l=0;l<vectorDimension;l++){
					nextMu.at(l) += pcGivenX[c].at(i) * data.at(i).at(l);
				}
			}
			for(int l=0;l<vectorDimension;l++){
				nextMu.at(l) /= Nc[c];
			}
			mu[c] = nextMu;
			//分散共分散行列の更新
			for(int i=0;i<data.size();i++){
				for(int k=0;k<vectorDimension;k++){
					for(int l=0;l<vectorDimension;l++){
						double yk = data.at(i).at(k) - nextMu.at(k);
						double yl = data.at(i).at(l) - nextMu.at(l);
						nextSigma.at(k).at(l) += pcGivenX[c].at(i) * yk * yl;
					}
				}
			}
			for(int k=0;k<vectorDimension;k++){
				for(int l=0;l<vectorDimension;l++){
					nextSigma.at(k).at(l) /= Nc[c];
				}
			}
			sigma[c] = nextSigma;
			//クラスタ出現確率の更新
			pc[c] = Nc[c] / SumAll;
		}
		//Q(theta;theta')->maxを解くことになるが、混合ガウス分布を仮定しているので、それぞれの更新式が導かれる。
		double nextQ = 0.0;
		for(int i=0;i<data.size();i++){
			for(int c=0;c<gmmNum;c++){
				nextQ += pcGivenX[c].at(i) * log( pc[c] * pxGivenC(data.at(i), mu[c], sigma[c]) );
			}
		}
		if(fabs(nextQ - Q) < EPSILON) break;
		//収束しなかったので、Qを更新
		Q = nextQ;
		turn++;
		if(turn % 100 == 0){
			printf("\n<turn = %d>\n", turn);
			for(int c=0;c<gmmNum;c++){
				printf("\ncluster%d\nmu = (", c + 1);
				for(int l=0;l<vectorDimension;l++){
					printf("%4.2f ", mu[c].at(l));
				}
				printf(")\nsigma\n");
				for(int k=0;k<vectorDimension;k++){
					for(int l=0;l<vectorDimension;l++){
						printf("%4.2f ", sigma[c].at(k).at(l));
					}
					printf("\n");
				}
			}
		}
	}
	//結果の出力
	printf("---------\nestimated\n---------");
	for(int c=0;c<gmmNum;c++){
		printf("\n<cluster%d>\nmu = (", c + 1);
		for(int l=0;l<vectorDimension;l++){
			printf("%4.2f ", mu[c].at(l));
		}
		printf(")\nsigma\n");
		for(int k=0;k<vectorDimension;k++){
			for(int l=0;l<vectorDimension;l++){
				printf("%4.2f ", sigma[c].at(k).at(l));
			}
			printf("\n");
		}
	}
	return 0;
}