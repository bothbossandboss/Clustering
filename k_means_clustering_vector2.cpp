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

#define OUTPUT1 "cluster1.dat"
#define OUTPUT2 "cluster2.dat"
#define VECTOR_DIMENSION 2
#define K 2

using namespace std;

/**
 * 要素同士の距離
 * ユークリッド距離を採用。
 */
double elementSimilarity(vector<double> x, vector<double> y){
	assert(x.size() == y.size());
	double sum = 0.0;
	for (int i=0;i<x.size();++i){
		double zi = x.at(i) - y.at(i);
		sum += zi * zi;
	}
	return sqrt(sum); //必要なければ2乗を採用でもいいか。
}
/**
 * 2つのクラスタが等しいか
 */
bool isEqual(vector< vector<double> > &previous, vector< vector<double> > &next){
	if(previous.size() != next.size()) return false;
	sort(previous.begin(), previous.end());
	sort(next.begin(), next.end());
	for(int i=0;i<previous.size();i++){
		for(int j=0;j<VECTOR_DIMENSION;j++){
			if(previous.at(i).at(j) != next.at(i).at(j)) return false;
		}
	}
	return true;
}
/**
 * 主要処理
 * 各データと平均ベクトルの距離を測定し、クラスタリングし直す。
 */
void kMeans(vector< vector<double> > &data, vector< vector<double> > nextCluster[K], vector<double> mu[K]){
	//データを分類
	double sim[K];
	for(int i=0;i<data.size();i++){
		//各クラスタとの類似度
		for(int l=0;l<K;l++){
			sim[l] = elementSimilarity(data.at(i), mu[l]);
		}
		int belongCluster;
		double max = -10.0;
		for(int l=0;l<K;l++){
			if(max < sim[l]){
				max = sim[l];
				belongCluster = l;
			}
		}
		nextCluster[belongCluster].push_back(data.at(i));
	}
	//各クラスタの平均ベクトルの更新
	for(int l=0;l<K;l++){
		double  sum[VECTOR_DIMENSION] = {0.0};
		for(int i=0;i<nextCluster[l].size();i++){
			for(int j=0;j<VECTOR_DIMENSION;j++){
				sum[j] += nextCluster[l].at(i).at(j);
			}
		}
		for(int j=0;j<VECTOR_DIMENSION;j++){
			sum[j] /= (double)nextCluster[l].size();
			mu[l].at(j) = sum[j];
		}
	}
}

int main(int argc, char *argv[]){
	return 0;
}