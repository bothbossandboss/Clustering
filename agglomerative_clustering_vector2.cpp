/**
 * 多次元データを凝集型クラスタリングで分類する。
 * まずは2次元データを扱う。
 * 多次元ガウス分布(共分散=0)から出力された正解のデータをクラスタリングすることを考える。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <assert.h>

#define SINGLE 1
#define COMPLETE 2
#define CENTROID 3
#define VECTOR_DIMENSION 2

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
 * クラスタの重心ベクトルを計算
 */
vector<double> centroid(vector<vector<double> > &cluster){
	vector<double> v = vector<double>(VECTOR_DIMENSION);
	for(int i=0;i<cluster.size();i++){
		for(int j=0;j<VECTOR_DIMENSION;j++){
			v.at(j) += cluster.at(i).at(j);
		}
	}
	for(int j=0;j<VECTOR_DIMENSION;j++){
		v.at(j) /= cluster.size();
	}
	return v;
}
/**
 * 単連結法によるクラスタ間の類似度
 */
double singleLinkSimilarity(vector<vector<double> > &ci, vector<vector<double> > &cj){
	double max = -10.0;
	for(int k=0;k<ci.size();k++){
		for(int l=0;l<cj.size();l++){
			double s = elementSimilarity(ci.at(k), cj.at(l));
			if(s > max) max = s;
		}
	}
	return max;
}
/**
 * 完全連結法によるクラスタ間の類似度
 */
double completeLinkSimilarity(vector<vector<double> > &ci, vector<vector<double> > &cj){
	double min = 10000.0;
	for(int k=0;k<ci.size();k++){
		for(int l=0;l<cj.size();l++){
			double s = elementSimilarity(ci[k], cj[l]);
			if(s < min) min = s;
		}
	}
	return min;
}
/**
 * 重心法によるクラスタ間の類似度
 */
double centroidSimilarity(vector<vector<double> > &ci, vector<vector<double> > &cj){
	return elementSimilarity(centroid(ci), centroid(cj));
}

int main(int argc, char *argv[]){
	return 0;
}