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

#define OUTPUT1 "cluster1.dat"
#define OUTPUT2 "cluster2.dat"
#define SINGLE 1
#define COMPLETE 2
#define CENTROID 3
#define VECTOR_DIMENSION 2
#define END_NUM 2

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
vector<double> centroid(vector< vector<double> > &cluster){
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
double singleLinkSimilarity(vector< vector<double> > &ci, vector< vector<double> > &cj){
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
double completeLinkSimilarity(vector< vector<double> > &ci, vector< vector<double> > &cj){
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
double centroidSimilarity(vector< vector<double> > &ci, vector< vector<double> > &cj){
	return elementSimilarity(centroid(ci), centroid(cj));
}
/**
 * 凝集型クラスタリング
 */
void agglomerative(vector< vector< vector<double> > > &clusters, int *I, int *J, int method){
	double maxS = -10.0;
	for(int i=0;i<clusters.size();i++){
		for(int j=0;j<clusters.size();j++){
			if(i == j)continue;
			double sim;
			switch(method){
				case SINGLE:
					sim  = singleLinkSimilarity(clusters[i], clusters[j]);
					break;
				case COMPLETE:
					sim  = completeLinkSimilarity(clusters[i], clusters[j]);
					break;
				case CENTROID:
					sim = centroidSimilarity(clusters[i], clusters[j]);
					break;
				default:
					break;
			}
			if(sim > maxS ||
				 (method == CENTROID && clusters[*I].size() + clusters[*J].size() > clusters[i].size() + clusters[j].size())){
				maxS = sim;
				*I = i;
				*J = j;
			}
		}
	}
}

int main(int argc, char *argv[]){
	char inputName[128];
	FILE *input, *output1, *output2;
	cout << "input file name : ";
	cin >> inputName;
	int method;
	cout << "method (single-link = 1, complete-link = 2, centroid = 3): ";
	cin >> method;
	if( (input = fopen(inputName, "r")) == NULL ){
		perror("open input file");
		return -1;
	}
	if( ( (output1 = fopen(OUTPUT1, "w")) == NULL )||( (output2 = fopen(OUTPUT2, "w")) == NULL ) ){
		perror("open output file");
		return -1;
	}

	vector< vector<double> > data;
	//ファイルから読み取る。
	vector< vector< vector<double> > > clusters;
	for(int c=0;c<data.size();c++){
		vector< vector<double> > v = vector< vector<double> >(1, vector<double>(VECTOR_DIMENSION, 0.0));
		v.at(0) = data.at(c);
		clusters.push_back(v);
	}
	//クラスタリング開始。クラスタ数が2になったら終了。
	while(clusters.size() > END_NUM){
		int I, J;
		I = J = -1;
		agglomerative(clusters, &I, &J, method);
		vector< vector<double> > tmp = clusters.at(J);
		clusters.at(I).insert(clusters.at(I).end(), tmp.begin(), tmp.end());
		clusters.erase(clusters.begin() + J);
	}
	for(int i=0;i<clusters.at(0).size();i++){
		for (int j=0;j<VECTOR_DIMENSION;++j){
			fprintf(output1, "%f ", clusters.at(0).at(i).at(j));
		}
		printf("\n");
	}
	for(int i=0;i<clusters.at(1).size();i++){
		for (int j=0;j<VECTOR_DIMENSION;++j){
			fprintf(output2, "%f ", clusters.at(1).at(i).at(j));
		}
		printf("\n");
	}
	fclose(input);
	fclose(output1);
	fclose(output2);
	return 0;
}