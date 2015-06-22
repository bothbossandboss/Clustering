/**
 * 多次元データをk-meansクラスタリングで分類する。
 * 多次元ガウス分布(共分散=0)から出力された正解のデータをクラスタリングすることを考える。
 * データのベクトルの次元とクラスタ数は標準入力から指定する。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <sstream>
#include <iomanip>

#define RESULT_DIR "result/"

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
bool isEqual(vector< vector<double> > previous, vector< vector<double> > next){
	if(previous.size() != next.size()) return false;
	sort(previous.begin(), previous.end());
	sort(next.begin(), next.end());
	for(int i=0;i<previous.size();i++){
		for(int j=0;j<previous.at(i).size();j++){
			if(previous.at(i).at(j) != next.at(i).at(j)) return false;
		}
	}
	return true;
}
/**
 * 主要処理
 * 各データと平均ベクトルの距離を測定し、クラスタリングし直す。
 * nextClusterは空のまま渡す。
 * dataとmuは値が既に入っている。
 */
void kMeans(vector< vector<double> > &data, vector<vector<vector<double> > > &nextCluster, vector<vector<double> > &mu){
	int clusterSize = (int)nextCluster.size();
	int vectorSize = data.at(0).size();
	//データを分類
	double sim[clusterSize];
	for(int i=0;i<data.size();i++){
		//各クラスタとの類似度
		for(int l=0;l<clusterSize;l++){
			sim[l] = elementSimilarity(data.at(i), mu[l]);
		}
		int belongCluster;
		double max = -10.0;
		for(int l=0;l<clusterSize;l++){
			if(max < sim[l]){
				max = sim[l];
				belongCluster = l;
			}
		}
		nextCluster[belongCluster].push_back(data.at(i));
	}
	//各クラスタの平均ベクトルの更新
	for(int l=0;l<clusterSize;l++){
		double  sum[vectorSize];
		for(int j=0;j<vectorSize;j++) sum[j] = 0.0;
		for(int i=0;i<nextCluster[l].size();i++){
			for(int j=0;j<vectorSize;j++){
				sum[j] += nextCluster[l].at(i).at(j);
			}
		}
		for(int j=0;j<vectorSize;j++){
			sum[j] /= (double)nextCluster[l].size();
			mu[l].at(j) = sum[j];
		}
	}
}

int main(int argc, char *argv[]){
	/**
	 * データ準備
	 */
	char inputName[128];
	int vectorDimension, clusterNum;
	FILE *input;
	cout << "input file name : ";
	cin >> inputName;
	cout << "vector dimension = ";
	cin >> vectorDimension;
	cout << "cluster num = ";
	cin >> clusterNum;
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
	fclose(input);

	/**
	 * クラスタリング
	 */
	srand((unsigned int)time(NULL));
	vector< vector< vector<double> > > previousCluster;
	vector< vector<double> > mu;
	//とりあえず初期値はランダムに選択。
	for(int l=0;l<clusterNum;l++){
		int tmp = rand() % (int)data.size();
		mu.push_back(data.at(tmp));
		vector< vector<double> > v;
		previousCluster.push_back(v);
	}
	kMeans(data, previousCluster, mu);
	bool flag = true;
	int turn = 0;
	while(flag){
		if(turn % 10 == 0) printf("turn = %d\n", turn++);
		//クラスタを更新
		vector< vector< vector<double> > > nextCluster;
		for(int l=0;l<clusterNum;l++){
			vector< vector<double> > v;
			nextCluster.push_back(v);
		}
		kMeans(data, nextCluster, mu);
		//クラスタが一致しているか否か
		flag = false;
		for(int i=0;i<clusterNum;i++){
			for(int j=0;j<i;j++){
				if(i == j) continue;
				if(!isEqual(previousCluster[i], nextCluster[j])) flag = true;
			}
		}
		for(int l=0;l<clusterNum;l++){
			previousCluster[l] = nextCluster[l];
		}
	}
	/**
	 * 結果出力
	 */
	for(int l=0;l<clusterNum;l++){
		printf("cluster%d : mu = (", l);
		for(int j=0;j<vectorDimension;j++){
			printf("%f ", mu[l].at(j));
		}
		printf(")\n");
	}

	for(int c=0;c<clusterNum;c++){
		ostringstream oss;
		oss << RESULT_DIR << "cluster" << c << ".dat";
		string str = oss.str();
		char *outputName = (char *)str.c_str();
		FILE *output;
		if( (output = fopen(outputName, "w")) == NULL ){
			perror("open output file");
			return -1;
		}
		for(int i=0;i<previousCluster[c].size();i++){
			for (int j=0;j<vectorDimension;++j){
				fprintf(output, "%f ", previousCluster[c].at(i).at(j));
			}
			fprintf(output, "\n");
		}
		fclose(output);
	}
	return 0;
}