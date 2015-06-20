/**
 * 1次元データを凝集型クラスタリングで分類する。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>

#define SINGLE 1
#define COMPLETE 2
#define CENTROID 3

using namespace std;

double elementSimilarity(double xk, double xl){
	return 1.0 / fabs(xk - xl);
}

double centroid(vector<double> &v){
	double sum = 0.0;
	for(int i=0;i<v.size();i++){
		sum += v.at(i);
	}
	return sum / (double)v.size();
}
/**
 * 単連結法によるクラスタ間の類似度
 */
double singleLinkSimilarity(vector<double> &ci, vector<double> &cj){
	double max = -10.0;
	for(int k=0;k<ci.size();k++){
		for(int l=0;l<cj.size();l++){
			double s = elementSimilarity(ci[k], cj[l]);
			if(s > max) max = s;
		}
	}
	return max;
}
/**
 * 完全連結法によるクラスタ間の類似度
 */
double completeLinkSimilarity(vector<double> &ci, vector<double> &cj){
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
double centroidSimilarity(vector<double> &ci, vector<double> &cj){
	return elementSimilarity(centroid(ci), centroid(cj));
}
/**
 * 凝集型クラスタリング
 */
void agglomerative(vector< vector<double> > &cluster, int *I, int *J, int method){
	double maxS = -10.0;
	for(int i=0;i<cluster.size();i++){
		for(int j=0;j<cluster.size();j++){
			if(i == j)continue;
			double sim;
			switch(method){
				case SINGLE:
					sim  = singleLinkSimilarity(cluster[i], cluster[j]);
					break;
				case COMPLETE:
					sim  = completeLinkSimilarity(cluster[i], cluster[j]);
					break;
				case CENTROID:
					sim = centroidSimilarity(cluster[i], cluster[j]);
					break;
				default:
					break;
			}
			if(sim > maxS ||
				 (method == CENTROID && cluster[*I].size() + cluster[*J].size() > cluster[i].size() + cluster[j].size())){
				maxS = sim;
				*I = i;
				*J = j;
			}
		}
	}
}

int main(int argc, char *argv[]){
	int method;
	cout << "method (single-link = 1, complete-link = 2, centroid = 3): ";
	cin >> method;
	double x[4] = {0.0, 1.0, 3.0, 5.5};
	vector<double> firstC[4];
	vector< vector<double> > cluster;
	for(int i=0;i<4;i++){
		firstC[i].push_back(x[i]);
		cluster.push_back(firstC[i]);
	}
	printf("first\n");
	for(int i=0;i<cluster.size();i++){
		for(int j=0;j<cluster[i].size();j++){
			printf("%f ", cluster[i].at(j));
		}
		printf("\n");
	}
	printf("-----------\n");
	while(cluster.size() >= 2){
		int I, J;
		I = J = -1;
		agglomerative(cluster, &I, &J, method);
		vector<double> tmp = cluster[J];
		cluster[I].insert(cluster[I].end(), tmp.begin(), tmp.end());
		cluster.erase(cluster.begin() + J);
		for(int i=0;i<cluster.size();i++){
			for(int j=0;j<cluster[i].size();j++){
				printf("%4.2f ", cluster[i].at(j));
			}
			printf("\n");
		}
		printf("-----------\n");
	}
	return 0;
}