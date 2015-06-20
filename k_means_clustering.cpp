/**
 * 1次元データをk-meansクラスタリングで分類する。
 * k = 2 (2つのクラスタに分類)とする。
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <algorithm>

#define K 2

using namespace std;

double elementSimilarity(double xk, double xl){
	return 1.0 / fabs(xk - xl);
}

bool isEqual(vector<double> &previous, vector<double> &next){
	if(previous.size() != next.size()) return false;
	sort(previous.begin(), previous.end());
	sort(next.begin(), next.end());
	for(int i=0;i<previous.size();i++){
		if(previous.at(i) != next.at(i)) return false;
	}
	return true;
}

void kMeans(vector<double> &d, vector<double> &v1, vector<double> &v2, double *m1, double *m2){
	double sim1, sim2, sum;
	for(int i=0;i<d.size();i++){
		//クラスタ1との類似度
		sim1 = elementSimilarity(d.at(i), *m1);
		//クラスタ2との類似度
		sim2 = elementSimilarity(d.at(i), *m2);
		if(sim1 > sim2){
			//d.at(i)はクラスタ1に所属
			v1.push_back(d.at(i));
		}else{
			//d.at(i)はクラスタ2に所属
			v2.push_back(d.at(i));
		}
	}
	sum = 0.0;
	for(int i=0;i<v1.size();i++){
		printf("%4.2f ", v1.at(i));
		sum += v1.at(i);
	}
	printf("\n");
	*m1 = sum / (double)v1.size();
	sum = 0.0;
	for(int i=0;i<v2.size();i++){
		printf("%4.2f ", v2.at(i));
		sum += v2.at(i);
	}
	printf("\n-----\n");
	*m2 = sum / (double)v2.size();
}

int main(int argc, char *argv[]){
	double m1, m2;
	cout << "m1 = ";
	cin >> m1;
	cout << "m2 = ";
	cin >> m2;
	double x[4] = {0.0, 1.0, 3.0, 5.5};
	double sim1, sim2, sum1, sum2;
	vector<double> d;
	for(int i=0;i<4;i++){
		d.push_back(x[i]);
	}
	vector<double> c1, c2;
	//最初のクラスタリング
	kMeans(d, c1, c2, &m1, &m2);
	while(1){
		//クラスタを更新
		vector<double> v1, v2;
		kMeans(d, v1, v2, &m1, &m2);
		if( (isEqual(c1, v1) && isEqual(c2, v2)) || (isEqual(c1, v2) && isEqual(c2, v1)) ) break;
		c1 = v1;
		c2 = v2;
	}
	return 0;
}