#include <iostream>
#include <cmath>
#include "culc.hpp"

using namespace std;

double limited_average(double *data,int data_len){
	double sum = 0;
	int count=0;
	for(int i=0;i<data_len;i++){
		if (data[i]!= 0)
		{
			sum+=data[i];
			count++;
		}
	}
	if(count!=0){
		return sum/count;
	}else{
		return 0;
	}
};
