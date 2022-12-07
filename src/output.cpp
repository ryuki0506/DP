#include <iostream>
#include <vector>
#include "output.hpp"

using namespace std;

void output_GR(double *data,int data_row,int data_column){
	for(int i=0;i<data_row;i++){
		for(int j=0;j<data_column;j++){
			if (j<data_column-1)
			{
				cout << data[data_row*i+j]<<", ";
			}else{
				cout << data[data_row*i+j]<<endl;
			}
		}
	}
};
