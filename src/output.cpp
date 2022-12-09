#include <iostream>
#include <vector>
#include "output.hpp"

using namespace std;

void show_field(double *data,int size){
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (j<size-1)
			{
				cout << data[size*i+j]<<",";
			}else{
				cout << data[size*i+j]<<endl;
			}			
		}
	}
};