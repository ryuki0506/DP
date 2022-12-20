#include <iostream>
#include <fstream>
#include "output.hpp"

using namespace std;

void show_field(double *data, int size, int max_size,bool Isshow)
{
	if (Isshow && (size==max_size))
	{
		for (int i = 0; i < size; i++)
		{
			for (int j = 0; j < size; j++)
			{
				if (j < size - 1)
				{
					cout << data[size * i + j] << ",";
				}
				else
				{
					cout << data[size * i + j] << endl;
				}
			}
		}
		cout<<endl;
	}
};
void output_settings(int noise_mode,int calc_mode){
	ofstream ofs("../result/settings.txt");
	ofs<<noise_mode<<endl;
	ofs<<calc_mode<<endl;
	ofs.close();
};