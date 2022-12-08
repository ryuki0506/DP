#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len):field_size(len){
	field=new double[field_size*field_size];
	num_of_least_energy_pathes=new double[field_size*field_size];
	//cout<<"Create the field whose length is "<<len<<endl;
}

Field::~Field(){
	delete[] field;
	delete[] num_of_least_energy_pathes;
}

void Field::set_size(int size){
	field_size=size;
}

int Field::get_size(){
	return field_size;
}

void Field::set_potential(double p){
	random_device seed_gen;
	default_random_engine engine(seed_gen());

	bernoulli_distribution dist(1-p);

	for (int i = 0; i < field_size; i++)
	{
		for (int j = 0; j < field_size; j++){
			field[field_size*i+j]=exp(-dist(engine));
		}
	}
}

double *Field::get_potential(){
	return field;
}

double Field::get_partition_function(){
	double *part=new double[field_size*field_size];//分配関数を作成
	for (int i = 0; i < field_size; i++)//分配関数を初期化
	{
		for (int j = 0; j < field_size; j++)
		{
			part[field_size*i+j]=0;
			num_of_least_energy_pathes[field_size*i+j]=0;
		}
	}
	
	for (int i = 0; i < field_size; i++)//分配関数を漸化式に従って計算
	{
		for (int j = 0; j <= i; j++)
		{
			if (i==0 and j==0)
			{
				part[field_size*(i-j)+j]=field[field_size*(i-j)+j];
				num_of_least_energy_pathes[field_size*(i-j)+j]=1;
			}else if(j==0){
				part[field_size*(i-j)+j]=field[field_size*(i-j)+j]*part[field_size*(i-j-1)+j];
				num_of_least_energy_pathes[field_size*(i-j)+j]=num_of_least_energy_pathes[field_size*(i-j-1)+j];
			}
			else if(j==i){	
				part[field_size*(i-j)+j]=field[field_size*(i-j)+j]*part[field_size*(i-j)+j-1];
				num_of_least_energy_pathes[field_size*(i-j)+j]=num_of_least_energy_pathes[field_size*(i-j)+j-1];
			}else{
				part[field_size*(i-j)+j]=field[field_size*(i-j)+j]*max(part[field_size*(i-j-1)+j],part[field_size*(i-j)+j-1]);
				if (part[field_size*(i-j-1)+j]<=part[field_size*(i-j)+j-1])
				{
					num_of_least_energy_pathes[field_size*(i-j)+j]=num_of_least_energy_pathes[field_size*(i-j)+j-1];
				}else{
					num_of_least_energy_pathes[field_size*(i-j)+j]=num_of_least_energy_pathes[field_size*(i-j-1)+j];
				}
			}
		}
	}

	double sum=0;
	for (int j = 0; j < field_size; j++)
	{
		sum+=num_of_least_energy_pathes[field_size*(field_size-j-1)+j]*part[field_size*(field_size-j-1)+j];
	}
	delete[] part;

	return sum;
}

double Field::get_growth_rate(){
	double Z=get_partition_function();
	if (Z!=0)
	{
		return log(Z)/field_size;
	}else{
		return 0;
	}
}