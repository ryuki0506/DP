#include <iostream>
#include <cmath>
#include <random>
#include "field.hpp"
using namespace std;

Field::Field(int len):field_size(len){
	field=new double[field_size*field_size];
	//cout<<"Create the field whose length is "<<len<<endl;
}

Field::~Field(){
	delete[] field;
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
			field[field_size*i+j]=dist(engine);
		}
	}
}

double *Field::get_potential(){
	return field;
}

double Field::get_partition_function(){
	double *part=new double[field_size*field_size];
	for (int i = 0; i < field_size; i++)
	{
		for (int j = 0; j < field_size; j++)
		{
			part[field_size*i+j]=0;
		}
	}
	
	for (int i = 0; i < field_size; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i==0 and j==0)
			{
				part[field_size*(i-j)+j]=1;
			}else if(j==0){
				if (field[field_size*(i-j)+j]==0 and field[field_size*(i-j-1)+j]==0)
				{
					part[field_size*(i-j)+j]=part[field_size*(i-j-1)+j];
				}
			}
			else if(j==i){
				if (field[field_size*(i-j)+j]==0 and field[field_size*(i-j)+j-1]==0)
				{
					part[field_size*(i-j)+j]=part[field_size*(i-j)+j-1];
				}
			}else{
				if (field[field_size*(i-j)+j]==0 and (field[field_size*(i-j-1)+j]==0 or field[field_size*(i-j)+j-1]==0))
				{
					part[field_size*(i-j)+j]=part[field_size*(i-j-1)+j]+part[field_size*(i-j)+j-1];
				}
			}
		}
	}

	double sum=0;
	for (int j = 0; j < field_size; j++)
	{
		sum+=part[field_size*(field_size-j-1)+j];
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