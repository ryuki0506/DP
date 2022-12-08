#ifndef __FIELD_H__
#define __FIELD_H__

class Field{
public:
	Field(int len);
	~Field();

	void set_size(int size);
	int get_size();

	void set_potential(double p);
	double *get_potential();

	void set_partition_function();
	double *get_partition_function();
	double calc_partition_function();
	double get_growth_rate();

private:
	int field_size;
	double *field;
	double *partition_function;
	double *num_of_least_energy_pathes;
};


#endif