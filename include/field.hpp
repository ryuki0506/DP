#ifndef __FIELD_H__
#define __FIELD_H__

class Field{
public:
	Field(int len);
	~Field();

	void set_size(int size);
	int get_size();

	void set_potential(double p,int mode);
	double *get_potential();
	void set_partition_function();
	double *SofE_all_path(int Emax,int pos,int depth);
	double *get_partition_function();
	double *get_num_of_least_energy_pathes();
	double calc_pysical_quantity(int calc_mode,bool parcolation,bool Isfixed);
	double calc_Emin(bool parcolation);
	double get_FPT(bool parcolation,bool Isfixed);
	double get_entropy(bool parcolation,bool Isfixed);

private:
	int field_size;
	double *field;
	double *partition_function;
	double *num_of_least_energy_pathes;
};


#endif