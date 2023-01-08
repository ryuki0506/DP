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
	void time_evolution();
	double *WofE_all_path(int Emax,int pos,int depth);
	double *WofE_all_path(int Emax);
	double *get_FPT();
	double *get_W_Emin();
	
	double calc_pysical_quantity(int calc_mode,bool parcolation,bool Isfixed);
	double calc_Emin(bool parcolation,bool Isfixed);
	double calc_FPT(bool parcolation,bool Isfixed);
	double calc_W_Emin(bool parcolation,bool Isfixed);
	double calc_entropy(bool parcolation,bool Isfixed);

private:
	int field_size;
	double *field;
	double *FPT;
	double *W_Emin;
};


#endif