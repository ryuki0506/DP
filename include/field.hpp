#ifndef __FIELD_H__
#define __FIELD_H__

class Field{
public:
	Field(int len,int E);
	~Field();

	void set_size(int size);
	int get_size();

	void set_potential(double p,int mode);
	double *get_potential();
	void time_evolution(double temp);
	double *get_FPT();
	double *get_W_Emin();

	void set_WofE();
	double *get_WofE();
	void set_Z(double temp);
	double get_Z();
	
	double calc_pysical_quantity(int calc_mode,bool parcolation,bool Isfixed);
	
	double *calc_WofE(int pos,int depth);
	double calc_Emin(bool parcolation,bool Isfixed);
	double calc_FPT(bool parcolation,bool Isfixed);
	double calc_W_Emin(bool parcolation,bool Isfixed);
	double calc_entropy(bool parcolation,bool Isfixed);

private:
	int field_size;
	int Emax;
	double *field;

	double *FPT;
	double *W_Emin;
	
	double Z;
	double *WofE;

};


#endif