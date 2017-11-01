#ifndef DIA_H
#define DIA_H
#include <math.h>
template<class T> class DIA;
template<class T> class Vector;

//template<class T> DIA<T> operator+(const DIA<T> &A1,const DIA<T> &A2);
//template<class T> DIA<T> operator-(const DIA<T> &A1,const DIA<T> &A2);
template<class T> void DIAMulVec(Vector<T> &y,const DIA<T> &A,const Vector<T> &v);
template<class T> Vector<T> operator*(const DIA<T> &A,const Vector<T> &v);
//template<class T> Vector<T> operator*<T>(const Vector<T> &v,const DIA<T> &A);
template<class T> DIA<T> operator+(const DIA<T> &A,T d);
template<class T> DIA<T> operator+(T d,const DIA<T> &A);
template<class T> DIA<T> operator-(const DIA<T> &A,T d);
template<class T> DIA<T> operator-(T d,const DIA<T> &A);
template<class T> DIA<T> operator*(const DIA<T> &A,T d);
template<class T> DIA<T> operator*(T d,const DIA<T> &A);
template<class T> DIA<T> operator/(const DIA<T> &A,T d);

template<class T>
class DIA
{
protected:
	T *data;
	int *offsets;
	int *beg_offsets;
	int *end_offsets;
	int dim_r;
	int dim_c;
	int stride;
	int dim_offsets;
	std::string specType;


public:
	//Konstruktoren
	DIA(int dim_r, int dim_c, int dim_offsets, const T* newData=NULL, const int* newOffsets=NULL, const char* newSpecType="None");
	DIA(const DIA<T> &A);

	//Destruktor
	~DIA();

	//Zugriffsmethoden
	int getrows() const {return dim_r;}
	int getcols() const {return dim_c;}
	int getstride() const {return stride;}
	int getoffset(int offIndex) const;
	int getoffIndex(int offset) const;
	void setData(const T *newData);
	//void setOffsets(const int *newOffsets);
	void setDiagonal(const T* diagonal, int offIndex);
	void addVecToDiagonal(const T* vec, int offIndex);
	void addScalToDiagonal(T d, int offIndex);
	void mulVecToDiagonal(const T* vec, int offIndex);
	void mulScalToDiagonal(T d, int offIndex);



	//Operatoren
	T &operator()(int i, int j); //Zeile i, Spalte j
	T &operator()(int i, int j) const; //Zeile i, Spalte j
	DIA<T> &operator=(const DIA<T> &A);
	DIA<T> &operator=(T d);
	DIA<T> &operator+=(T d);
	DIA<T> &operator-=(T d);
	DIA<T> &operator*=(T d);
	DIA<T> &operator/=(T d);
	
	//friend DIA<T> operator+<T>(const DIA<T> &A1,const DIA<T> &A2);
	//friend DIA<T> operator-<T>(const DIA<T> &A1,const DIA<T> &A2);
	friend void DIAMulVec<T>(Vector<T> &y,const DIA<T> &A,const Vector<T> &v);
	friend Vector<T> operator*<T>(const DIA<T> &A,const Vector<T> &v);
	//friend Vector<T> operator*<T>(const Vector<T> &v,const DIA<T> &A);
	friend DIA<T> operator+<T>(const DIA<T> &A,T d);
	friend DIA<T> operator+<T>(T d,const DIA<T> &A);
	friend DIA<T> operator-<T>(const DIA<T> &A,T d);
	friend DIA<T> operator-<T>(T d,const DIA<T> &A);
	friend DIA<T> operator*<T>(const DIA<T> &A,T d);
	friend DIA<T> operator*<T>(T d,const DIA<T> &A);
	friend DIA<T> operator/<T>(const DIA<T> &A,T d);
	
	//LES_Solver
	//Sequentielle Version von CG
	int CG(Vector<T> &xj, const Vector<T> &b, int MAXIT, T tol) const;

	//Ausgabemethoden
	void print() const;
	void printoffsets() const;
	void printdata() const;
	/*T checksum() const{
		T summe = 0.0;
		for(int i=0;i<dim_offsets*stride;i++)
			summe += data[i]*data[i];
		return summe;
	}
	int checksumoffsets() const{
		int summe = 0;
		for(int i=0;i<dim_offsets;i++)
			summe += offsets[i]*offsets[i];
		return summe;
	}*/
};
#endif /* DIA_H */
