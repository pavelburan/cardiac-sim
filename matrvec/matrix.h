#ifndef MATRIX_H
#define MATRIX_H
#define IND(i,j,stride) ((j)*(stride)+(i))

template<class T> class Matrix;
template<class T> class Vector;

template<class T> Matrix<T> operator+(const Matrix<T> &A1,const Matrix<T> &A2);
template<class T> Matrix<T> operator-(const Matrix<T> &A1,const Matrix<T> &A2);
template<class T> Matrix<T> operator*(const Matrix<T> &A1,const Matrix<T> &A2);
template<class T> Matrix<T> operator/(const Matrix<T> &A1,const Matrix<T> &A2);
template<class T> Vector<T> operator*(const Matrix<T> &A,const Vector<T> &v);
template<class T> Vector<T> operator*(const Vector<T> &v,const Matrix<T> &A);
template<class T> Matrix<T> operator+(const Matrix<T> &A,T d);
template<class T> Matrix<T> operator+(T d,const Matrix<T> &A);
template<class T> Matrix<T> operator-(const Matrix<T> &A,T d);
template<class T> Matrix<T> operator-(T d,const Matrix<T> &A);
template<class T> Matrix<T> operator*(const Matrix<T> &A,T d);
template<class T> Matrix<T> operator*(T d,const Matrix<T> &A);
template<class T> Matrix<T> operator/(const Matrix<T> &A,T d);
template<class T> Matrix<T> matPro(const Matrix<T> &A, const Matrix<T> &B);


template<class T>
class Matrix
{
public:
	//Datenstruktur
	T* data;     //Eintraege der Matrix
	int dim_r;      //Anzahl der Zeilen
	int dim_c;      //Anzahl der Spalten
	int stride;		//Abstand der Spalten im Speicher
	char type;
	//Es gibt folgende Typen von Matrizen
	//'r'->hat eigenen Speicher
	//'l'->hat keinen eigenen Speicher, sondern zeigt auf Speicher von anderer Matrix

public:
	//Konstruktoren
	Matrix(int r=0, int c=0);
	Matrix(char art, int dim);
	Matrix(T value, int r=0, int c=0);
	Matrix(const Matrix<T> &A);
	Matrix(Matrix<T> &A, int beg_r, int end_r, int beg_c, int end_c);
	Matrix(T* feld, int r, int c, char type = 'r');

	//Zugriffsfunktionen
	int getrows() const {return dim_r;}
	int getcols() const {return dim_c;}
	int getstride() const {return stride;}

	//Destruktor
	~Matrix();

	//Operatoren
	T &operator()(int i, int j) {return data[IND(i,j,dim_r)];} //Zeile i, Spalte j
	T &operator()(int i, int j) const {return data[IND(i,j,dim_r)];} //Zeile i, Spalte j
	Matrix<T> operator()(int beg_r, int end_r, int beg_c, int end_c);
	Matrix<T> &operator=(const Matrix<T> &A);
	Matrix<T> &operator+=(const Matrix<T> &A);
	Matrix<T> &operator-=(const Matrix<T> &A);
	Matrix<T> &operator*=(const Matrix<T> &A);
	Matrix<T> &operator/=(const Matrix<T> &A);
	Matrix<T> &operator=(T d);
	Matrix<T> &operator+=(T d);
	Matrix<T> &operator-=(T d);
	Matrix<T> &operator*=(T d);
	Matrix<T> &operator/=(T d);
	friend Matrix<T> operator+<T>(const Matrix<T> &A1,const Matrix<T> &A2);
	friend Matrix<T> operator-<T>(const Matrix<T> &A1,const Matrix<T> &A2);
	friend Matrix<T> operator*<T>(const Matrix<T> &A1,const Matrix<T> &A2);
	friend Matrix<T> operator/<T>(const Matrix<T> &A1,const Matrix<T> &A2);
	friend Vector<T> operator*<T>(const Matrix<T> &A,const Vector<T> &v);
	friend Vector<T> operator*<T>(const Vector<T> &v,const Matrix<T> &A);
	friend Matrix<T> operator+<T>(const Matrix<T> &A,T d);
	friend Matrix<T> operator+<T>(T d,const Matrix<T> &A);
	friend Matrix<T> operator-<T>(const Matrix<T> &A,T d);
	friend Matrix<T> operator-<T>(T d,const Matrix<T> &A);
	friend Matrix<T> operator*<T>(const Matrix<T> &A,T d);
	friend Matrix<T> operator*<T>(T d,const Matrix<T> &A);
	friend Matrix<T> operator/<T>(const Matrix<T> &A,T d);

	//Dateizugriffsmethoden und Ausgabe
	void print() const;
	const Matrix<T>& save(const char* filename) const;
	Matrix<T>& load(const char* filename);

	//Bearbeitungsmethoden
	Matrix<T> &insert(int r, int c, const Matrix<T> &A, int ins_r, int ins_c,int anz_r, int anz_c);
	Matrix<T> &insert(int r, int c, char direction ,const Vector<T> &v,int ins_x,int length); //direction: 'h'->horizontal 'v'->vertical
	Vector<T> extract(int r, int c, char direction, int length);

	//Methoden
	Matrix<T> transpose();
	Matrix<T> op(double (*F)(double x)) const;
	Matrix<T> pot(T exponent) const;
	//Matrix<T> pot(int exponent) const;
	Matrix<T> &addNoise(T min=-1, T max=1);
	friend Matrix<T> matPro<T>(const Matrix<T> &A, const Matrix<T> &B);
	Matrix<T> &matPro(const Matrix<T> &A);

	friend class Vector<T>;
};

#endif //MATRIX_H

