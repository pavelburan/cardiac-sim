#ifndef VECTOR_H
#define VECTOR_H

template<class T> class Matrix;
template<class T> class Vector;
template<class T> Vector<T> operator+(const Vector<T> &v1,const Vector<T> &v2);
template<class T> Vector<T> operator-(const Vector<T> &v1,const Vector<T> &v2);
template<class T> Vector<T> operator*(const Vector<T> &v1,const Vector<T> &v2);
template<class T> Vector<T> operator/(const Vector<T> &v1,const Vector<T> &v2);
template<class T> Vector<T> operator*(const Matrix<T> &A,const Vector<T> &v);
template<class T> Vector<T> operator*(const Vector<T> &v,const Matrix<T> &A);
template<class T> Vector<T> operator+(const Vector<T> &v,T d);
template<class T> Vector<T> operator+(T d, const Vector<T> &v);
template<class T> Vector<T> operator-(const Vector<T> &v,T d);
template<class T> Vector<T> operator-(T d, const Vector<T> &v);
template<class T> Vector<T> operator*(const Vector<T> &v,T d);
template<class T> Vector<T> operator*(T d, const Vector<T> &v);
template<class T> Vector<T> operator/(const Vector<T> &v,T d);
template<class T> void VecAddVecMulScal(Vector<T> &erg,const Vector<T> &a,T b, const Vector<T> &c);
template<class T> T scalPro(const Vector<T> &v1, const Vector<T> &v2);


template<class T>
class Vector
{
public:
	//Datenstruktur
	T* data;  //Eintrage des Vektors
	int dim;     //Dimension
	char type;
	//Es gibt folgende Typen von Vektoren
	//'r'->hat eigenen Speicher
	//'l'->hat keinen eigenen Speicher, sondern zeigt auf Speicher von anderem Vector oder Zeile einer Matrix
public:
	//Konstruktoren
	Vector(int dim=0);
	Vector(T value,int dim=0);
	Vector(const Vector<T> &v);
	Vector(T* feld, int dim, char type='r');
	Vector(T a,T b,int N);
	Vector(T a, T dx,T b);
	Vector(const Vector<T> &v, int ins_x, int length);
	Vector(const Matrix<T> &A, int ins_r, int ins_c, int length);
	
	//Zugriffsfunktionen
	int getlen() const {return dim;}
	int size() const { return dim;}
	char gettype() const {return type;}
	void settype(char newtype);
	
	//Destruktor
	~Vector();
	
	//Operatoren
	T &operator()(int i) const {return data[i];};
	T &operator()(int i) {return data[i];};
	T &operator[](int i) const {return data[i];};
	T &operator[](int i) {return data[i];};
	//T &begin() {return data[0];};
	//T &begin() const {return data[0];};
	//T &end() {return data[dim-1];};
	//T &end() const {return data[dim-1];};


	Vector<T> &operator=(const Vector<T> &v);
	Vector<T> &operator+=(const Vector<T> &v);
	Vector<T> &operator-=(const Vector<T> &v);
	Vector<T> &operator*=(const Vector<T> &v);
	Vector<T> &operator/=(const Vector<T> &v);
	Vector<T> &operator=(T d);
	Vector<T> &operator+=(T d);
	Vector<T> &operator-=(T d);
	Vector<T> &operator*=(T d);
	Vector<T> &operator/=(T d);
	friend Vector<T> operator+<T>(const Vector<T> &v1,const Vector<T> &v2);
	friend Vector<T> operator-<T>(const Vector<T> &v1,const Vector<T> &v2);
	friend Vector<T> operator*<T>(const Vector<T> &v1,const Vector<T> &v2);
	friend Vector<T> operator/<T>(const Vector<T> &v1,const Vector<T> &v2);
	friend Vector<T> operator*<T>(const Matrix<T> &A,const Vector<T> &v);
	friend Vector<T> operator*<T>(const Vector<T> &v,const Matrix<T> &A);
	friend Vector<T> operator+<T>(const Vector<T> &v,T d);
	friend Vector<T> operator+<T>(T d, const Vector<T> &v);
	friend Vector<T> operator-<T>(const Vector<T> &v,T d);
	friend Vector<T> operator-<T>(T d, const Vector<T> &v);
	friend Vector<T> operator*<T>(const Vector<T> &v,T d);
	friend Vector<T> operator*<T>(T d, const Vector<T> &v);
	friend Vector<T> operator/<T>(const Vector<T> &v,T d);
	//Saxpy-Operator(erg=b*a+c)
	friend void VecAddVecMulScal<T>(Vector<T> &erg,const Vector<T> &a,T b, const Vector<T> &c);
	//Dateizugriffsmethoden und Ausgabe
	void print(char direction='v') const;  //direction: 'h'->horizontal 'v'->vertical
	const Vector<T>& save(const char* filename, char direction='v') const;
	Vector<T>& load(const char* filename);

	//Bearbeitungsmethoden
	void resize(int new_dim);
	Vector<T> &insert(int xBeg, int xEnd, int step, T value);
	Vector<T> &insert(int x, const Vector<T> &v,int ins_x, int length);
	Vector<T> &insert(int x, const Matrix<T> &A,int ins_r, int ins_c, char direction, int length);
	void linspace(T a,T b, int N);

	//Methoden
	T norm(int p=2) const;
	T norm2() const;
	T norm2square() const;

	Vector<T> op(double (*F)(double x)) const;
	Vector<T> pot(T exponent) const;
	//Vector<T> pot(int exponent) const;
	Vector<T> &addNoise(T min=-1, T max=1);
	friend T scalPro<T>(const Vector<T> &v1, const Vector<T> &v2);

	friend class Matrix<T>;
};

#endif //VECTOR_H
