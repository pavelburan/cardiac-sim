template<class T>
Vector<T>::Vector(int dim):data(new T[dim]),dim(dim),type('r') {
}

template<class T>
Vector<T>::Vector(T value, int dim):data(new T[dim]),dim(dim),type('r') {
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]=value;
}

template<class T>
Vector<T>::Vector(const Vector<T> &v):data(new T[v.dim]),dim(v.dim),type('r') {
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]=pvdata[i];
}

template<class T>
Vector<T>::Vector(T* feld, int dim, char type):data(NULL),dim(dim),type(type){
	#ifdef _ERROR_
	if(type!='r' && type!='l'){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Invalid type for vector"<<std::endl;
		exit(1);
	}
	#endif
	if(type=='l')
		data=feld;
	else{
		data = new T[dim];
		T* __restrict__ pdata = data;
		T* __restrict__ pfeld = feld;
		#pragma omp parallel for
		for(int i=0;i<dim;i++)
			pdata[i]=pfeld[i];
	}
}

template<class T>
Vector<T>::Vector(T a, T b, int N):data(new T[N]),dim(N),type('r'){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0; i<dim;i++)
		pdata[i]=a+(b-a)/(N-1)*i;
}

template<class T>
Vector<T>::Vector(T a, T dx, T b):data(NULL),dim(0),type('r'){
	dim = ceil((b-a)/dx + 1);
	data = new T[dim];
	T xi = a;
	for(int i=0;i<(dim-1);i++){
		data[i] = xi;
		xi += dx;
	}
	data[dim-1] = b;
}

template<class T>
Vector<T>::Vector(const Vector<T> &v, int ins_x, int length):data(v.data+ins_x),dim(length),type('l'){
	#ifdef _ERROR_
	if((ins_x+length)>v.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Insertvector is to small"<<std::endl;
		exit(1);
	}
	#endif
}

template<class T>
Vector<T>::Vector(const Matrix<T> &A, int ins_r, int ins_c, int length):data(A.data+IND(ins_r,ins_c,A.dim_r)),dim(length),type('l'){
	#ifdef _ERROR_
	if((ins_r+length)>A.dim_r || ins_c>=A.dim_c){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Insertmatrix is to small"<<std::endl;
		exit(1);
	}
	#endif
}

template<class T>
Vector<T>::~Vector(){
	if(type=='r')
		delete[] data;
}

template<class T>
void Vector<T>::settype(char newtype){
	#ifdef _ERROR_
	if(type!='l' && type!='r'){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Invalid type for vector"<<std::endl;
	    exit(1);
	}
	#endif

	if(type!='l' && newtype=='l')
		delete[] data;
	if(type=='l' && newtype!='l'){
		T* __restrict__ olddata=data;
		data = new T[dim];
		T* __restrict__ pdata = data;
		for(int i=0; i<dim; i++)
			pdata[i]=olddata[i];
	}
	type=newtype;
}

template<class T>
Vector<T>& Vector<T>::operator=(const Vector<T> &v){
	if(v.dim!=dim){
		#ifdef _ERROR_
		if(type=='l'){
			std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"Vector dimensions must agree, if the type is 'l'"<<std::endl;
			exit(1);
		}
		#endif
		delete[] data;
		data = new T[v.dim];
		dim = v.dim;
	}
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]=pvdata[i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator=(T d){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]=d;
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(const Vector<T> &v){
    #ifdef _ERROR_
	if(dim!=v.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
	   exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]+=pvdata[i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator+=(T d){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]+=d;
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator-=(const Vector<T> &v){
	#ifdef _ERROR_
	if(dim!=v.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]-=pvdata[i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator-=(T d){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]-=d;
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(const Vector<T> &v){
	#ifdef _ERROR_
	if(dim!=v.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]*=pvdata[i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator*=(T d){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]*=d;
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator/=(const Vector<T> &v){
	#ifdef _ERROR_
	if(dim!=v.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]/=pvdata[i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::operator/=(T d){
	*this*=1.0/d;
	/*T* __restrict__ pdata = data;
	for(int i=0;i<dim;i++)
		pdata[i]/=d;*/
	return *this;
}

template<class T>
Vector<T> operator+(const Vector<T> &v1, const Vector<T> &v2){
	Vector<T> erg(v1);
	erg+=v2;
	return erg;
}

template<class T>
Vector<T> operator-(const Vector<T> &v1, const Vector<T> &v2){
	Vector<T> erg(v1);
	erg-=v2;
	return erg;
}

template<class T>
Vector<T> operator*(const Vector<T> &v1, const Vector<T> &v2){
	Vector<T> erg(v1);
	erg*=v2;
	return erg;
}

template<class T>
Vector<T> operator/(const Vector<T> &v1, const Vector<T> &v2){
	Vector<T> erg(v1);
	erg/=v2;
	return erg;
}

template<class T>
Vector<T> operator+(const Vector<T> &v,T d){
	Vector<T> erg(v);
	erg+=d;
	return erg;
}

template<class T>
Vector<T> operator+(T d, const Vector<T> &v){
	Vector<T> erg(v);
	erg+=d;
	return erg;
}

template<class T>
Vector<T> operator-(const Vector<T> &v,T d){
	Vector<T> erg(v);
	erg-=d;
	return erg;
}

template<class T>
Vector<T> operator-(T d, const Vector<T> &v){
	Vector<T> erg(v);
	(erg-=d)*=-1;
	return erg;
}

template<class T>
Vector<T> operator*(const Vector<T> &v,T d){
	Vector<T> erg(v);
	erg*=d;
	return erg;
}

template<class T>
Vector<T> operator*(T d, const Vector<T> &v){
	Vector<T> erg(v);
	erg*=d;
	return erg;
}

template<class T>
Vector<T> operator/(const Vector<T> &v,T d){
	Vector<T> erg(v);
	erg/=d;
	return erg;
}

template<class T>
void VecAddVecMulScal(Vector<T> &erg,const Vector<T> &a,T b, const Vector<T> &c){
	#ifdef _ERROR_
	if(erg.dim!=a.dim || erg.dim != c.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pergdata = erg.data;
	T* __restrict__ padata = a.data;
	T* __restrict__ pcdata = c.data;
	#pragma omp parallel for
	for(int i=0;i<erg.dim;i++)
		pergdata[i]=b*padata[i]+pcdata[i];
}

template<class T>
void Vector<T>::resize(int new_dim){
	T* new_data = new T[new_dim];
	int n = std::max(dim, new_dim);
	T* __restrict__ pdata = data;
	T* __restrict__ pnew_data = new_data;
	#pragma omp parallel for
	for(int i=0;i<n;i++)
		pnew_data[i] = pdata[i];
	delete[] data;
	dim = new_dim;
	data = new_data;
	
	
}


template<class T>
void Vector<T>::print(char direction) const{
	char c=(direction=='v')?'\n':' ';
	for(int i=0;i<dim;i++)
		printf("%e%c",data[i],c);
	printf("\n");
}

template<class T>
const Vector<T>& Vector<T>::save(const char* filename, char direction) const{
	std::ofstream file( filename );
	if( file.is_open() ){
		char c=(direction=='v')?'\n':' ';
		file.precision(16);
		file.setf(std::ios::scientific);
		for(int i=0;i<dim;i++)
			file<<data[i]<<c;
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<filename<<" could not be created"<<std::endl;
		exit(1);
	}
	file.close();
	return *this;
}

template<class T>
Vector<T>& Vector<T>::load(const char* filename){
	std::ifstream file( filename );
	if( file.is_open() ){
		int rows = 1;
		int cols = 0;
		char zeile[8192];
		char *token;
		file.getline(zeile, 8192);
		token = strtok(zeile," ");  // Zerlege zeile in einzelne Elemente (durch Leerzeichen getrennt)
		while (token != 0 )
		{
			++cols;
			token = strtok(0, " ");
		}
		while( !file.eof() ){
			file.getline(zeile, 8192);
			if(strlen(zeile) > 3)
				++rows;
		}
		int new_dim = std::max(rows,cols);

		if(dim != new_dim){
			#ifdef _ERROR_
			if(type=='l'){
				std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
				std::cerr<<"Vector dimensions must agree, if the type is 'l'"<<std::endl;
				exit(1);
			}
			#endif
			delete[] data;
			dim = new_dim;
			data = new T[dim];
		}
		//Dateizeiger an den Anfang
		file.clear();
		file.seekg(0L,std::ios::beg);
		for(int i=0;i<rows;i++){
			file.getline(zeile, 8192);
			token = strtok(zeile," ");
			for(int j=0;j<cols;j++){
				data[i+j] = atof(token);
				token = strtok(0, " ");
			}
		}
	}
	else{
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"File "<<filename<<" could not be opened"<<std::endl;
		exit(1);
	}
	file.close();
	return *this;
}

template<class T>
T Vector<T>::norm(int p) const{
	int i;
	T value=0.0;
	T* __restrict__ pdata = data;
	if(p){
		#pragma omp parallel for reduction(+:value)
		for(i=0;i<dim;i++)
			value+=pow(pdata[i],p);
		return pow(value,1.0/p);
	}
	else{
		for(i=0;i<dim;i++)
			if(data[i]>value)
				value=pdata[i];
		return value;
	}
}

template<class T>
T Vector<T>::norm2() const{
	int i;
	T value=0.0;
	T* __restrict__ pdata = data;
	#pragma omp parallel for reduction(+:value)
	for(i=0;i<dim;i++)
		value += pdata[i]*pdata[i];
	return sqrt(value);
}

template<class T>
T Vector<T>::norm2square() const{
	int i;
	T value=0.0;
	T* __restrict__ pdata = data;
	#pragma omp parallel for reduction(+:value)
	for(i=0;i<dim;i++)
		value += pdata[i]*pdata[i];
	return value;
}

template<class T>
Vector<T> Vector<T>::op(double (*F)(double x)) const{
	Vector<T> v(dim,type);
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pvdata[i]=(T)F((double)pdata[i]);
	return v;
}

template<class T>
Vector<T> Vector<T>::pot(T exponent) const {
	Vector<T> v=Vector<T>(dim,type);
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pvdata[i]=pow(pdata[i],exponent);
	return v;
}

/*
template<class T>
Vector<T> Vector<T>::pot(int exponent) const {
	Vector<T> v(1.0,dim);
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	for(int i=0; i<abs(exponent); ++i)
		v*=*this;
	if( exponent < 0 )
		#pragma omp parallel for
		for(int i=0;i<dim;i++)
			pvdata[i]=1.0/pvdata[i];
	return v;
}*/

template<class T>
Vector<T>& Vector<T>::addNoise(T min, T max){
	#ifdef _ERROR_
	if(max < min){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"max must be bigger than min"<<std::endl;
		exit(1);
	}
	#endif
	
	srand ( time(NULL) );
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim;i++)
		pdata[i]+=min+(rand()*(max-min))/RAND_MAX;
	return *this;
}

template<class T>
T scalPro(const Vector<T> &v1, const Vector<T> &v2){
	#ifdef _ERROR_
	if(v1.dim!=v2.dim){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T value=0.0;
	#pragma omp parallel for reduction(+:value)
	for(int i=0;i<v1.dim;i++)
		value+=v1.data[i]*v2.data[i];
	return value;
}

template<class T>
Vector<T>& Vector<T>::insert(int xBeg, int xEnd, int step, T value){
	#ifdef _ERROR_
	if(xEnd>dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector is to small in"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=xBeg;i<xEnd;i+=step)
		pdata[i] = value;
	return *this;	
}

template<class T>
Vector<T>& Vector<T>::insert(int x, const Vector<T> &v,int ins_x, int length){
	int i;
	#ifdef _ERROR_
	if((x+length)>dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector is to small in"<<std::endl;
		exit(1);
	}
	if((ins_x+length)>v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Insertdatavector is to small"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	#pragma omp parallel for
	for(i=0;i<length;i++)
		pdata[x+i]=pvdata[ins_x+i];
	return *this;
}

template<class T>
Vector<T>& Vector<T>::insert(int x, const Matrix<T> &A,int ins_r, int ins_c, char direction, int length){
	int i;
	#ifdef _ERROR_
	if(direction!='h' && direction!='v'){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Invalid direction in"<<std::endl;
	    exit(1);
	}
	if((x+length)>dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector is to small"<<std::endl;
	    exit(1);
	}
	if((((ins_c+length)>A.dim_c || ins_r>=A.dim_r ) && direction=='h') || (((ins_r+length)>A.dim_r || ins_c>=A.dim_c) && direction=='v')){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Insertmatrix is to small"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	if(type=='v')
		#pragma omp parallel for
		for(i=0;i<length;i++){
			pdata[x+i]=pAdata[IND(ins_r+i,ins_c,A.stride)];
		}
	else
		#pragma omp parallel for
		for(i=0;i<length;i++){
			pdata[x+i]=pAdata[IND(ins_r,ins_c+i,A.stride)];
		}
	return *this;
}

template<class T>
void Vector<T>::linspace(T a,T b, int N){
	#ifdef _ERROR_
	if(type!='r'){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Invalid type for vector"<<std::endl;
		exit(1);
	}
	#endif
	if(N!=dim){
		delete[] data;
		data = new T[N];
	}
	dim=N;
	type='r';
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0; i<dim;i++)
		pdata[i]=a+(b-a)/(N-1)*i;
}

