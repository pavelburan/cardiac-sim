template<class T>
Matrix<T>::Matrix(int r, int c):data(new T[r*c]),dim_r(r),dim_c(c),stride(r),type('r') {
}

template<class T>
Matrix<T>::Matrix(char art, int dim):data(new T[dim*dim]),dim_r(dim),dim_c(dim),stride(dim),type('r') {
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)] = (i==j)?1:0;
	}
}

template<class T>
Matrix<T>::Matrix(T value,int r, int c):data(new T[r*c]),dim_r(r),dim_c(c),stride(r),type('r') {
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)] = value;
	}
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &A):data(new T[A.dim_r*A.dim_c]),dim_r(A.dim_r),dim_c(A.dim_c),stride(A.dim_r),type('r') {
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)] = pAdata[IND(i,j,A.stride)];
	}
}

//Noch zu aendern, nicht konsistent mit normalen c-feldern
template<class T>
Matrix<T>::Matrix(T* feld, int r, int c, char type):data(feld),dim_r(r),dim_c(c),stride(r),type(type) {
	if( type == 'r'){
		data = new T[r*c];
		T* __restrict__ pdata = data;
		const T* __restrict__ pfeld = feld;
		for(int j=0;j<dim_c;j++){
			#pragma omp parallel for
			for(int i=0;i<dim_r;i++)
				pdata[IND(i,j,stride)] = pfeld[IND(i,j,dim_r)];
		}
	}
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &A, int beg_r, int end_r, int beg_c, int end_c):data(A.data+IND(beg_r,beg_c,A.stride)),dim_r(end_r-beg_r+1),
																	  dim_c(end_c-beg_c+1),stride(A.stride),type('l'){
	#ifdef _ERROR_
	if(end_c>=A.dim_c || end_r>=A.dim_r){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Referenced Matrix is to small"<<std::endl;
	   exit(1);
	}
	if(beg_r<0 || beg_c<0){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"beg_r and beg_c have to be greater than zero"<<std::endl;
		exit(1);
	}
	#endif
}

template<class T>
Matrix<T>::~Matrix(){
	if(type=='r')
		delete[] data;
}

template<class T>
Matrix<T> Matrix<T>::operator()(int beg_r, int end_r, int beg_c, int end_c){
	#ifdef _ERROR_
	if(end_c>=dim_c || end_r>=dim_r){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix is to small"<<std::endl;
	   exit(1);
	}
	if(beg_r<0 || beg_c<0){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"beg_r and beg_c have to be greater than zero"<<std::endl;
		exit(1);
	}
	#endif
	return Matrix<T>(*this, beg_r, end_r, beg_c, end_c);
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &A){
	if(A.dim_r*A.dim_c!=dim_r*dim_c){
		#ifdef _ERROR_
		if(type=='l'){
			std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"Matrix dimensions must agree, if the type is 'l'"<<std::endl;
			exit(1);
		}
		#endif
		delete[] data;
		dim_r = A.dim_r;
		dim_c = A.dim_c;
		stride = dim_r;
		data = new T[stride*dim_c];
	}
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)] = pAdata[IND(i,j,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(T d){
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++)
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]=d;
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T> &A){
	#ifdef _ERROR_
	if(dim_r!=A.dim_r || dim_c!=A.dim_c){
	   std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
	   std::cerr<<"Matrix dimensions must agree"<<std::endl;
	   exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]+=pAdata[IND(i,j,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator+=(T d){
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]+=d;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T> &A){
	#ifdef _ERROR_
	if(dim_r!=A.dim_r || dim_c!=A.dim_c){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix dimensions must agree"<<std::endl;
	   exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]-=pAdata[IND(i,j,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator-=(T d){
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]-=d;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T> &A){
	#ifdef _ERROR_
	if(dim_r!=A.dim_r || dim_c!=A.dim_c){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix dimensions must agree"<<std::endl;
	   exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]*=pAdata[IND(i,j,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator*=(T d){
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]*=d;
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T> &A){
	#ifdef _ERROR_
	if(dim_r!=A.dim_r || dim_c!=A.dim_c){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix dimensions must agree"<<std::endl;
	   exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]/=pAdata[IND(i,j,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::operator/=(T d){
	*this *= 1.0/d; 
	return *this;
}

template<class T>
Matrix<T> operator+(const Matrix<T> &A1, const Matrix<T> &A2){
	return Matrix<T>(A1)+=A2;
}

template<class T>
Matrix<T> operator-(const Matrix<T> &A1, const Matrix<T> &A2){
	return Matrix<T>(A1)-=A2;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &A1, const Matrix<T> &A2){
	return Matrix<T>(A1)*=A2;
}

template<class T>
Matrix<T> operator/(const Matrix<T> &A1, const Matrix<T> &A2){
	return Matrix<T>(A1)/=A2;
}

template<class T>
Vector<T> operator*(const Matrix<T> &A,const Vector<T> &v){
	int i,j;
	Vector<T> w(0.0,A.dim_r);
	#ifdef _ERROR_
	if(A.dim_c!=v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pwdata = w.data;
	T* __restrict__ pAdata = A.data;
	T* __restrict__ pvdata = v.data;
	for(i=0;i<A.dim_r;i++){
		T sum = 0.0;
		#pragma omp parallel for reduction(+:sum)
		for(j=0;j<A.dim_c;j++)
			sum+=pAdata[IND(i,j,A.stride)]*pvdata[j];
		pwdata[i] = sum;
	}
	return w;
}

template<class T>
Vector<T> operator*(const Vector<T> &v,const Matrix<T> &A){
	int i,j;
	Vector<T> w(A.dim_c);
	#ifdef _ERROR_
	if(A.dim_r!=v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pwdata = w.data;
	T* __restrict__ pAdata = A.data;
	T* __restrict__ pvdata = v.data;
	for(i=0;i<A.dim_c;i++){
		T sum = 0.0;
		#pragma omp parallel for reduction(+:sum)
		for(j=0;j<A.dim_r;j++)
			sum += pAdata[IND(j,i,A.stride)]*pvdata[j];
		pwdata[i] = sum;
	}
	return w;
}

template<class T>
Matrix<T> operator+(const Matrix<T> &A, T d){
	return Matrix<T>(A)+=d;
}

template<class T>
Matrix<T> operator+(T d, const Matrix<T> &A){
	return Matrix<T>(A)+=d;
}

template<class T>
Matrix<T> operator-(const Matrix<T> &A, T d){
	return Matrix<T>(A)-=d;
}

template<class T>
Matrix<T> operator-(T d, const Matrix<T> &A){
	return (Matrix<T>(A)-=d)*=-1;
}

template<class T>
Matrix<T> operator*(const Matrix<T> &A, T d){
	return Matrix<T>(A)*=d;
}

template<class T>
Matrix<T> operator*(T d, const Matrix<T> &A){
	return Matrix<T>(A)*=d;
}

template<class T>
Matrix<T> operator/(const Matrix<T> &A, T d){
	return Matrix<T>(A)/=d;
}

template<class T>
Matrix<T> Matrix<T>::transpose(){
	#ifdef _ERROR_
	if(type=='l' && dim_r != dim_c){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix dimensions must be square, if the type is 'l'"<<std::endl;
		exit(1);
	}
	#endif
	Matrix<T> A(*this);
	if( dim_r != dim_c)
		stride = dim_c;
	dim_r = A.dim_c;
	dim_c	= A.dim_r;
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)] = pAdata[IND(j,i,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T> Matrix<T>::op(double (*F)(double x)) const{
	Matrix<T> A(dim_r,dim_c);
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pAdata[IND(i,j,A.stride)]=(T)F((double)(pdata[IND(i,j,stride)]));
	}
	return A;
}

template<class T>
Matrix<T> Matrix<T>::pot(T exponent) const {
	Matrix<T> A(dim_r,dim_c);
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int j=0;j<dim_c;j++){
		#pragma omp parallel for
		for(int i=0;i<dim_r;i++)
			pAdata[IND(i,j,A.stride)]=pow(pdata[IND(i,j,stride)], exponent);
	}
	return A;
}

/*Matrix<T> Matrix<T>::pot(int exponent) const {
	Matrix<T> A(1.0,dim_r,dim_c);
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(int i=0; i<abs(exponent); ++i)
		A*=*this;
	if( exponent < 0 )
		for(int j=0;j<dim_c;j++)
			for(int i=0;i<dim_r;i++)
				pAdata[IND(i,j,A.stride)]=1.0/pAdata[IND(i,j,stride)];
	return A;
}*/

template<class T>
Matrix<T>& Matrix<T>::matPro(const Matrix<T> &A){
	T *feld = new T[dim_r*A.dim_c];
	#ifdef _ERROR_
	if(dim_c!=A.dim_r){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix dimensions must agree"<<std::endl;
	    exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	T* __restrict__ pfeld = feld;
	for(int i=0;i<A.dim_c;i++){
		for(int j=0;j<dim_r;j++){
			pfeld[IND(i,j,dim_r)]=0.0;
			#pragma omp parallel for
			for(int k=0;k<dim_c;k++)
				pfeld[IND(i,j,dim_r)]+=pdata[IND(i,k,stride)]*pAdata[IND(k,j,A.stride)];
		}
	}

	delete[] data;
	dim_r=A.dim_r;
	dim_c=A.dim_c;
	stride = dim_r;
	data=feld;
	
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::addNoise(T min, T max){
	#ifdef _ERROR_
	if(max < min){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"max must be bigger than min"<<std::endl;
		exit(1);
	}
	#endif
	
	srand ( time(NULL) );
	T* __restrict__ pdata = data;
	for(int j=0;j<dim_c;j++)
		for(int i=0;i<dim_r;i++)
			pdata[IND(i,j,stride)]+=min+(rand()*(max-min))/RAND_MAX;
	return *this;
}

template<class T>
Matrix<T> matPro(const Matrix<T> &A, const Matrix<T> &B){
	return Matrix<T>(A).matPro(B);
}

template<class T>
void Matrix<T>::print() const{
	int i,j;
	for(i=0;i<dim_r;i++){
		for(j=0;j<dim_c;j++)
			printf("%.16e ",data[IND(i,j,stride)]);
		printf("\n");
	}
}

template<class T>
const Matrix<T>& Matrix<T>::save(const char* filename) const{
	std::ofstream file( filename );
	if( file.is_open() ){
		file.precision(16);
		file.setf(std::ios::scientific);
		for(int i=0;i<dim_r;i++){
			for(int j=0;j<dim_c;j++)
				file<<data[IND(i,j,stride)]<<" ";
			file<<"\n";
		}
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
Matrix<T>& Matrix<T>::load(const char* filename){
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

		if(dim_r*dim_c != rows*cols){
			#ifdef _ERROR_
			if(type=='l'){
				std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
				std::cerr<<"Matrix dimensions must agree, if the type is 'l'"<<std::endl;
				exit(1);
			}
			#endif
			delete[] data;
			dim_r = rows;
			dim_c = cols;
			stride = dim_r;
			data = new T[stride*dim_c];
		}
		//Dateizeiger an den Anfang
		file.clear();
		file.seekg(0L,std::ios::beg);
		for(int i=0;i<dim_r;i++){
			file.getline(zeile, 8192);
			token = strtok(zeile," ");
			for(int j=0;j<dim_c;j++){
				data[IND(i,j,stride)] = atof(token);
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
Matrix<T>& Matrix<T>::insert(int r, int c,const Matrix<T> &A,int ins_r, int ins_c,int anz_r, int anz_c){
	int i,j;
	#ifdef _ERROR_
	if((c+anz_c)>dim_c || (r+anz_r)>dim_r){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix is to small"<<std::endl;
	    exit(1);
	}
	if((ins_c+anz_c)>A.dim_c || (ins_r+anz_r)>A.dim_r){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Insertmatrix is to small"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pAdata = A.data;
	for(i=0;i<anz_r;i++){
		#pragma omp parallel for
		for(j=0;j<anz_c;j++)
			pdata[IND(i+r,j+c,stride)]=pAdata[IND(i+ins_r,j+ins_c,A.stride)];
	}
	return *this;
}

template<class T>
Matrix<T>& Matrix<T>::insert(int r, int c, char direction ,const Vector<T> &v,int ins_x,int length){
	int i;
	#ifdef _ERROR_
	if(direction!='h' && direction!='v'){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Invalid direction"<<std::endl;
	    exit(1);
	}
	if((((c+length)>dim_c || r>=dim_r ) && direction=='h') || (((r+length)>dim_r || c>=dim_c) && direction=='v')){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Matrix is to small"<<std::endl;
	    exit(1);
	}
	if((ins_x+length)>v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"InsertVector is to small"<<std::endl;
		exit(1);
	}
	#endif
	T* __restrict__ pdata = data;
	T* __restrict__ pvdata = v.data;
	if(direction=='v'){
		#pragma omp parallel for
		for(i=0;i<length;i++)
			pdata[IND(i+r,c,stride)]=pvdata[ins_x+i];
	}
	else{
		#pragma omp parallel for
		for(i=0;i<length;i++)
			pdata[IND(r,i+c,stride)]=pvdata[ins_x+i];
	}
	return *this;
}

template<class T>
Vector<T> Matrix<T>::extract(int r, int c, char direction, int length){
	return Vector<T>(length,direction).insert(0,*this,r,c,direction,length);
}

