template<class T>
DIA<T>::DIA(int dim_r, int dim_c, int dim_offsets, const T* newData, const int* newOffsets, const char* newSpecType):dim_r(dim_r),dim_c(dim_c),stride((dim_r>dim_c)?dim_c:dim_r),
																																							dim_offsets(dim_offsets),
																																							data(new T[dim_offsets*((dim_r>dim_c)?dim_c:dim_r)]),
																																							offsets(new int[dim_offsets]),
																																							beg_offsets(new int[dim_offsets]),
																																							end_offsets(new int[dim_offsets]),
																																							specType(newSpecType){		
	if(newOffsets == NULL){
		for(int k=0;k<dim_offsets;k++)
			offsets[k] = k-0*dim_offsets/2;
	}
	else{
		for(int k=0;k<dim_offsets;k++)
			offsets[k] = newOffsets[k];
	}

	for(int k=0;k<dim_offsets;k++){
		if(offsets[k] < 0)
			beg_offsets[k] = -1*offsets[k];
		else
			beg_offsets[k] = 0.0;
		if(offsets[k] > 0)
			end_offsets[k] = stride - offsets[k];
		else
			end_offsets[k] = stride;
	}
	
	if(specType=="lex3Dppp" || specType=="lex3Dppc" || specType=="lex3Dpcp" || specType=="lex3Dpcc" || specType=="lex3Dcpp" || specType=="lex3Dcpc" || specType=="lex3Dccp" || specType=="lex3Dccc"){
		#ifdef _ERROR_
		if(offsets[6] != -1*offsets[0] || offsets[6] <= offsets[0] || offsets[5] != -1*offsets[1] || offsets[5] <= offsets[1] || offsets[4] != 1 || offsets[2] != -1 || offsets[3] != 0){
			std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"Diamatrix has no "<<specType<<" structure!"<<std::endl;
			exit(1);
		}
		#endif
		if(specType=="lex3Dppp"){
			beg_offsets[0] = 0;
			end_offsets[0] = stride;
			beg_offsets[1] = 0;
			end_offsets[1] = stride;
			beg_offsets[2] = 0;
			end_offsets[2] = stride;
			beg_offsets[4] = 0;
			end_offsets[4] = stride;
			beg_offsets[5] = 0;
			end_offsets[5] = stride;
			beg_offsets[6] = 0;
			end_offsets[6] = stride;
		}
		else if(specType=="lex3Dppc"){
			beg_offsets[1] = 0;
			end_offsets[1] = stride;
			beg_offsets[2] = 0;
			end_offsets[2] = stride;
			beg_offsets[4] = 0;
			end_offsets[4] = stride;
			beg_offsets[5] = 0;
			end_offsets[5] = stride;
		}
		else if(specType=="lex3Dpcp"){
			beg_offsets[0] = 0;
			end_offsets[0] = stride;
			beg_offsets[2] = 0;
			end_offsets[2] = stride;
			beg_offsets[4] = 0;
			end_offsets[4] = stride;
			beg_offsets[6] = 0;
			end_offsets[6] = stride;
		}
		else if(specType=="lex3Dpcc"){
			beg_offsets[2] = 0;
			end_offsets[2] = stride;
			beg_offsets[4] = 0;
			end_offsets[4] = stride;
		}
		else if(specType=="lex3Dcpp"){
			beg_offsets[0] = 0;
			end_offsets[0] = stride;
			beg_offsets[1] = 0;
			end_offsets[1] = stride;
			beg_offsets[5] = 0;
			end_offsets[5] = stride;
			beg_offsets[6] = 0;
			end_offsets[6] = stride;
		}
		else if(specType=="lex3Dcpc"){
			beg_offsets[1] = 0;
			end_offsets[1] = stride;
			beg_offsets[5] = 0;
			end_offsets[5] = stride;
		}
		else if(specType=="lex3Dccp"){
			beg_offsets[0] = 0;
			end_offsets[0] = stride;
			beg_offsets[6] = 0;
			end_offsets[6] = stride;
		}
	}
	
	if(specType=="lex2Dpp" || specType=="lex2Dcp" || specType=="lex2Dpc" || specType=="lex2Dcc"){
		#ifdef _ERROR_
		if(offsets[4] != -1*offsets[0] || offsets[4] <= offsets[0]|| offsets[3] != 1 || offsets[1] != -1 || offsets[2] != 0){
			std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"Diamatrix has no "<<specType<<" structure!"<<std::endl;
			exit(1);
		}
		#endif
		if(specType=="lex2Dpp"){
			for(int k=0;k<dim_offsets;k++){
				beg_offsets[k] = 0;
				end_offsets[k] = stride;
			}
		}
		else if(specType=="lex2Dcp"){
			beg_offsets[0] = 0;
			end_offsets[0] = stride;
			beg_offsets[4] = 0;
			end_offsets[4] = stride;
		}
		else if(specType=="lex2Dpc"){
			beg_offsets[1] = 0;
			end_offsets[1] = stride;
			beg_offsets[3] = 0;
			end_offsets[3] = stride;
		}
	}
	
	if(specType=="lex1Dp" || specType=="lex1Dc"){
		#ifdef _ERROR_
		if(offsets[2] != 1 || offsets[0] != -1 || offsets[1] != 0){
			std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
			std::cerr<<"Diamatrix has no "<<specType<<" structure!"<<std::endl;
			exit(1);
		}
		#endif
		if(specType=="lex1Dp"){
			for(int k=0;k<dim_offsets;k++){
				beg_offsets[k] = 0;
				end_offsets[k] = stride;
			}
		}
	}
	
	
	setData(newData);
	
}

template<class T>
DIA<T>::DIA(const DIA<T> &A):dim_r(A.dim_r),dim_c(A.dim_c),stride(A.stride),dim_offsets(A.dim_offsets),data(new T[A.dim_offsets*A.stride]),offsets(new int[A.dim_offsets]),beg_offsets(new int[A.dim_offsets]),end_offsets(new int[A.dim_offsets]){
	for(int k=0;k<dim_offsets;k++){
		offsets[k] = A.offsets[k];
		beg_offsets[k] = A.beg_offsets[k];
		end_offsets[k] = A.end_offsets[k];
	}

	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim_offsets*stride;i++){
		pdata[i] = A.data[i];
	}
}

template<class T>
DIA<T>::~DIA(){
	delete[] data;
	delete[] offsets;
	delete[] beg_offsets;
	delete[] end_offsets;
}

template<class T>
int DIA<T>::getoffset(int offIndex) const{
	#ifdef _ERROR_
	if(offIndex < 0 || offIndex > dim_offsets){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Offindex="<<offIndex<<" doesn't exist!"<<std::endl;
		exit(1);
	}
	#endif
	return offsets[offIndex];
}

template<class T>
int DIA<T>::getoffIndex(int offset) const{
	int offIndex=2*stride;
	for(int k=0;k<dim_offsets;k++){
		if(offsets[k]==offset)
			offIndex=k;
	}
	#ifdef _ERROR_
	if(offIndex==2*stride){
		std::cerr<<"Error in"<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Offset="<<offset<<" doesn't exist!"<<std::endl;
		exit(1);
	}
	#endif
	return offIndex;
}

template<class T>
void DIA<T>::setData(const T *newData){
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim_offsets*stride;i++){
		pdata[i] = 0.0;
	}
	if(newData != NULL){
		for(int k=0;k<dim_offsets;k++)
			setDiagonal(newData+k*stride,k);
	}
}

/*template<class T>
void DIA<T>::setOffsets(const int *newOffsets){
	if(newOffsets == NULL){
		for(int i=0;i<dim_offsets;i++)
			offsets[i] = i-0*dim_offsets/2;
	}
	else{
		for(int k=0;k<dim_offsets;k++)
			offsets[k] = newOffsets[k];
	}
}*/

template<class T>
void DIA<T>::setDiagonal(const T* diagonal, int offIndex){
	T* __restrict__ pdata = data + offIndex*stride;
	#pragma omp parallel for
	for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
		pdata[i] = diagonal[i];
}

template<class T>
void DIA<T>::addVecToDiagonal(const T* vec, int offIndex){
	T* __restrict__ pdata = data + offIndex*stride;
	#pragma omp parallel for
	for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
		pdata[i] += vec[i];
}

template<class T>
void DIA<T>::addScalToDiagonal(T d, int offIndex){
	T* __restrict__ pdata = data + offIndex*stride;
	#pragma omp parallel for
	for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
		pdata[i] += d;
}

template<class T>
void DIA<T>::mulVecToDiagonal(const T* vec, int offIndex){
	T* __restrict__ pdata = data + offIndex*stride;
	#pragma omp parallel for
	for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
		pdata[i] *= vec[i];
}

template<class T>
void DIA<T>::mulScalToDiagonal(T d, int offIndex){
	T* __restrict__ pdata = data + offIndex*stride;
	#pragma omp parallel for
	for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
		pdata[i] *= d[i];
}

template<class T>
T &DIA<T>::operator()(int i, int j){ //Zeile i, Spalte j
	int offIndex = getoffIndex(j-i);
	return data[offIndex*stride+i];
}

template<class T>
T &DIA<T>::operator()(int i, int j) const{ //Zeile i, Spalte j
	int offIndex = getoffIndex(j-i);
	return data[offIndex*stride+i];
}

template<class T>
DIA<T> &DIA<T>::operator=(const DIA<T> &A){
	delete[] data;
	delete[] offsets;
	delete[] beg_offsets;
	delete[] end_offsets;
	
	dim_r = A.dim_r;
	dim_c = A.dim_c;
	stride = A.stride;
	dim_offsets = A.dim_offsets;
	data = new T[dim_offsets*stride];
	offsets = new int[dim_offsets];
	beg_offsets(new int[dim_offsets]);
	end_offsets(new int[dim_offsets]);
	
	for(int k=0;k<dim_offsets;k++){
		offsets[k] = A.offsets[k];
		beg_offsets[k] = A.beg_offsets[k];
		end_offsets[k] = A.end_offsets[k];
	}
	T* __restrict__ pdata = data;
	#pragma omp parallel for
	for(int i=0;i<dim_offsets*stride;i++){
		pdata[i] = A.data[i];
		//std::cerr<<A.data[i]<<std::endl;
	}
	
	return *this; 	
}

template<class T>
DIA<T> &DIA<T>::operator=(T d){
	for(int offIndex=0;offIndex<dim_offsets;offIndex++){
		T* __restrict__ pdata = data + offIndex*stride;
		#pragma omp parallel for
		for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
			pdata[i] = d;
	}
	
	return *this;
}

template<class T>
DIA<T> &DIA<T>::operator+=(T d){
	for(int offIndex=0;offIndex<dim_offsets;offIndex++){
		T* __restrict__ pdata = data + offIndex*stride;
		#pragma omp parallel for
		for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
			pdata[i] += d;
	}
	
	return *this;
}

template<class T>
DIA<T> &DIA<T>::operator-=(T d){
	for(int offIndex=0;offIndex<dim_offsets;offIndex++){
		T* __restrict__ pdata = data + offIndex*stride;
		#pragma omp parallel for
		for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
			pdata[i] -= d;
	}
	
	return *this;
}

template<class T>
DIA<T> &DIA<T>::operator*=(T d){
	for(int offIndex=0;offIndex<dim_offsets;offIndex++){
		T* __restrict__ pdata = data + offIndex*stride;
		#pragma omp parallel for
		for(int i=beg_offsets[offIndex];i<end_offsets[offIndex];i++)
			pdata[i] *= d;
	}
	
	return *this;
}

template<class T>
DIA<T> &DIA<T>::operator/=(T d){
	*this *= 1.0/d;
	return *this;
}

template<class T>
void DIAMulVec(Vector<T> &y,const DIA<T> &A,const Vector<T> &v){
	#ifdef _ERROR_
	if(y.dim!=A.dim_r || A.dim_c!=v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector DIA Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	int stride = A.stride;
	T* __restrict__ pydata = y.data;
	const T* __restrict__ pAdata = A.data;
	const T* __restrict__ pvdata = v.data;
	if(A.specType == "lex1Dp"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
	
		int nx=A.dim_r;
		int ind = 0;
		pydata[ind] = pAdata[ind]*pvdata[ind+offset0+nx];
		pydata[ind] += pAdata[ind+1*stride]*pvdata[ind+offset1];
		pydata[ind] += pAdata[ind+2*stride]*pvdata[ind+offset2];
		#pragma omp parallel for
		for(int i=1;i<(nx-1);i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0];
			pydata[i] += pAdata[i+1*stride]*pvdata[i+offset1];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
		}
		ind = nx-1;
		pydata[ind] = pAdata[ind]*pvdata[ind+offset0];
		pydata[ind] += pAdata[ind+1*stride]*pvdata[ind+offset1];
		pydata[ind] += pAdata[ind+2*stride]*pvdata[ind+offset2-nx];
	}
	else if(A.specType == "lex1Dc"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
		#pragma omp parallel for
		for(int i=0;i<y.dim;i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0];
			pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
		}
	}
	else if(A.specType == "lex2Dpp" || A.specType == "lex2Dpc" || A.specType == "lex2Dcp"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
		int offset3 = A.offsets[3];
		int offset4 = A.offsets[4];
	
		int nx=offset4;
		int ny=A.dim_r/nx;
		int ind = 0;
		#pragma omp parallel for
		for(int i=nx;i<nx*(ny-1);i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
			pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
		}
		#pragma omp parallel for
		for(int i=0;i<nx;i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0+nx*ny];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
			pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
			pydata[i+nx*(ny-1)] = pAdata[i+nx*(ny-1)]*pvdata[i+nx*(ny-1)+offset0];
			pydata[i+nx*(ny-1)] += pAdata[i+nx*(ny-1)+2*stride]*pvdata[i+nx*(ny-1)+offset2];
			pydata[i+nx*(ny-1)] += pAdata[i+nx*(ny-1)+4*stride]*pvdata[i+nx*(ny-1)+offset4-nx*ny];
		}
		#pragma omp parallel for
		for(int j=0;j<ny;j++){
			ind = j*nx;
			pydata[ind] += pAdata[ind+stride]*pvdata[ind+offset1+nx];
			pydata[ind] += pAdata[ind+3*stride]*pvdata[ind+offset3];
			for(int i=j*nx+1;i<(j+1)*nx-1;i++){
				pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
				pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
			}
			ind = (j+1)*nx - 1;
			pydata[ind] += pAdata[ind+stride]*pvdata[ind+offset1];
			pydata[ind] += pAdata[ind+3*stride]*pvdata[ind+offset3-nx];
		}
	}
	else if(A.specType == "lex2Dcc"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
		int offset3 = A.offsets[3];
		int offset4 = A.offsets[4];
		#pragma omp parallel for
		for(int i=0;i<y.dim;i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0];
			pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
			pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
			pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
		}
	}
	else if(A.specType == "lex3Dppp" || A.specType == "lex3Dppc" || A.specType == "lex3Dpcp" || A.specType == "lex3Dpcc" || A.specType == "lex3Dcpp" || A.specType == "lex3Dcpc" || A.specType == "lex3Dccp"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
		int offset3 = A.offsets[3];
		int offset4 = A.offsets[4];
		int offset5 = A.offsets[5];
		int offset6 = A.offsets[6];
	
		int n = A.dim_r;
		int nxy = offset6;
		int nx = offset5;
		int ny = nxy / nx;
		int nz = n / nxy;
		if(A.specType == "lex3Dppp" || A.specType == "lex3Dpcp" || A.specType == "lex3Dcpp" || A.specType == "lex3Dccp"){
			#pragma omp parallel for
			for(int i=nxy;i<nxy*(nz-1);i++){
				pydata[i] = pAdata[i]*pvdata[i+offset0];
				pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
				pydata[i] += pAdata[i+6*stride]*pvdata[i+offset6];
			}
			#pragma omp parallel for
			for(int i=0;i<nxy;i++){
				pydata[i] = pAdata[i]*pvdata[i+(offset0+n)];
				pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
				pydata[i] += pAdata[i+6*stride]*pvdata[i+offset6];
				pydata[i+nxy*(nz-1)] = pAdata[i+nxy*(nz-1)]*pvdata[i+nxy*(nz-1)+offset0];
				pydata[i+nxy*(nz-1)] += pAdata[i+nxy*(nz-1)+3*stride]*pvdata[i+nxy*(nz-1)+offset3];
				pydata[i+nxy*(nz-1)] += pAdata[i+nxy*(nz-1)+6*stride]*pvdata[i+nxy*(nz-1)+offset6-n];
			}
		}
		else{
			#pragma omp parallel for
			for(int i=0;i<n;i++){
				pydata[i] = pAdata[i]*pvdata[i+offset0];
				pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
				pydata[i] += pAdata[i+6*stride]*pvdata[i+offset6];
			}
		}

		if(A.specType == "lex3Dppp" || A.specType == "lex3Dppc" || A.specType == "lex3Dcpp" || A.specType == "lex3Dcpc"){
			#pragma omp parallel for
			for(int k=0;k<nz;k++){
				for(int i=k*nxy+nx;i<(k+1)*nxy-nx;i++){
					pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
					pydata[i] += pAdata[i+5*stride]*pvdata[i+offset5];
				}
				for(int i=k*nxy;i<k*nxy+nx;i++){
					pydata[i] += pAdata[i+stride]*pvdata[i+offset1+nxy];
					pydata[i] += pAdata[i+5*stride]*pvdata[i+offset5];
					pydata[i+nx*(ny-1)] += pAdata[i+nx*(ny-1)+stride]*pvdata[i+nx*(ny-1)+offset1];
					pydata[i+nx*(ny-1)] += pAdata[i+nx*(ny-1)+5*stride]*pvdata[i+nx*(ny-1)+offset5-nxy];
				}		
			}
		}
		else{
			#pragma omp parallel for
			for(int i=0;i<n;i++){
				pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
				pydata[i] += pAdata[i+5*stride]*pvdata[i+offset5];
			}
		}
		
		if(A.specType == "lex3Dppp" || A.specType == "lex3Dppc" || A.specType == "lex3Dpcp" || A.specType == "lex3Dpcc"){
			#pragma omp parallel for
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					for(int i=k*nxy+j*nx+1;i<k*nxy+(j+1)*nx-1;i++){
						pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
						pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
					}
					int i = k*nxy+j*nx;
					pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2+nx];
					pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
					i = k*nxy+(j+1)*nx-1;
					pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
					pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4-nx];
				}
			}
		}
		else{
			#pragma omp parallel for
			for(int i=0;i<n;i++){
				pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
				pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
			}
		}
	}
	else if(A.specType == "lex3Dccc"){
		int offset0 = A.offsets[0];
		int offset1 = A.offsets[1];
		int offset2 = A.offsets[2];
		int offset3 = A.offsets[3];
		int offset4 = A.offsets[4];
		int offset5 = A.offsets[5];
		int offset6 = A.offsets[6];
		int n = A.dim_r;
		#pragma omp parallel for
		for(int i=0;i<n;i++){
			pydata[i] = pAdata[i]*pvdata[i+offset0];
			pydata[i] += pAdata[i+stride]*pvdata[i+offset1];
			pydata[i] += pAdata[i+2*stride]*pvdata[i+offset2];
			pydata[i] += pAdata[i+3*stride]*pvdata[i+offset3];
			pydata[i] += pAdata[i+4*stride]*pvdata[i+offset4];
			pydata[i] += pAdata[i+5*stride]*pvdata[i+offset5];
			pydata[i] += pAdata[i+6*stride]*pvdata[i+offset6];
		}
	}
	else{
		#pragma omp parallel for
		for(int i=0;i<y.dim;i++)
			pydata[i] = 0.0;
		for(int offIndex=0;offIndex<A.dim_offsets;offIndex++){
			pAdata = A.data+offIndex*A.stride;
			pvdata = v.data+A.offsets[offIndex];
			#pragma omp parallel for
			for(int i=A.beg_offsets[offIndex];i<A.end_offsets[offIndex];i++)
				pydata[i] += pAdata[i]*pvdata[i];
		}
	}
}
/*
template<class T>
void DIAMulVec(Vector<T> &y,const DIA<T> &A,const Vector<T> &v){
	#ifdef _ERROR_
	if(y.dim!=A.dim_r || A.dim_c!=v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"Vector DIA Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	int stride = A.stride;
	T* __restrict__ pydata = y.data;
	const T* __restrict__ pAdata = A.data;
	const T* __restrict__ pvdata = v.data;
	//#pragma omp parallel for
	for(int i=0;i<y.dim;i++)
		pydata[i] = 0.0;
	for(int offIndex=0;offIndex<A.dim_offsets;offIndex++){
		pAdata = A.data+offIndex*A.stride;
		pvdata = v.data+A.offsets[offIndex];
		//#pragma omp parallel for
		for(int i=0;i<y.dim;i++){
			if(i==176){
				std::cerr<<"pydata[i]="<<pydata[i]<<" += "<<"pAdata[i]="<<pAdata[i]<<" * "<<"pvdata[i]="<<pvdata[i]<<std::endl;
			}
			pydata[i] += pAdata[i]*pvdata[i];
			if(i==176){
				std::cerr<<"epydata[i]="<<pydata[i]<<std::endl;
			}
		}
	}	
}
*/

template<class T> 
Vector<T> operator*(const DIA<T> &A,const Vector<T> &v){
	#ifdef _ERROR_
	if(A.dim_c!=v.dim){
		std::cerr<<"Error in "<< __FUNCTION__ << " in " << __FILE__ << " at line " << __LINE__ << std::endl;
		std::cerr<<"DIA Vector dimensions must agree"<<std::endl;
		exit(1);
	}
	#endif
	Vector<T> erg(0.0,A.dim_r);
	DIAMulVec(erg,A,v);
	return erg;
}

template<class T>
DIA<T> operator+(const DIA<T> &A,T d){
	return DIA<T>(A)+=d;
}

template<class T>
DIA<T> operator+(T d,const DIA<T> &A){
	return DIA<T>(A)+=d;
}

template<class T>
DIA<T> operator-(const DIA<T> &A,T d){
	return DIA<T>(A)-=d;
}

template<class T>
DIA<T> operator-(T d,const DIA<T> &A){
	return DIA<T>(A)-=d;
}

template<class T>
DIA<T> operator*(const DIA<T> &A,T d){
	return DIA<T>(A)*=d;
}

template<class T>
DIA<T> operator*(T d,const DIA<T> &A){
	return DIA<T>(A)*=d;
}

template<class T>
DIA<T> operator/(const DIA<T> &A,T d){
	return DIA<T>(A)/=d;
}
#include <stdlib.h>
template<class T>
int DIA<T>::CG(Vector<T> &xj, const Vector<T> &b, int MAXIT, T tol) const{
	Vector<T> r(b-(*this)*xj);
	T rTr = r.norm2square();
	T nr0 = sqrt(rTr);
	Vector<T> pbig(0.0,r.dim+2*dim_c);
	Vector<T> p(pbig,dim_c,r.dim);
	p=r;
	//Vector<T> p(r);
	Vector<T> Ap(xj.dim);
	T gamma, rTr_neu, beta;
	//std::cerr<<"xj="<<xj.norm()<<" b="<<b.norm()<<"A="<<static_cast<const void *>(data)<<std::endl;
	for (int j=0;j<MAXIT;j++){
		//std::cerr<<j<<" "<<rTr<<std::endl;
		//Ap=A*p;
		DIAMulVec(Ap,*this,p);
		gamma = rTr/scalPro(p,Ap);
		//std::cerr<<xj.norm2square()<<" "<<r.norm2square()<<" "<<p.norm2square()<<" "<<Ap.norm2square()<<" "<<gamma<<std::endl;		
		//xj+=gamma*p;
		VecAddVecMulScal(xj,p,gamma,xj);
		//r-=gamma*Ap;
		VecAddVecMulScal(r,Ap,-1*gamma,r);
		rTr_neu = r.norm2square();
		if (sqrt(rTr_neu/nr0) < tol)
			return j; //Abbruch nach Erreichen der tol
		beta = rTr_neu / rTr;
		rTr = rTr_neu;
		//p=r+beta*p;
		VecAddVecMulScal(p,p,beta,r);
	}
	//std::cerr<<"xj="<<xj.norm()<<" b="<<b.norm()<<std::endl;
	return MAXIT;
}

template<class T>
void DIA<T>::print() const{
	int offset,k,zero;
	for(int i=0;i<dim_r;i++){
		for(int j=0;j<dim_c;j++){
			offset=j-i;
			zero=1;
			for(k=0;k<dim_offsets;k++){
				if(offset==offsets[k]){
					std::cout<<std::scientific<<data[i+k*stride]<<" ";
					zero=0;
					break;
				}
			}
			if(zero)
				std::cout<<0.0<<" ";
		}
		std::cout<<std::endl;
	}
}

template<class T>
void DIA<T>::printoffsets() const{
	for(int i=0;i<dim_offsets;i++)
		std::cout<<offsets[i]<<" ";
	std::cout<<std::endl;
}

template<class T>
void DIA<T>::printdata() const{
	for(int i=0;i<dim_r;i++){
		for(int j=0;j<dim_offsets;j++)
			std::cout<<std::scientific<<data[i+j*stride]<<" ";
		std::cout<<std::endl;
	}
}

