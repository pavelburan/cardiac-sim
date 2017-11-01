#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H
#include "system.h"
#include "grids/grid.h"
#include "models/model.h"
#include <time.h>
#include <string.h>
/*omp*/#include <omp.h>
#ifdef PETSC
#include <petscksp.h>
#endif

#define VM_VEC_MODELVARS_VEC_ORDERED //VM_LIN_MODELVARS_VEC_ORDERED VM_LIN_MODELVARS_LIN_ORDERED

#ifdef VM_VEC_MODELVARS_VEC_ORDERED
	#define VM_VEC_ORDERED true
	#define MODELVARS_VEC_ORDERED true
#elif VM_LIN_MODELVARS_VEC_ORDERED
	#define VM_VEC_ORDERED false
	#define MODELVARS_VEC_ORDERED true
#else
	#define VM_VEC_ORDERED false
	#define MODELVARS_VEC_ORDERED false
#endif

void performDiffusionStepFE(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	Vm = rhsVm + CoeffMatrixRhs*Vm;
}

void performDiffusionStepBE(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	CoeffMatrix.CG(Vm, rhsVm, 1000, 0.01);
}

void performDiffusionStepCN(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	rhsVm += CoeffMatrixRhs*Vm;
	CoeffMatrix.CG(Vm, rhsVm, 1000, 0.01);
}

void performDiffusionStepBDF2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	rhsVm += (4.0/3.0)*Vm - (1.0/3.0)*Vmprev;
	CoeffMatrix.CG(Vm, rhsVm, 1000, 0.01);
}

void performDiffusionStepAM2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	rhsVm += CoeffMatrixRhs*((2.0/3.0)*Vm-(1.0/12.0)*Vmprev);
	CoeffMatrix.CG(Vm, rhsVm, 1000, 0.01);
}

void performDiffusionStepMS2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm){
	rhsVm += CoeffMatrixRhs*((4.0/3.0)*Vm+(1.0/3.0)*Vmprev) + Vmprev;
	CoeffMatrix.CG(Vm, rhsVm, 1000, 0.01);
}

#ifdef PETSC

#endif

/*
void timeIntMono_nonadaptive(System &sys){
	double tc;
	const Grid &grid = sys.getGrid();
	const Model &model = grid.getModel(0);
	std::string timeIntSchemeDiffusion = sys.getTimeIntSchemeDiffusion();
	
	bool VmVecOrdered = VM_VEC_ORDERED;
	bool modelVarsVecOrdered = MODELVARS_VEC_ORDERED;
	sys.setStateStructure(VmVecOrdered, modelVarsVecOrdered);
	int n = sys.getn();
	int N = sys.getN();
	int kStride = sys.getStateIndexStride();
	int varStride = sys.getStateVarStride();
	int stateSize = sys.getStateSize();
	int HHBeg = 0;
	int HHEnd = HHBeg + model.getnHH();
	int NLBeg = HHEnd;
	int NLEnd = NLBeg + model.getnNL();
			
	Vector<double> y(0.0,stateSize);
	sys.calcInitialCondition(y.data);
	Vector<double> yprev(y);
	Vector<double> a(0.0,stateSize);
	Vector<double> b(0.0,stateSize);
	Vector<double> Vm(y, 0, n);
	Vector<double> Vmprev(yprev,0,n);
	Vector<double> rhsVm(0.0, n);
	Vector<double> temp(0.0,stateSize);
	
	double tb = sys.gettb(); 
	double te = sys.gette();
	double dt = sys.getdt();
	
	double delta_double = 1e-1;
	
	double *restrict py = y.data;
	double *restrict py_Vm = (VmVecOrdered && modelVarsVecOrdered) ? py+N : py;
	double *restrict py_vars = (VmVecOrdered && modelVarsVecOrdered) ? (py+N+VECLEN) : py+N;
	double *restrict pyprev = yprev.data;
	double *restrict pa = a.data;
	double *restrict pa_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pa+N+VECLEN) : pa+N;
	double *restrict pb = b.data;
	double *restrict pf_Vm = (VmVecOrdered && modelVarsVecOrdered) ? pb+N : pb;
	double *restrict pb_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pb+N+VECLEN) : pb+N;
	double *restrict pf_vars = pb_vars;
	double *restrict pVm = Vm.data;
	double *restrict prhsVm = rhsVm.data;
	double *restrict pdVmdt = pb;
	double *restrict ptemp = temp.data;
	
	bool addVmToRhs;
	double CoeffMatrixAddDiagScal = 0.0;
	double CoeffMatrixRhsAddDiagScal = 0.0;
	double CoeffMatrixFac = 0.0;
	double CoeffMatrixRhsFac = 0.0;
	void (*performDiffusionStep)(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm);
	if(timeIntSchemeDiffusion == "FE"){
		addVmToRhs = false;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepFE;
	}
	else if(timeIntSchemeDiffusion == "BE"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0*dt;
		performDiffusionStep = performDiffusionStepBE;
	}
	else if(timeIntSchemeDiffusion == "CN"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -0.5*dt;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = 0.5*dt;
		performDiffusionStep = performDiffusionStepCN;	
	}
	else if(timeIntSchemeDiffusion == "BDF2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 3.0/2.0;
		CoeffMatrixFac = -1.0*dt;
		performDiffusionStep = performDiffusionStepBDF2;
	}
	else if(timeIntSchemeDiffusion == "AM2"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -5.0/12.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepAM2;	
	}
	else if(timeIntSchemeDiffusion == "MS2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0/3.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepMS2;	
	}
	
	DIA<double> CoeffMatrix = grid.getCoeffMatrixMonoDIA(CoeffMatrixAddDiagScal, CoeffMatrixFac);
	DIA<double> CoeffMatrixRhs = grid.getCoeffMatrixMonoDIA(CoeffMatrixRhsAddDiagScal, CoeffMatrixRhsFac);
		
	//Vec pets_Vm, pets_rhsVm;
	//VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, NULL, &pets_Vm);
	//VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, NULL, &pets_rhsVm);
	//VecPlaceArray(pets_Vm, pVm);
	//VecPlaceArray(pets_rhsVm, prhsVm);
	
	//Mat pets_CoeffMatrix;
	//grid.getCoeffMatrixMonoPetsc(pets_CoeffMatrix, 1.0, -1.0*dt);
	
	//KSP ksp;
	//KSPCreate(PETSC_COMM_WORLD, &ksp);
	//KSPSetOperators(ksp, pets_CoeffMatrix, pets_CoeffMatrix);
	//KSPSetTolerances(ksp, 1.e-6, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT);
	//KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	//KSPSetFromOptions(ksp);
	
	double finished = false;
	double t=tb;
	if(!sys.resume()){
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){
			model.f_Vm_vec(&pVm[k], &py_vars[k*kStride], &pdVmdt[k], t, varStride);
		}
				
		finished = finished || sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
	}
	for(t=tb+dt;t<(te+dt/10)&&!finished;t+=dt){
		std::cout<<"t="<<t<<" te="<<te<<std::endl;
		
		tc=omp_get_wtime();
		yprev = y;
		//double tc2=omp_get_wtime();
		grid.getRhsMono(pVm, prhsVm, t, dt, addVmToRhs);

		//int n = varStride;
		//int pos = N;//(int(N/2.17)/VECLEN)*kStride*VECLEN + int(N/2.17)%VECLEN + N;
		//std::cerr<<py[pos]<<" "<<y[pos+n]<<" "<<y[pos+2*n]<<" "<<y[pos+3*n]<<" "<<y[pos+4*n]<<" "<<y[pos+5*n]<<" "<<y[pos+6*n]<<" "<<y[pos+7*n]<<" "<<prhsVm[0]<<std::endl;
		performDiffusionStep(CoeffMatrix, CoeffMatrixRhs, Vm, Vmprev, rhsVm);
		//std::cerr<<pVm[0]<<" "<<Vmprev[0]<<" "<<prhsVm[0]<<std::endl;
		//KSPSolve(ksp, pets_rhsVm, pets_Vm);
		//int its;
		//KSPGetIterationNumber(ksp,&its);
		//tc2=omp_get_wtime()-tc2;
		//std::cerr<<"ImplStep:"<<tc2*1000<<"ms its:"<<its<<std::endl;
		grid.correctVm(pVm);
		
		yprev = y;
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){		
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				py_Vm[k*kStride+l] = pVm[k+l];
			}			
			
			model.ab_HH_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pa_vars[k*kStride], &pb_vars[k*kStride], t, varStride);
			for(int m=HHBeg;m<HHEnd;m++){
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					double bDiva;
					bDiva = pb_vars[k*kStride+m*varStride+l]/pa_vars[k*kStride+m*varStride+l];
					py_vars[k*kStride+m*varStride+l] = exp(dt*pa_vars[k*kStride+m*varStride+l])*(py_vars[k*kStride+m*varStride+l] + bDiva)-bDiva;
				}
			}
			
			model.f_Vm_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pf_Vm[k*kStride], t, varStride);
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				py_Vm[k*kStride+l] += dt*pf_Vm[k*kStride+l];
			}
			
			model.f_NL_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pf_vars[k*kStride], t, varStride);
			for(int m=NLBeg;m<NLEnd;m++){
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					py_vars[k*kStride+m*varStride+l] += dt*pf_vars[k*kStride+m*varStride+l];
				}
			}
			
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				pVm[k+l] = py_Vm[k*kStride+l];
				pdVmdt[k+l] = pf_Vm[k*kStride+l];
			}
		}
		
		tc=omp_get_wtime()-tc;
		std::cerr<<"ImplStep:"<<tc*1000<<"ms"<<std::endl;

		finished = sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
	}

	t = finished ? t : te;
	sys.obs_end(pyprev, py, pdVmdt, ptemp, t, dt);
	
	//VecResetArray(pets_Vm);
	//VecResetArray(pets_rhsVm);
	//VecDestroy(&pets_Vm);
	//VecDestroy(&pets_rhsVm);
	//MatDestroy(&pets_CoeffMatrix);
	//KSPDestroy(&ksp);
}*/


void timeIntMono_nonadaptive(System &sys){
	double tc;
	const Grid &grid = sys.getGrid();
	const Model &model = grid.getModel(0);
	std::string timeIntSchemeDiffusion = sys.getTimeIntSchemeDiffusion();
	
	bool VmVecOrdered = VM_VEC_ORDERED;
	bool modelVarsVecOrdered = MODELVARS_VEC_ORDERED;
	sys.setStateStructure(VmVecOrdered, modelVarsVecOrdered);
	int n = sys.getn();
	int N = sys.getN();
	int kStride = sys.getStateIndexStride();
	int varStride = sys.getStateVarStride();
	int stateSize = sys.getStateSize();
	int HHBeg = 0;
	int HHEnd = HHBeg + model.getnHH();
	int NLBeg = HHEnd;
	int NLEnd = NLBeg + model.getnNL();
			
	Vector<double> y(0.0,stateSize);
	sys.calcInitialCondition(y.data);
	Vector<double> yprev(y);
	Vector<double> a(0.0,stateSize);
	Vector<double> b(0.0,stateSize);
	Vector<double> Vm(y, 0, n);
	Vector<double> Vmprev(yprev,0,n);
	Vector<double> rhsVm(0.0, n);
	Vector<double> temp(0.0,stateSize);
	
	double tb = sys.gettb(); 
	double te = sys.gette();
	double dt = sys.getdt();
	
	double delta_double = 1e-1;
	
	double *restrict py = y.data;
	double *restrict py_Vm = (VmVecOrdered && modelVarsVecOrdered) ? py+N : py;
	double *restrict py_vars = (VmVecOrdered && modelVarsVecOrdered) ? (py+N+VECLEN) : py+N;
	double *restrict pyprev = yprev.data;
	double *restrict pa = a.data;
	double *restrict pa_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pa+N+VECLEN) : pa+N;
	double *restrict pb = b.data;
	double *restrict pf_Vm = (VmVecOrdered && modelVarsVecOrdered) ? pb+N : pb;
	double *restrict pb_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pb+N+VECLEN) : pb+N;
	double *restrict pf_vars = pb_vars;
	double *restrict pVm = Vm.data;
	double *restrict prhsVm = rhsVm.data;
	double *restrict pdVmdt = pb;
	double *restrict ptemp = temp.data;
	
	bool addVmToRhs;
	double CoeffMatrixAddDiagScal = 0.0;
	double CoeffMatrixRhsAddDiagScal = 0.0;
	double CoeffMatrixFac = 0.0;
	double CoeffMatrixRhsFac = 0.0;
	void (*performDiffusionStep)(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm);
	if(timeIntSchemeDiffusion == "FE"){
		addVmToRhs = false;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepFE;
	}
	else if(timeIntSchemeDiffusion == "BE"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0*dt;
		performDiffusionStep = performDiffusionStepBE;
	}
	else if(timeIntSchemeDiffusion == "CN"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -0.5*dt;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = 0.5*dt;
		performDiffusionStep = performDiffusionStepCN;	
	}
	else if(timeIntSchemeDiffusion == "BDF2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -2.0/3.0*dt;
		performDiffusionStep = performDiffusionStepBDF2;
	}
	else if(timeIntSchemeDiffusion == "AM2"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -5.0/12.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepAM2;	
	}
	else if(timeIntSchemeDiffusion == "MS2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0/3.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepMS2;	
	}
	
	DIA<double> CoeffMatrix = grid.getCoeffMatrixMonoDIA(CoeffMatrixAddDiagScal, CoeffMatrixFac);
	DIA<double> CoeffMatrixRhs = grid.getCoeffMatrixMonoDIA(CoeffMatrixRhsAddDiagScal, CoeffMatrixRhsFac);

	double finished = false;
	double t=tb;
	if(!sys.resume()){
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){
			model.f_Vm_vec(&pVm[k], &py_vars[k*kStride], &pdVmdt[k], t, varStride);
		}
				
		finished = finished || sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
	}
	for(t=tb+dt;t<(te+dt/10)&&!finished;t+=dt){
		std::cout<<"t="<<t<<" te="<<te<<std::endl;
		
		/*omp*/tc=omp_get_wtime();
		yprev = y;
		grid.getRhsMono(pVm, prhsVm, t, dt, addVmToRhs);	
		performDiffusionStep(CoeffMatrix, CoeffMatrixRhs, Vm, Vmprev, rhsVm);
		grid.correctVm(pVm);
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){		
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				py_Vm[k*kStride+l] = pVm[k+l];
			}
			model.f_Vm_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pf_Vm[k*kStride], t, varStride);
			model.ab_HH_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pa_vars[k*kStride], &pb_vars[k*kStride], t, varStride);
			model.f_NL_vec(&py_Vm[k*kStride], &py_vars[k*kStride], &pf_vars[k*kStride], t, varStride);
			
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				py_Vm[k*kStride+l] += dt*pf_Vm[k*kStride+l];
			}
	
			for(int m=HHBeg;m<HHEnd;m++){
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					double bDiva;
					bDiva = pb_vars[k*kStride+m*varStride+l]/pa_vars[k*kStride+m*varStride+l];
					py_vars[k*kStride+m*varStride+l] = exp(dt*pa_vars[k*kStride+m*varStride+l])*(py_vars[k*kStride+m*varStride+l] + bDiva)-bDiva;
				}
			}
						
			for(int m=NLBeg;m<NLEnd;m++){
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					py_vars[k*kStride+m*varStride+l] += dt*pf_vars[k*kStride+m*varStride+l];
				}
			}
			
			#pragma omp simd
			#pragma ivdep
			for(int l=0;l<VECLEN;l++){
				pVm[k+l] = py_Vm[k*kStride+l];
				pdVmdt[k+l] = pf_Vm[k*kStride+l];
			}
		}
		
		finished = sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
		/*omp*/tc=omp_get_wtime()-tc;
		/*omp*/std::cerr<<"ImplStep:"<<tc*1000<<"ms"<<std::endl;
	}

	t = finished ? t : te;
	sys.obs_end(pyprev, py, pdVmdt, ptemp, t, dt);
	
	/*VecResetArray(pets_Vm);
	VecResetArray(pets_rhsVm);
	VecDestroy(&pets_Vm);
	VecDestroy(&pets_rhsVm);
	MatDestroy(&pets_CoeffMatrix);
	KSPDestroy(&ksp);*/
}

void timeIntMono_adaptive(System &sys){
	double tc;
	const Grid &grid = sys.getGrid();
	const Model &model = grid.getModel(0);
	std::string timeIntSchemeDiffusion = sys.getTimeIntSchemeDiffusion();
	
	bool VmVecOrdered = VM_VEC_ORDERED;
	bool modelVarsVecOrdered = MODELVARS_VEC_ORDERED;
	sys.setStateStructure(VmVecOrdered, modelVarsVecOrdered);
	int n = sys.getn();
	int N = sys.getN();
	int kStride = sys.getStateIndexStride();
	int varStride = sys.getStateVarStride();
	int stateSize = sys.getStateSize();
	int HHBeg = 0;
	int HHEnd = HHBeg + model.getnHH();
	int NLBeg = HHEnd;
	int NLEnd = NLBeg + model.getnNL();
			
	Vector<double> y(0.0,stateSize);
	Vector<double> y2(0.0,stateSize);
	Vector<double> y_last(0.0,stateSize);
	sys.calcInitialCondition(y.data);
	Vector<double> yprev(y);
	Vector<double> a(0.0,stateSize);
	Vector<double> b(0.0,stateSize);
	Vector<double> a_last(0.0,stateSize);
	Vector<double> b_last(0.0,stateSize);
	Vector<double> Vm(y, 0, n);
	Vector<double> Vmprev(yprev,0,n);
	Vector<double> rhsVm(0.0, n);
	Vector<double> temp(0.0,stateSize);
	
	double tb = sys.gettb(); 
	double te = sys.gette();
	double dt = sys.getdt();
	
	double adp_tol = 1e-2;
	double delta_double = 1e-2;
	double min_dtAdap = 0.01*dt;
	double max_dtAdap = dt;
	std::vector<double> dt0Adap(N/VECLEN, max_dtAdap);
	
	double *restrict py = y.data;
	double *restrict py_Vm = (VmVecOrdered && modelVarsVecOrdered) ? py+N : py;
	double *restrict py_vars = (VmVecOrdered && modelVarsVecOrdered) ? (py+N+VECLEN) : py+N;
	double *restrict py2 = y2.data;
	double *restrict py_Vm2 = (VmVecOrdered && modelVarsVecOrdered) ? py2+N : py;
	double *restrict py_vars2 = (VmVecOrdered && modelVarsVecOrdered) ? (py2+N+VECLEN) : py+N;
	double *restrict py_last = y_last.data;
	double *restrict py_Vm_last = (VmVecOrdered && modelVarsVecOrdered) ? py_last+N : py;
	double *restrict py_vars_last = (VmVecOrdered && modelVarsVecOrdered) ? (py_last+N+VECLEN) : py+N;
	double *restrict pyprev = yprev.data;
	double *restrict pa = a.data;
	double *restrict pa_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pa+N+VECLEN) : pa+N;
	double *restrict pa_last = a_last.data;
	double *restrict pa_vars_last = (VmVecOrdered && modelVarsVecOrdered) ? (pa_last+N+VECLEN) : pa+N;
	double *restrict pb = b.data;
	double *restrict pf_Vm = (VmVecOrdered && modelVarsVecOrdered) ? pb+N : pb;
	double *restrict pb_vars = (VmVecOrdered && modelVarsVecOrdered) ? (pb+N+VECLEN) : pb+N;
	double *restrict pf_vars = pb_vars;
	double *restrict pb_last = b_last.data;
	double *restrict pf_Vm_last = (VmVecOrdered && modelVarsVecOrdered) ? pb_last+N : pb;
	double *restrict pb_vars_last = (VmVecOrdered && modelVarsVecOrdered) ? (pb_last+N+VECLEN) : pb+N;
	double *restrict pf_vars_last = pb_vars_last;
	double *restrict pVm = Vm.data;
	double *restrict prhsVm = rhsVm.data;
	double *restrict pdVmdt = pb;
	double *restrict ptemp = temp.data;
	
	bool addVmToRhs;
	double CoeffMatrixAddDiagScal = 0.0;
	double CoeffMatrixRhsAddDiagScal = 0.0;
	double CoeffMatrixFac = 0.0;
	double CoeffMatrixRhsFac = 0.0;
	void (*performDiffusionStep)(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vmprev, Vector<double> &rhsVm);
	if(timeIntSchemeDiffusion == "FE"){
		addVmToRhs = false;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepFE;
	}
	else if(timeIntSchemeDiffusion == "BE"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0*dt;
		performDiffusionStep = performDiffusionStepBE;
	}
	else if(timeIntSchemeDiffusion == "CN"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -0.5*dt;
		CoeffMatrixRhsAddDiagScal = 1.0;
		CoeffMatrixRhsFac = 0.5*dt;
		performDiffusionStep = performDiffusionStepCN;	
	}
	else if(timeIntSchemeDiffusion == "BDF2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -2.0/3.0*dt;
		performDiffusionStep = performDiffusionStepBDF2;
	}
	else if(timeIntSchemeDiffusion == "AM2"){
		addVmToRhs = true;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -5.0/12.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepAM2;	
	}
	else if(timeIntSchemeDiffusion == "MS2"){
		addVmToRhs = false;
		CoeffMatrixAddDiagScal = 1.0;
		CoeffMatrixFac = -1.0/3.0*dt;
		CoeffMatrixRhsAddDiagScal = 0.0;
		CoeffMatrixRhsFac = dt;
		performDiffusionStep = performDiffusionStepMS2;	
	}
	
	DIA<double> CoeffMatrix = grid.getCoeffMatrixMonoDIA(CoeffMatrixAddDiagScal, CoeffMatrixFac);
	DIA<double> CoeffMatrixRhs = grid.getCoeffMatrixMonoDIA(CoeffMatrixRhsAddDiagScal, CoeffMatrixRhsFac);

	
	/*Vec pets_Vm, pets_rhsVm;
	VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, NULL, &pets_Vm);
	VecCreateSeqWithArray(PETSC_COMM_SELF, 1, n, NULL, &pets_rhsVm);
	VecPlaceArray(pets_Vm, pVm);
	VecPlaceArray(pets_rhsVm, prhsVm);
	
	Mat pets_CoeffMatrix;
	grid.getCoeffMatrixMonoPetsc(pets_CoeffMatrix, 1.0, -1.0*dt);
	
	KSP ksp;
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetOperators(ksp, pets_CoeffMatrix, pets_CoeffMatrix);
	KSPSetTolerances(ksp, 1.e-6, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT);
	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
	KSPSetFromOptions(ksp);*/
	
	double finished = false;
	double t=tb;
	if(!sys.resume()){
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){
			model.f_Vm_vec(&pVm[k], &py_vars[k*kStride], &pdVmdt[k], t, varStride);
		}
				
		finished = finished || sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
	}
	for(t=tb+dt;t<(te+dt/10)&&!finished;t+=dt){
		std::cout<<"t="<<t<<" te="<<te<<std::endl;
		
		/*omp*/tc=omp_get_wtime();
		yprev = y;
		//double tc2=omp_get_wtime();
		grid.getRhsMono(pVm, prhsVm, t, dt, addVmToRhs);
		performDiffusionStep(CoeffMatrix, CoeffMatrixRhs, Vm, Vmprev, rhsVm);
		//KSPSolve(ksp, pets_rhsVm, pets_Vm);
		//int its;
		//KSPGetIterationNumber(ksp,&its);
		///*omp*/tc2=omp_get_wtime()-tc2;
		///*omp*/std::cerr<<"ImplStep:"<<tc2*1000<<"ms its:"<<its<<std::endl;
		grid.correctVm(pVm);
		
		yprev = y;
		#pragma omp parallel for
		#pragma ivdep
		for(int k=0;k<N;k+=VECLEN){
			const int ind_k = k*kStride;
			
			#ifdef VM_VEC_ORDERED
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					py_Vm[k*kStride+l] = pVm[k+l];
				}
			#endif
			
			//Calculate first right hand side
			model.f_Vm_vec(&py_Vm[ind_k], &py_vars[ind_k], &pf_Vm[ind_k], t, varStride);
			model.ab_HH_vec(&py_Vm[ind_k], &py_vars[ind_k], &pa_vars[ind_k], &pb_vars[ind_k], t, varStride);
			model.f_NL_vec(&py_Vm[ind_k], &py_vars[ind_k], &pf_vars[ind_k], t, varStride);
			double tAdap = t;
			double dtAdap = dt0Adap[k/VECLEN];
			while(tAdap + Eps::t() < t+dt){
				//second righ hand side of last step is first right hand side of actuel step		
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					py_Vm_last[ind_k+l] = py_Vm[ind_k+l];
					pf_Vm_last[ind_k+l] = pf_Vm[ind_k+l];
				}
				for(int m=HHBeg;m<HHEnd;m++){
					const int ind_km = k*kStride+m*varStride;	
					#pragma omp simd
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						py_vars_last[ind_km+l] = py_vars[ind_km+l];
						pa_vars_last[ind_km+l] = pa_vars[ind_km+l];
						pb_vars_last[ind_km+l] = pb_vars[ind_km+l];
					}
				}				
				for(int m=NLBeg;m<NLEnd;m++){
					const int ind_km = k*kStride+m*varStride;	
					#pragma omp simd
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						py_vars_last[ind_km+l] = py_vars[ind_km+l];
						pf_vars_last[ind_km+l] = pf_vars[ind_km+l];
					}
				}
				
				while(true){
					//Calculate Y1
					#pragma omp simd
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						py_Vm[ind_k+l] = py_Vm_last[ind_k+l] + dtAdap*pf_Vm_last[ind_k+l];
					}	
					for(int m=HHBeg;m<HHEnd;m++){
						const int ind_km = k*kStride+m*varStride;	
						#pragma omp simd
						#pragma ivdep
						for(int l=0;l<VECLEN;l++){
							double bDiva = pb_vars_last[ind_km+l]/pa_vars_last[ind_km+l];
							py_vars[ind_km+l] = (fabs(pa_vars_last[ind_km+l])<Eps::rel()) ? py_vars_last[ind_km+l] + dtAdap*(pa_vars_last[ind_km+l]*py_vars_last[ind_km+l] + pb_vars_last[ind_km+l])
																														: exp(dtAdap*pa_vars_last[ind_km+l])*(py_vars_last[ind_km+l] + bDiva)-bDiva;
						}
					}					
					for(int m=NLBeg;m<NLEnd;m++){
						const int ind_km = k*kStride+m*varStride;
						#pragma omp simd
						#pragma ivdep
						for(int l=0;l<VECLEN;l++){
							py_vars[ind_km+l] = py_vars_last[ind_km+l] + dtAdap*pf_vars_last[ind_km+l];
						}
					}
					
					//Calculate Y2
					model.f_Vm_vec(&py_Vm[ind_k], &py_vars[ind_k], &pf_Vm[ind_k], tAdap+dtAdap, varStride);
					model.ab_HH_vec(&py_Vm[ind_k], &py_vars[ind_k], &pa_vars[ind_k], &pb_vars[ind_k], tAdap+dtAdap, varStride);
					model.f_NL_vec(&py_Vm[ind_k], &py_vars[ind_k], &pf_vars[ind_k], tAdap+dtAdap, varStride);
					
					#pragma omp simd
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						py_Vm2[ind_k+l] = py_Vm_last[ind_k+l] + 0.5*dtAdap*(pf_Vm_last[ind_k+l] + pf_Vm[ind_k+l]);
					}	
					for(int m=HHBeg;m<HHEnd;m++){
						const int ind_km = k*kStride+m*varStride;
						#pragma omp simd
						#pragma ivdep
						for(int l=0;l<VECLEN;l++){
							double a_m = 0.5*(pa_vars_last[ind_km+l]+pa_vars[ind_km+l]);
							double b_m = 0.5*(pb_vars_last[ind_km+l]+pb_vars[ind_km+l]);
							double bDiva = b_m/a_m;
							py_vars2[ind_km+l] = (fabs(a_m)<Eps::rel()) ? py_vars_last[ind_km+l] + dtAdap*0.5*(pa_vars_last[ind_km+l]*py_vars_last[ind_km+l] + pb_vars_last[ind_km+l] + pa_vars[ind_km+l]*py_vars[ind_km+l] + pb_vars[ind_km+l])
																					   : exp(dtAdap*a_m)*(py_vars_last[ind_km+l] + bDiva)-bDiva;
						}
					}					
					for(int m=NLBeg;m<NLEnd;m++){
						const int ind_km = k*kStride+m*varStride;
						#pragma omp simd
						#pragma ivdep
						for(int l=0;l<VECLEN;l++){
							py_vars2[ind_km+l] = py_vars_last[ind_km+l] + 0.5*dtAdap*(pf_vars_last[ind_km+l] + pf_vars[ind_km+l]);
						}
					}
							
					//Calculate maximal error between Y1 and Y2
					double error[VECLEN];
					#pragma omp simd
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						error[l] = (py_Vm2[ind_k+l]>delta_double) ? fabs((py_Vm2[ind_k+l]-py_Vm[ind_k+l]) / py_Vm2[ind_k+l])
																  : fabs((py_Vm2[ind_k+l]-py_Vm[ind_k+l]) / delta_double);
					}
					for(int m=HHBeg;m<NLEnd;m++){
						const int ind_km = k*kStride+m*varStride;
						#pragma omp simd
						#pragma ivdep
						for(int l=0;l<VECLEN;l++){
							double tempError = (py_vars2[ind_km+l]>delta_double) ? fabs((py_vars2[ind_km+l]-py_vars[ind_km+l]) / py_vars[ind_km+l])
																				 : fabs((py_vars2[ind_km+l]-py_vars[ind_km+l]) / delta_double);
							error[l] = (error[l] > tempError) ? error[l] : tempError;
						}
					}
					double maxError = 0.0;
					//#pragma omp simd reduction(maximum:maxError)
					#pragma ivdep
					for(int l=0;l<VECLEN;l++){
						maxError = (maxError > error[l]) ? maxError : error[l];
					}
					
					//calculate new dtAdap
					if(maxError < 2.0*adp_tol || dtAdap - min_dtAdap < Eps::t()){
						double new_dtAdap = dtAdap*sqrt(adp_tol/maxError);
						tAdap += dtAdap;
						new_dtAdap = (new_dtAdap > min_dtAdap) ? new_dtAdap : min_dtAdap;
						if(tAdap + new_dtAdap > t+dt+Eps::t()){
							dt0Adap[k/VECLEN] = (new_dtAdap < max_dtAdap) ? new_dtAdap : max_dtAdap;
							dtAdap = t+dt-tAdap;
						}
						else
							dtAdap = new_dtAdap;
						break;
					}
					else{
						dtAdap = (dtAdap > 0.5*dtAdap) ? 0.5*dtAdap : min_dtAdap;
					}
				}
				//if(k==0) std::cerr<<"End tAdap="<<tAdap<<" t+dt="<<t+dt<<" dtAdap="<<dtAdap<<std::endl;
			}
			
			#ifdef VM_VEC_ORDERED
				#pragma omp simd
				#pragma ivdep
				for(int l=0;l<VECLEN;l++){
					pVm[k+l] = py_Vm[k*kStride+l];
					pdVmdt[k+l] = pf_Vm[k*kStride+l];
				}
			#endif	
		}
		/*omp*/tc=omp_get_wtime()-tc;
		/*omp*/std::cerr<<"ImplStep:"<<tc*1000<<"ms"<<std::endl;

		finished = sys.obs(pyprev, py, pdVmdt, ptemp, t, dt);
	}
	t = finished ? t : te;
	sys.obs_end(pyprev, py, pdVmdt, ptemp, t, dt);
	
	/*VecResetArray(pets_Vm);
	VecResetArray(pets_rhsVm);
	VecDestroy(&pets_Vm);
	VecDestroy(&pets_rhsVm);
	MatDestroy(&pets_CoeffMatrix);
	KSPDestroy(&ksp);*/
}

#endif //TIMEINTEGRATION_H

