#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H

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

void performDiffusionStepFE(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void performDiffusionStepBE(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void performDiffusionStepCN(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void performDiffusionStepBDF2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void performDiffusionStepAM2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void performDiffusionStepMS2(DIA<double> &CoeffMatrix, DIA<double> &CoeffMatrixRhs, Vector<double> &Vm, Vector<double> &Vm_prev, Vector<double> &rhsVm);

void timeInt_nonadaptive(const Model &model, double *py, double *pa, double *pb, double t, double dt);

void timeIntMono_nonadaptive(System &sys);

void timeIntMono_adaptive(System &sys);

#endif //TIMEINTEGRATION_H
