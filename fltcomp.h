#ifndef FLTCOMP_H
#define FLTCOMP_H
#include <math.h>
namespace FltCmp{
	inline bool aEqualRel(float a, float b, int maxUlps=4)
	{
	    int32_t aInt = *(int32_t*)&a;
	    // Make aInt lexicographically ordered as a twos-complement int
	    if(aInt < 0) aInt = 0x80000000 - aInt;
	    // Make bInt lexicographically ordered as a twos-complement int
	    int32_t bInt = *(int32_t*)&b;
	    if(bInt < 0) bInt = 0x80000000 - bInt;
	    int32_t intDiff = fabs(aInt - bInt);
	   
	    return ( fabs(aInt - bInt) <= maxUlps ) ? true : false;
	    //return ( abs(*(int32_t*)&a - *(int32_t*)&b) <= max Ulps ) ? true : false;
	}
	
	inline bool aEqualRel(double a, double b, int maxUlps=4)
	{
	    int64_t aInt = *(int64_t*)&a;
	    // Make aInt lexicographically ordered as a twos-complement int
	    if(aInt < 0) aInt = 0x8000000000000000 - aInt;
	    // Make bInt lexicographically ordered as a twos-complement int
	    int64_t bInt = *(int64_t*)&b;
	    if (bInt < 0) bInt = 0x8000000000000000 - bInt;
	    
	    return ( fabs(aInt - bInt) <= maxUlps ) ? true : false;
	    //return ( abs(*(int64_t*)&a - *(int64_t*)&b) <= max Ulps ) ? true : false;    
	}
	
	inline bool aGreaterRel(float a, float b, int maxUlps=4)
	{
		return (a >= b) ? true : aEqualRel(a, b, maxUlps);
	}
	
	inline bool aGreaterRel(double a, double b, int maxUlps=4)
	{
		return (a >= b) ? true : aEqualRel(a, b, maxUlps);
	}
	
	inline bool aLessRel(float a, float b, int maxUlps=4)
	{
		return (a <= b) ? true : aEqualRel(a, b, maxUlps);
	}
	
	inline bool aLessRel(double a, double b, int maxUlps=4)
	{
		return (a <= b) ? true : aEqualRel(a, b, maxUlps);
	}
	
	inline bool dGreaterRel(float a, float b, int maxUlps=4)
	{
		return (a > b && !aEqualRel(a, b, maxUlps)) ? true : false;
	}
	
	inline bool dGreaterRel(double a, double b, int maxUlps=4)
	{
		return (a > b && !aEqualRel(a, b, maxUlps)) ? true : false;
	}
	
	inline bool dLessRel(float a, float b, int maxUlps=4)
	{
		return (a < b && !aEqualRel(a, b, maxUlps)) ? true : false;
	}
	
	inline bool dLessRel(double a, double b, int maxUlps=4)
	{
		return (a < b && !aEqualRel(a, b, maxUlps)) ? true : false;
	}
	
	inline bool aEqualAbs(float a, float b, float absEps)
	{
		return ( fabs(a-b) <= absEps ) ? true : false;
	}
	
	inline bool aEqualAbs(double a, double b, double absEps)
	{
		return ( fabs(a-b) <= absEps ) ? true : false;
	}
	
	inline bool aGreaterAbs(float a, float b, float absEps)
	{
		return (a >= b-absEps) ? true : false;
	}
	
	inline bool aGreaterAbs(double a, double b, double absEps)
	{
		return (a >= b-absEps) ? true : false;
	}
	
	inline bool aLessAbs(float a, float b, float absEps)
	{
		return (a <= b+absEps) ? true : false;
	}
	
	inline bool aLessAbs(double a, double b, double absEps)
	{
		return (a <= b+absEps) ? true : false;
	}
	
	inline bool dGreaterAbs(float a, float b, float absEps)
	{
		return (a-absEps > b) ? true : false;
	}
	
	inline bool dGreaterAbs(double a, double b, double absEps)
	{
		return (a-absEps > b) ? true : false;
	}
	
	inline bool dLessAbs(float a, float b, float absEps)
	{
		return (a+absEps < b) ? true : false;
	}
	
	inline bool dLessAbs(double a, double b, double absEps)
	{
		return (a+absEps < b) ? true : false;
	}
}

#endif //FLTCOMP_H

