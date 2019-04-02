#pragma once
#include"DataStructure.h"
#include"include.h"
namespace ClassInter {

	class Interpolate
	{
	public:
		Interpolate(){}
		virtual ~Interpolate(){}
		virtual  Point interpolate(vector<Point> data, Time t, int staff) = 0;
	};


	class LarrangeInterpolate :public Interpolate
	{
	public:
		LarrangeInterpolate(){}
		virtual ~LarrangeInterpolate(){}
		virtual Point interpolate(vector<Point> data, Time t, int staff);
	};

	class ChebyshevInterpolate :public Interpolate
	{
	public:
		ChebyshevInterpolate(){}
		virtual ~ChebyshevInterpolate(){}
		virtual Point interpolate(vector<Point> data, Time t, int staff);
	};

	class Spline3 :public Interpolate
	{
	public:
		Spline3() {}
		virtual ~Spline3() {}
		virtual Point interpolate(vector<Point> data, Time t, int staff);
	};

	class Spline3v2 :public Interpolate
	{
	public:
		Spline3v2() {}
		virtual ~Spline3v2() {}
		virtual Point interpolate(vector<Point> data, Time t, int staff);
	};

	//class BSpline :public Interpolate
	//{
	//private:
	//	vector<double> U;
	//private:
	//	void SetU(vector<double>& U, int knotsize, int slope);
	//	double BaseFunction(int i,int k,double u, const vector<double>& U);	//¼ÆËã»ùº¯Êý
	//public:
	//	BSpline(){}
	//	virtual ~BSpline() {}
	//	virtual Point interpolate(vector<Point> data, Time t, int slope);
	//};

}