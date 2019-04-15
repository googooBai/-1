#pragma once
#include"DataStructure.h"
#include"include.h"
namespace ClassInter {

	class Interpolate
	{
	public:
		Interpolate(){}
		virtual ~Interpolate(){}
		virtual  Point interpolate(Time t) = 0;
	};


	class LarrangeInterpolate :public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
	public:
		explicit LarrangeInterpolate(const vector<Point>&data,int slope)
			:Slope(slope),knownp(data){}
		virtual ~LarrangeInterpolate(){}
		virtual Point interpolate(Time t);
	};

	class ChebyshevInterpolate:public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
		MatrixXd **T;
		MatrixXd **L;
		MatrixXd **A;
		MatrixXd **V;
		double *error[3];
	public:
		explicit ChebyshevInterpolate(const vector<Point>& data, int slope);
		virtual ~ChebyshevInterpolate(){}
		virtual Point interpolate(Time t);
	};

	//class Spline3 :public Interpolate
	//{
	//public:
	//	Spline3() {}
	//	virtual ~Spline3() {}
	//	virtual Point interpolate(vector<Point> data, Time t, int staff);
	//};

	//class Spline3v2 :public Interpolate
	//{
	//public:
	//	Spline3v2() {}
	//	virtual ~Spline3v2() {}
	//	virtual Point interpolate(vector<Point> data, Time t, int staff);
	//};

	//class BSpline
	//{
	//private:
	//	MatrixXd *V;	//控制点 x y z
	//	vector<Point> Q;	//型值点
	//	MatrixXd *Qmat;		//型值点矩阵
	//private:
	//	double Base(Time t, int i);
	//	//void SetU(vector<double>& U, int knotsize, int slope);
	//	//double BaseFunction(int i,int k,double u, const vector<double>& U);	//计算基函数
	//public:
	//	explicit BSpline(const vector<Point>& data);
	//	virtual ~BSpline() {}
	//	virtual Point interpolate(Time t);
	//};

	class Polynomial :public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
		MatrixXd** T;
		MatrixXd** A;
		MatrixXd** L;
		MatrixXd** V;
		double *error[3];
	private:
	public:
		explicit Polynomial(const vector<Point>& data, int slope);
		virtual ~Polynomial(){}
		virtual Point interpolate(Time t);
	};		
}