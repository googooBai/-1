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
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		) = 0;		
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
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		);
	};

	class ChebyshevInterpolate:public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
		MatrixXd **A;
		Time* start;
		Time* end;
		double *error[3];
		bool *is_knownp_error;
	public:
		explicit ChebyshevInterpolate(const vector<Point>& data, int slope);
		virtual ~ChebyshevInterpolate() {}
		virtual Point interpolate(Time t);
		Point interpolate(Time t, ostream& os);
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		);
	};

	class Spline3v2 :public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
		MatrixXd **B;
		MatrixXd **M;
		MatrixXd **d;
		vector<double> **h;
	public:
		explicit Spline3v2(const vector<Point>& data,int slope);
		virtual ~Spline3v2() {}
		virtual Point interpolate(Time t);
		Point interpolate(Time t, ostream& os);
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		);
		//virtual Point interpolate(vector<Point> data, Time t, int staff);
	};


	class Polynomial :public Interpolate
	{
	private:
		int Slope;
		vector<Point> knownp;
		MatrixXd** A;
		double *error[3];
		bool *is_knownp_error;
	private:
	public:
		explicit Polynomial(const vector<Point>& data, int slope);
		virtual ~Polynomial() {}
		virtual Point interpolate(Time t);
		Point interpolate(Time t, ostream& os);
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		);
	};	
}
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