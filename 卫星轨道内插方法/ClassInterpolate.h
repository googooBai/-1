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
		MatrixXd **T;
		MatrixXd **L;
		MatrixXd **A;
		MatrixXd **V;
		MatrixXd **P[3];
		double *error[3];
		int *time_knownp;
	    int *time_knownp_error;
		bool *is_knownp_error;
	public:
		explicit ChebyshevInterpolate(const vector<Point>& data, int slope);
		virtual ~ChebyshevInterpolate() {}
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

	//class Spline3 :public Interpolate
	//{
	//public:
	//	Spline3() {}
	//	virtual ~Spline3() {}
	//	virtual Point interpolate(vector<Point> data, Time t, int staff);
	//};

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
		virtual void CmpError(
			vector<SourceData>& data,
			const vector<Point>& ob,
			ostream& os,
			ostream& os2,
			ostream& os3,
			ostream& os4
		){}
		//virtual Point interpolate(vector<Point> data, Time t, int staff);
	};

	//class BSpline
	//{
	//private:
	//	MatrixXd *V;	//���Ƶ� x y z
	//	vector<Point> Q;	//��ֵ��
	//	MatrixXd *Qmat;		//��ֵ�����
	//private:
	//	double Base(Time t, int i);
	//	//void SetU(vector<double>& U, int knotsize, int slope);
	//	//double BaseFunction(int i,int k,double u, const vector<double>& U);	//���������
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
		MatrixXd **P[3];
		double *error[3];
		int *time_knownp;
		int *time_knownp_error;
		bool *is_knownp_error;
	private:
	public:
		explicit Polynomial(const vector<Point>& data, int slope);
		virtual ~Polynomial() {}
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
}