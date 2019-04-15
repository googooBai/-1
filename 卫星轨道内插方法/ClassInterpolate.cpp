#include"ClassInterpolate.h"
#include"Function.h"

//namespace ClassInter {
//
//	Point Spline3::interpolate(vector<Point> data, Time t, int staff)
//	{
//		int posfirst = 0, posend = data.size() - 1;
//		int index = FindIndex(posfirst, posend, t, data);
//
//		Point result;
//		result.time = t;
//
//		//处理边界
//		if (index < staff / 2)
//			index = 0;
//		else if (index >(data.size() - (staff - staff / 2) - 1))
//			index = data.size() - staff - 1;
//		else
//			index = index - staff / 2;
//
//		//复制已知点
//		vector<Point> knownp(staff + 1);
//		copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());
//
//		//计算步长
//		vector<double> h(staff);
//		for (int i = 1; i <= staff; i++)
//		{
//			h[i - 1] = knownp[i].time - knownp[i - 1].time;
//		}
//
//		//计算S''(x);矩阵运算B*m=l
//		MatrixXd B=Eigen::MatrixXd::Zero(staff -1, staff -1);
//		MatrixXd l(staff -1, 3);
//		MatrixXd m(staff + 1, 3);
//		m(0, 0) = m(0, 1) = m(0, 2) =m(staff,0)= m(staff ,1)= m(staff,2)= 0;
//		for (int i = 0; i < staff - 1; i++)
//		{
//			if (i == 0)
//			{
//				B(i, i) = 2;
//				B(i, i + 1) = h[i + 1] / (h[i] + h[i + 1]);
//			}
//			else if (i == staff - 2)
//			{
//				B(i, i) = 2;
//				B(i, i - 1) = h[i] / (h[i] + h[i + 1]);
//			}
//			else
//			{
//				B(i, i) = 2;
//				B(i, i + 1) = h[i + 1] / (h[i] + h[i + 1]);
//				B(i, i - 1) = h[i] / (h[i] + h[i + 1]);
//			}
//			l(i, 0) = 6 * ((knownp[i + 2].x - knownp[i + 1].x) / h[i + 1] - (knownp[i + 1].x - knownp[i].x) / h[i]) / (h[i + 1] + h[i]);
//			l(i, 1) = 6 * ((knownp[i + 2].y - knownp[i + 1].y) / h[i + 1] - (knownp[i + 1].y - knownp[i].y) / h[i]) / (h[i + 1] + h[i]);
//			l(i, 2) = 6 * ((knownp[i + 2].z - knownp[i + 1].z) / h[i + 1] - (knownp[i + 1].z - knownp[i].z) / h[i]) / (h[i + 1] + h[i]);
//		}
//		m.middleRows(1,staff-1) = B.inverse()*l;	//[m0,m1,m2,m3,...,mn]'
//
//		//确定S(x)
//		index = FindIndex(0, knownp.size() - 1, t, knownp);	
//		double a1, a2, a3, a4;
//		//x
//		a1 = (m(index + 1, 0) - m(index, 0)) / (6 * h[index]);
//		a2 = 0.5*m(index, 0);
//		a3 = (knownp[index + 1].x - knownp[index].x) / h[index] - h[index] * m(index, 0) / 2 - h[index] * (m(index + 1, 0) - m(index, 0)) / 6;
//
//		a4 = knownp[index].x;
//		result.x =
//			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
//			a2*(t - knownp[index].time)*(t - knownp[index].time) +
//			a3*(t - knownp[index].time) +
//			a4;
//		//y
//		a1 = (m(index + 1, 1) - m(index, 1)) / (6 * h[index]);
//		a2 = 0.5*m(index, 1);
//		a3 = (knownp[index + 1].y - knownp[index].y) / h[index] - h[index] * m(index, 1) / 2 - h[index] * (m(index + 1, 1) - m(index, 1)) / 6;
//
//		a4 = knownp[index].y;
//		result.y =
//			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
//			a2*(t - knownp[index].time)*(t - knownp[index].time) +
//			a3*(t - knownp[index].time) +
//			a4;
//		//z
//		a1 = (m(index + 1, 2) - m(index, 2)) / (6 * h[index]);
//		a2 = 0.5*m(index, 2);
//		a3 = (knownp[index + 1].z - knownp[index].z) / h[index] - h[index] * m(index, 2) / 2 - h[index] * (m(index + 1, 2) - m(index, 2)) / 6;
//
//		a4 = knownp[index].z;
//		result.z =
//			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
//			a2*(t - knownp[index].time)*(t - knownp[index].time) +
//			a3*(t - knownp[index].time) +
//			a4;
//
//		return result;
//	}
//
//	Point Spline3v2::interpolate(vector<Point> data, Time t, int staff)
//	{
//		int posfirst = 0, posend = data.size() - 1;
//		int index = FindIndex(posfirst, posend, t, data);
//
//		Point result;
//		result.time = t;
//
//		//处理边界
//		if (index < staff / 2)
//			index = 0;
//		else if (index >(data.size() - (staff - staff / 2) - 1))
//			index = data.size() - staff - 1;
//		else
//			index = index - staff / 2;
//
//		//复制已知点
//		vector<Point> knownp(staff + 1);
//		copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());
//
//		//计算步长
//		vector<double> h(staff);
//		for (int i = 1; i <= staff; i++)
//		{
//			h[i - 1] = knownp[i].time - knownp[i - 1].time;
//		}
//
//		//BM = d
//		MatrixXd B = Eigen::MatrixXd::Zero(staff + 1, staff + 1);
//		MatrixXd M(staff + 1, 3);
//		MatrixXd d(staff + 1, 3);
//
//		for (int i = 0; i < staff + 1; i++)
//		{
//			if (i == 0||i==staff)
//			{
//				B(i, i) = 2;
//				d(i, 0) = d(i, 1) = d(i, 2) = 0;
//			}
//			else
//			{
//				B(i, i - 1) = h[i - 1] / (h[i - 1] + h[i]);
//				B(i, i) = 2;
//				B(i, i+1) = h[i] / (h[i - 1] + h[i]);
//				
//				d(i, 0) =
//					6 * (
//					(knownp[i + 1].x - knownp[i].x) / h[i] - (knownp[i].x - knownp[i - 1].x) / h[i - 1]
//						)
//					/ (h[i - 1] + h[i]);
//				d(i, 1) =
//					6 * (
//					(knownp[i + 1].y - knownp[i].y) / h[i] - (knownp[i].y - knownp[i - 1].y) / h[i - 1]
//						)
//					/ (h[i - 1] + h[i]);
//				d(i, 2) =
//					6 * (
//					(knownp[i + 1].z - knownp[i].z) / h[i] - (knownp[i].z - knownp[i - 1].z) / h[i - 1]
//						)
//					/ (h[i - 1] + h[i]);
//			}
//		}
//
//		M = B.inverse()*d;
//		index = FindIndex(0, knownp.size() - 1, t, knownp);
//		Point res;
//		res.time = t;
//		res.x
//			= M(index, 0)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
//			+ M(index + 1, 0)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
//			+ (knownp[index].x - M(index, 0)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
//			+ (knownp[index + 1].x - M(index + 1, 0)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
//		res.y
//			= M(index, 1)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
//			+ M(index + 1, 1)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
//			+ (knownp[index].y - M(index, 1)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
//			+ (knownp[index + 1].y - M(index + 1, 1)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
//		res.z
//			= M(index, 2)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
//			+ M(index + 1, 2)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
//			+ (knownp[index].z - M(index, 2)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
//			+ (knownp[index + 1].z - M(index + 1, 2)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
//
//		return res;
//	}
//	
//	BSpline::BSpline(const vector<Point>& data):Q(data) {
//		Qmat = new MatrixXd(Q.size()+2, 3);
//		(*Qmat)(0, 0) = (*Qmat)(0, 1) = (*Qmat)(0, 2) = (*Qmat)(Q.size() + 1, 0) = (*Qmat)(Q.size() + 1, 1) = (*Qmat)(Q.size() + 1, 2) = 0;
//		for (int i = 0; i < Q.size(); i++) {
//			(*Qmat)(i+1, 0) = Q[i].x;
//			(*Qmat)(i+1, 1) = Q[i].y;
//			(*Qmat)(i+1, 2) = Q[i].z;
//		}
//		
//		V = new MatrixXd(Q.size()+2, 3);
//		MatrixXd A(Q.size() + 2, Q.size() + 2);
//		for (int i = 0; i < Q.size() + 2; i++)
//		{
//			if (i == 0)
//			{
//				A(i, i) = 1;
//				A(i, i + 1) = -1;
//			}
//			else if (i == Q.size() + 1)
//			{
//				A(i, i) = 1;
//				A(i, i - 1) = -1;
//			}
//			else
//			{
//				A(i, i) = 4;
//				A(i, i - 1) = A(i, i + 1) = 1;
//			}
//		}
//		*V = A.inverse()*(*Qmat);
//	}
//	
//	Point BSpline::interpolate(Time t)
//	{
//		int index = FindIndex(0, Q.size(), t, Q);
//		Point res{ 0,0,0,0 };
//		if (index < Q.size() - 1)
//		{
//			res.time = t;
//			t = t - (int)t / 1.0;
//			res.x =
//				(*V)(index,0)*Base(t, 0) +
//				(*V)(index+1, 0)*Base(t, 1) +
//				(*V)(index+2, 0)*Base(t, 2) +
//				(*V)(index+3, 0)*Base(t, 3);
//			res.y =
//				(*V)(index, 1)*Base(t, 0) +
//				(*V)(index + 1, 1)*Base(t, 1) +
//				(*V)(index + 2, 1)*Base(t, 2) +
//				(*V)(index + 3, 1)*Base(t, 3);
//			res.z =
//				(*V)(index, 2)*Base(t, 0) +
//				(*V)(index + 1, 2)*Base(t, 1) +
//				(*V)(index + 2, 2)*Base(t, 2) +
//				(*V)(index + 3, 2)*Base(t, 3);
//			return res;
//		}
//		return res;	//异常处理
//	}
//
//	double BSpline::Base(Time t, int i)
//	{
//		switch (i)
//		{
//		case 0:
//			return (1 - t)*(1 - t)*(1 - t) / 6;
//		case 1:
//			return (3 * t*t*t - 6 * t*t + 4) / 6;
//		case 2:
//			return (-3 * t*t*t + 3 * t*t + 3*t + 1) / 6;
//		case 3:
//			return t*t*t / 6;
//		default:
//			break;
//		}
//		return 0.0;
//	}
//		
//}

