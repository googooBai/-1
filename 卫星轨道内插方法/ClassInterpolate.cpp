#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	Point LarrangeInterpolate::interpolate(vector<Point> data, Time t, int staff)
	{
		//查找内插点的位置
		int posfirst = 0, posend = data.size() - 1;
		int index = FindIndex(posfirst, posend, t, data);
		//处理边界
		if (index < staff / 2)
			index = 0;
		else if (index > (data.size() - (staff - staff /2) - 1))
			index = data.size() - staff-1 ;
		else
			index = index - staff /2;
		Point result{ 0.0,0.0,0.0,0.0 };
		result.time = t;
		//内插过程
		for (auto k = index; k <= index + staff; k++)
		{
			double tempx = 1, tempy = 1, tempz = 1;
			for (auto i = index; i <= index + staff; i++)
			{
				if (k != i)
				{
					tempx *= (t - data[i].time) / (data[k].time - data[i].time);
					tempy *= (t - data[i].time) / (data[k].time - data[i].time);
					tempz *= (t - data[i].time) / (data[k].time - data[i].time);
				}
			}
			result.x += tempx*data[k].x;
			result.y += tempy*data[k].y;
			result.z += tempz*data[k].z;
		}
		return result;
	}

	Point ChebyshevInterpolate::interpolate(vector<Point> data, Time t, int staff)
	{
		//查找内插点位置
		int posfirst = 0, posend = data.size() - 1;
		int index = FindIndex(posfirst, posend, t, data);


		//内插起始点与终结点
		if (index < staff / 2)
			index = 0;
		else if (index >(data.size() - (staff - staff / 2) - 1))
			index = data.size() - staff - 1;
		else
			index = index - staff / 2;
		//矩阵规模（（staff+1）*（staff*1））
		MatrixXd MatT(staff + 1, staff + 1);

		for (auto i = index; i <= index + staff; i++)
		{
			double tao = 2 * (data[i].time - data[index].time) / (data[index + staff].time - data[index].time) - 1;
			for (auto j = 0; j <= staff; j++)
			{
				if (j == 0)
					MatT(i - index, j) = 1;
				else if (j == 1)
					MatT(i - index, j) = tao;
				else
					MatT(i - index, j) = 2 * tao*MatT(i - index, j - 1) - MatT(i - index, j - 2);
			}
		}
		MatrixXd MatL(staff + 1, 3);
		for (auto i = index; i <= index + staff; i++)
		{
			MatL(i - index, 0) = data[i].x;
			MatL(i - index, 1) = data[i].y;
			MatL(i - index, 2) = data[i].z;
		}

		MatrixXd parameter = (MatT.transpose()*MatT).inverse()*MatT.transpose()*MatL;

		Point result;
		double tao = 2 * (t - data[index].time) / (data[index + staff].time - data[index].time) - 1;
		VectorXd T(staff + 1);
		for (auto j = 0; j <= staff; j++)
		{
			if (j == 0)
				T(j) = 1;
			else if (j == 1)
				T(j) = tao;
			else
				T(j) = 2 * tao*T(j - 1) - T(j - 2);
		}

		MatrixXd coor = T.transpose()*parameter;

		result.x = coor(0);
		result.y = coor(1);
		result.z = coor(2);
		result.time = t;
		return result;
	}

	Point Spline3::interpolate(vector<Point> data, Time t, int staff)
	{
		int posfirst = 0, posend = data.size() - 1;
		int index = FindIndex(posfirst, posend, t, data);

		Point result;
		result.time = t;

		//处理边界
		if (index < staff / 2)
			index = 0;
		else if (index >(data.size() - (staff - staff / 2) - 1))
			index = data.size() - staff - 1;
		else
			index = index - staff / 2;

		//复制已知点
		vector<Point> knownp(staff + 1);
		copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());

		//计算步长
		vector<double> h(staff);
		for (int i = 1; i <= staff; i++)
		{
			h[i - 1] = knownp[i].time - knownp[i - 1].time;
		}

		//计算S''(x);矩阵运算B*m=l
		MatrixXd B=Eigen::MatrixXd::Zero(staff -1, staff -1);
		MatrixXd l(staff -1, 3);
		MatrixXd m(staff + 1, 3);
		m(0, 0) = m(0, 1) = m(0, 2) =m(staff,0)= m(staff ,1)= m(staff,2)= 0;
		for (int i = 0; i < staff - 1; i++)
		{
			if (i == 0)
			{
				B(i, i) = 2;
				B(i, i + 1) = h[i + 1] / (h[i] + h[i + 1]);
			}
			else if (i == staff - 2)
			{
				B(i, i) = 2;
				B(i, i - 1) = h[i] / (h[i] + h[i + 1]);
			}
			else
			{
				B(i, i) = 2;
				B(i, i + 1) = h[i + 1] / (h[i] + h[i + 1]);
				B(i, i - 1) = h[i] / (h[i] + h[i + 1]);
			}
			l(i, 0) = 6 * ((knownp[i + 2].x - knownp[i + 1].x) / h[i + 1] - (knownp[i + 1].x - knownp[i].x) / h[i]) / (h[i + 1] + h[i]);
			l(i, 1) = 6 * ((knownp[i + 2].y - knownp[i + 1].y) / h[i + 1] - (knownp[i + 1].y - knownp[i].y) / h[i]) / (h[i + 1] + h[i]);
			l(i, 2) = 6 * ((knownp[i + 2].z - knownp[i + 1].z) / h[i + 1] - (knownp[i + 1].z - knownp[i].z) / h[i]) / (h[i + 1] + h[i]);
		}
		m.middleRows(1,staff-1) = B.inverse()*l;	//[m0,m1,m2,m3,...,mn]'

		//确定S(x)
		index = FindIndex(0, knownp.size() - 1, t, knownp);	
		double a1, a2, a3, a4;
		//x
		a1 = (m(index + 1, 0) - m(index, 0)) / (6 * h[index]);
		a2 = 0.5*m(index, 0);
		a3 = (knownp[index + 1].x - knownp[index].x) / h[index] - h[index] * m(index, 0) / 2 - h[index] * (m(index + 1, 0) - m(index, 0)) / 6;

		a4 = knownp[index].x;
		result.x =
			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
			a2*(t - knownp[index].time)*(t - knownp[index].time) +
			a3*(t - knownp[index].time) +
			a4;
		//y
		a1 = (m(index + 1, 1) - m(index, 1)) / (6 * h[index]);
		a2 = 0.5*m(index, 1);
		a3 = (knownp[index + 1].y - knownp[index].y) / h[index] - h[index] * m(index, 1) / 2 - h[index] * (m(index + 1, 1) - m(index, 1)) / 6;

		a4 = knownp[index].y;
		result.y =
			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
			a2*(t - knownp[index].time)*(t - knownp[index].time) +
			a3*(t - knownp[index].time) +
			a4;
		//z
		a1 = (m(index + 1, 2) - m(index, 2)) / (6 * h[index]);
		a2 = 0.5*m(index, 2);
		a3 = (knownp[index + 1].z - knownp[index].z) / h[index] - h[index] * m(index, 2) / 2 - h[index] * (m(index + 1, 2) - m(index, 2)) / 6;

		a4 = knownp[index].z;
		result.z =
			a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
			a2*(t - knownp[index].time)*(t - knownp[index].time) +
			a3*(t - knownp[index].time) +
			a4;

		return result;
	}

	Point Spline3v2::interpolate(vector<Point> data, Time t, int staff)
	{
		int posfirst = 0, posend = data.size() - 1;
		int index = FindIndex(posfirst, posend, t, data);

		Point result;
		result.time = t;

		//处理边界
		if (index < staff / 2)
			index = 0;
		else if (index >(data.size() - (staff - staff / 2) - 1))
			index = data.size() - staff - 1;
		else
			index = index - staff / 2;

		//复制已知点
		vector<Point> knownp(staff + 1);
		copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());

		//计算步长
		vector<double> h(staff);
		for (int i = 1; i <= staff; i++)
		{
			h[i - 1] = knownp[i].time - knownp[i - 1].time;
		}

		//BM = d
		MatrixXd B = Eigen::MatrixXd::Zero(staff + 1, staff + 1);
		MatrixXd M(staff + 1, 3);
		MatrixXd d(staff + 1, 3);

		for (int i = 0; i < staff + 1; i++)
		{
			if (i == 0||i==staff)
			{
				B(i, i) = 2;
				d(i, 0) = d(i, 1) = d(i, 2) = 0;
			}
			else
			{
				B(i, i - 1) = h[i - 1] / (h[i - 1] + h[i]);
				B(i, i) = 2;
				B(i, i+1) = h[i] / (h[i - 1] + h[i]);
				
				d(i, 0) =
					6 * (
					(knownp[i + 1].x - knownp[i].x) / h[i] - (knownp[i].x - knownp[i - 1].x) / h[i - 1]
						)
					/ (h[i - 1] + h[i]);
				d(i, 1) =
					6 * (
					(knownp[i + 1].y - knownp[i].y) / h[i] - (knownp[i].y - knownp[i - 1].y) / h[i - 1]
						)
					/ (h[i - 1] + h[i]);
				d(i, 2) =
					6 * (
					(knownp[i + 1].z - knownp[i].z) / h[i] - (knownp[i].z - knownp[i - 1].z) / h[i - 1]
						)
					/ (h[i - 1] + h[i]);
			}
		}

		M = B.inverse()*d;
		index = FindIndex(0, knownp.size() - 1, t, knownp);
		Point res;
		res.time = t;
		res.x
			= M(index, 0)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
			+ M(index + 1, 0)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
			+ (knownp[index].x - M(index, 0)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
			+ (knownp[index + 1].x - M(index + 1, 0)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
		res.y
			= M(index, 1)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
			+ M(index + 1, 1)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
			+ (knownp[index].y - M(index, 1)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
			+ (knownp[index + 1].y - M(index + 1, 1)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
		res.z
			= M(index, 2)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
			+ M(index + 1, 2)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
			+ (knownp[index].z - M(index, 2)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
			+ (knownp[index + 1].z - M(index + 1, 2)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];

		return res;
	}

	//void BSpline::SetU(vector<double>& U,int knotsize,int slope)
	//{
	//	int m = knotsize + slope;	//m=n+1+p
	//	for (int i = 0; i <= slope; i++)
	//		U.push_back(0.0);
	//	for (int i = slope + 1; i <= m - slope - 1; i++)
	//		U.push_back(U[i - 1] + 1 / (m - 2 * slope));
	//	for (int i = m - slope; i <= m; i++)
	//		U.push_back(1.0);
	//}

	//double BSpline::BaseFunction(int i, int k, double u, const vector<double>& U)
	//{
	//	if (k == 0)
	//		if (u >= U[i]&&u < U[i + 1])
	//			return 1.0;
	//		else
	//			return 0.0;
	//	double length1 = U[i + k] - U[i];
	//	double length2 = U[i + k + 1] - U[i + 1];
	//	if (length1 == 0.0)
	//		length1 = 1.0;
	//	if (length2 == 0.0)
	//		length2 = 0.0;
	//	return (u - U[i])*BaseFunction(i, k - 1, u, U) / length1
	//		+ (U[i + k + 1] - u) * BaseFunction(i + 1, k - 1, u, U) / length2;
	//}

	//Point BSpline::interpolate(vector<Point> data,Time t,int slope)
	//{
	//	int index = FindIndex(0, data.size(), t, data);
	//	Point res{ t,0.0,0.0,0.0 };
	//	for (int i = index; i <= index + slope; i++)
	//	{
	//		res.x += data[i].x*BaseFunction(i, slope, t, data);
	//		res.y += data[i].y*BaseFunction(i, slope, t, data);
	//		res.z += data[i].z*BaseFunction(i, slope, t, data);
	//	}
	//	return res;
	//}
}

