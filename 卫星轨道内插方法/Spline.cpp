#include"Function.h"

Point Spline(vector<Point> data, Time t, int staff)
{
	int posfirst = 0, posend = data.size() - 1;
	int index = FindIndex(posfirst, posend, t, data);

	//处理边界
	if (index < staff / 2)
		index = 0;
	else if (index >(data.size() - (staff - staff / 2) - 1))
		index = data.size() - staff - 1;
	else
		index = index - staff / 2;
	//复制已知点
	vector<Point> knownp(staff + 1);
	copy(data.begin() + index, data.begin() + index + staff+1, knownp.begin());
	//计算步长
	vector<double> h(staff);
	for (int i = 1; i <= staff; i++)
	{
		h[i - 1] = knownp[i].time - knownp[i - 1].time;
	}
	
	//B*M=D
	MatrixXd B= Eigen::MatrixXd::Zero(staff -1, staff -1);
	MatrixXd d(staff - 1,3);
	MatrixXd M(staff + 1, 3);
	M(0, 0) = M(0, 1) = M(0, 2) = M(staff, 0) = M(staff, 1) = M(staff, 2) = 0;
	for (int i = 0; i < staff-1 ; i++)
	{
		if (i == 0)
		{
			B(i, i) = 2;
			B(i, i + 1) = h[i+1] / (h[i] + h[i+1]);
		}
		else if (i == staff -2 )
		{
			B(i, i) = 2;
			B(i, i - 1) = h[i] / (h[i] + h[i+1]);
		}
		else
		{
			B(i, i) = 2;
			B(i, i + 1) = h[i + 1] / (h[i] + h[i + 1]);
			B(i, i - 1) = h[i] / (h[i] + h[i + 1]);
		}
	    d(i,0) =6*((knownp[i + 2].x - knownp[i+1].x) / h[i+1] - (knownp[i+1].x - knownp[i ].x) / h[i]) / (h[i+1] + h[i]);
		d(i, 1) = 6 * ((knownp[i + 2].y - knownp[i + 1].y) / h[i + 1] - (knownp[i + 1].y - knownp[i].y) / h[i]) / (h[i + 1] + h[i]);
		d(i, 2) = 6 * ((knownp[i + 2].z - knownp[i + 1].z) / h[i + 1] - (knownp[i + 1].z - knownp[i].z) / h[i]) / (h[i + 1] + h[i]);
	}

	M.middleRows(1,staff-1) = B.inverse()*d;
	index = FindIndex(0, knownp.size() - 1, t, knownp);
	//S(x)=a(x-xk)^3+b(x-xk)^2+c(x-xk)+d
	//系数
	double a1, a2, a3,a4 ;
	Point result;
	result.time = t;
	//x
	a1 = (M(index + 1, 0) - M(index, 0)) / (6 * h[index]);
	a2 = 0.5*M(index, 0);
	a3 = (knownp[index + 1].x - knownp[index].x) / h[index] - h[index] * M(index, 0) / 2 - h[index] * (M(index + 1, 0) - M(index, 0)) / 6;
																							
	a4 = knownp[index].x;
	result.x =
		a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
		a2*(t - knownp[index].time)*(t - knownp[index].time) +
		a3*(t - knownp[index].time) +
		a4;		
	//y
	a1 = (M(index + 1, 1) - M(index, 1)) / (6 * h[index]);
	a2 = 0.5*M(index, 1);
	a3 = (knownp[index + 1].y - knownp[index].y) / h[index] - h[index] * M(index, 1) / 2 - h[index] * (M(index + 1, 1) - M(index, 1)) / 6;
																							
	a4 = knownp[index].y;
	result.y =
		a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
		a2*(t - knownp[index].time)*(t - knownp[index].time) +
		a3*(t - knownp[index].time) +
		a4;
	//z
	a1 = (M(index + 1, 2) - M(index, 2)) / (6 * h[index]);
	a2 = 0.5*M(index, 2);
	a3 = (knownp[index + 1].z - knownp[index].z) / h[index] - h[index] * M(index, 2) / 2 -h[index] * (M(index + 1, 2) - M(index, 2)) / 6;
																							
	a4 = knownp[index].z;
	result.z =
		a1*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) +
		a2*(t - knownp[index].time)*(t - knownp[index].time) +
		a3*(t - knownp[index].time) +
		a4;
	return result;
}