#include"Function.h"

Point Spline(vector<Point> data, Time t, int staff)
{
	int posfirst = 0, posend = data.size() - 1;
	int index = FindIndex(posfirst, posend, t, data);

	//处理边界
	if (index > (data.size() - staff - 1))
		index = data.size() - staff - 1;
	//复制已知点
	vector<Point> knownp(staff + 1);
	copy(data.begin() + index, data.begin() + index + staff+1, knownp.begin());
	//计算步长
	vector<double> h(staff);
	for (int i = 1; i <= staff; i++)
	{
		h[i - 1] = knownp[i].time - knownp[i - 1].time;
	}
	//计算S''(x);矩阵运算B*m=l
	MatrixXd B(staff + 1, staff + 1);
	MatrixXd l(staff + 1,3);
	for (int i = 0; i <= staff; i++)
	{
		if (i == 0 || i == staff) //首项或末项
		{
			//设置B
			B(i, i) = 1;
			//设置l
			l(i,0) =l(i,1)=l(i,2)= 0; 
		}
		else
		{	
			//设置B
			B(i, i) = 2 * (h[i - 1] + h[i]);
			B(i, i - 1) = h[i - 1];
			B(i, i + 1) = h[i];
			//设置l
			l(i, 0) = (knownp[i + 1].x - knownp[i].x) / h[i] -
				(knownp[i].x - knownp[i - 1].x) / h[i - 1];	    //x
			l(i, 1) = (knownp[i + 1].y - knownp[i].y) / h[i] -
				(knownp[i].y - knownp[i - 1].y) / h[i - 1];		//y
			l(i, 2) = (knownp[i + 1].z - knownp[i].z) / h[i] -
				(knownp[i].z - knownp[i - 1].z) / h[i - 1];		//z
		}
	}
	l *= 6;
	MatrixXd m(staff + 1,3);
	m = B.inverse()*l;	//[m0,m1,m2,m3,...,mn]T
	
	//计算系数
	MatrixXd Xabcd(staff, 4); //[ai,bi,ci,di]
	MatrixXd Yabcd(staff, 4);
	MatrixXd Zabcd(staff, 4);
	for (int i = 0; i < staff; i++)
	{
		//X
		Xabcd(i, 0) = knownp[i].x;
		Xabcd(i, 1) = (knownp[i + 1].x - knownp[i].x) / h[i] 
			- h[i] * m(i,0) / 2 
			- h[i] * (m(i+1,0) - m(i,0))/ 6;
		Xabcd(i, 2) = m(i,0) / 2;
		Xabcd(i, 3) = (m(i+1,0) - m(i,0)) / (6 * h[i]);

		//Y
		Yabcd(i, 0) = knownp[i].y;
		Yabcd(i, 1) = (knownp[i + 1].y - knownp[i].y) / h[i]
			- h[i] * m(i, 1) / 2
			- h[i] * (m(i + 1, 1) - m(i, 1)) / 6;
		Yabcd(i, 2) = m(i, 1) / 2;
		Yabcd(i, 3) = (m(i + 1, 1) - m(i, 1)) / (6 * h[i]);

		//Z
		Zabcd(i, 0) = knownp[i].z;
		Zabcd(i, 1) = (knownp[i + 1].z - knownp[i].z) / h[i]
			- h[i] * m(i, 2) / 2
			- h[i] * (m(i + 1, 2) - m(i, 2)) / 6;
		Zabcd(i, 2) = m(i, 2) / 2;
		Zabcd(i, 3) = (m(i + 1, 2) - m(i, 2)) / (6 * h[i]);
	}
	//计算坐标
	index = FindIndex(0, knownp.size() - 1, t, knownp);
	Point result;
	result.time = t;
	result.x 
		= Xabcd(index, 0)
		+ Xabcd(index, 1)*(t - knownp[index].time)
		+ Xabcd(index, 2)*(t - knownp[index].time)*(t - knownp[index].time)
		+ Xabcd(index, 3)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time);
	result.y 
		= Yabcd(index, 0)
		+ Yabcd(index, 1)*(t - knownp[index].time)
		+ Yabcd(index, 2)*(t - knownp[index].time)*(t - knownp[index].time)
		+ Yabcd(index, 3)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time);
	result.z 
		= Zabcd(index, 0)
		+ Zabcd(index, 1)*(t - knownp[index].time)
		+ Zabcd(index, 2)*(t - knownp[index].time)*(t - knownp[index].time)
		+ Zabcd(index, 3)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time);

	return result;
}