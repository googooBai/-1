#include"ClassInterpolate.h"
#include"Function.h"

Point LarrangeInterpolate::interpolate(vector<Point> data, Time t, int staff)
{
	//�����ڲ���λ��
	int posfirst = 0, posend = data.size() - 1;
	int index = FindIndex(posfirst, posend, t, data);
	//����߽�
	if (index < staff / 2)
		index = 0;
	else if (index > (data.size() - (staff - staff /2) - 1))
		index = data.size() - staff-1 ;
	else
		index = index - staff /2;
	Point result{ 0.0,0.0,0.0,0.0 };
	result.time = t;
    //��һ��
	//t = Normalization(data, index, staff, t);
	//�ڲ����
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
	//�����ڲ��λ��
	int posfirst = 0, posend = data.size() - 1;
	int index = FindIndex(posfirst, posend, t, data);


	//�ڲ���ʼ�����ս��
	if (index < staff / 2)
		index = 0;
	else if (index >(data.size() - (staff - staff / 2) - 1))
		index = data.size() - staff - 1;
	else
		index = index - staff / 2;
	//�����ģ����staff+1��*��staff*1����
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
	//����߽�
	if (index < staff / 2)
		index = 0;
	else if (index >(data.size() - (staff - staff / 2) - 1))
		index = data.size() - staff - 1;
	else
		index = index - staff / 2;
	//������֪��
	vector<Point> knownp(staff + 1);
	copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());
	//t=Normalization(knownp, 0, staff, t);
	//���㲽��
	vector<double> h(staff);
	for (int i = 1; i <= staff; i++)
	{
		h[i - 1] = knownp[i].time - knownp[i - 1].time;
	}
	//����S''(x);��������B*m=l
	MatrixXd B=Eigen::MatrixXd::Zero(staff + 1, staff + 1);
	MatrixXd l(staff + 1, 3);
	for (int i = 0; i <= staff; i++)
	{
		if (i == 0 || i == staff) //�����ĩ��
		{
			//����B
			B(i, i) = 1;
			//����l
			l(i, 0) = l(i, 1) = l(i, 2) = 0;
		}
		else
		{
			//����B
			B(i, i) = 2 * (h[i - 1] + h[i]);
			B(i, i - 1) = h[i - 1];
			B(i, i + 1) = h[i];
			//����l
			l(i, 0) = (knownp[i + 1].x - knownp[i].x) / h[i] -
				(knownp[i].x - knownp[i - 1].x) / h[i - 1];	    //x
			l(i, 1) = (knownp[i + 1].y - knownp[i].y) / h[i] -
				(knownp[i].y - knownp[i - 1].y) / h[i - 1];		//y
			l(i, 2) = (knownp[i + 1].z - knownp[i].z) / h[i] -
				(knownp[i].z - knownp[i - 1].z) / h[i - 1];		//z
		}
	}
	l *= 6;
	MatrixXd m(staff + 1, 3);
	m = B.inverse()*l;	//[m0,m1,m2,m3,...,mn]'

						//����ϵ��
	MatrixXd Xabcd(staff, 4); //[ai,bi,ci,di]
	MatrixXd Yabcd(staff, 4);
	MatrixXd Zabcd(staff, 4);
	for (int i = 0; i < staff; i++)
	{
		//X
		Xabcd(i, 0) = knownp[i].x;
		Xabcd(i, 1) = (knownp[i + 1].x - knownp[i].x) / h[i]
			- h[i] * m(i, 0) / 2
			- h[i] * (m(i + 1, 0) - m(i, 0)) / 6;
		Xabcd(i, 2) = m(i, 0) / 2;
		Xabcd(i, 3) = (m(i + 1, 0) - m(i, 0)) / (6 * h[i]);

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
	//��������
	index = FindIndex(0, knownp.size() - 1, t, knownp);
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
