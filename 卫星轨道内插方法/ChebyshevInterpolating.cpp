#include<vector>
#include<iostream>
#include"DataStructure.h"
#include"Eigen\Dense"

using namespace Eigen;
using namespace std;

Point ChebyshevInterpolating(vector<SourceData> data,Time t)
{
	//�����ڲ��λ��
	int posfirst = 0, posend = data.size() - 1;
	int index = (posfirst + posend) / 2;
	while (posfirst < posend)
	{
		if (t >= data[index].time&&t < data[index + 1].time)
			break;
		if (t < data[index].time)
		{
			posend = index;
			index = (posfirst + posend) / 2;
		}
		if (t >= data[index + 1].time)
		{
			posfirst = index + 1;
			index = (posfirst + posend) / 2;
		}
	}
	int staffnum = 9;
	//�ڲ���ʼ�����ս��
	int start, end;
	if (index - staffnum / 2 < 0)
		start = 0;
	else
		start = index - staffnum / 2;
	if (index + staffnum / 2 > data.size() - 1)
		end = data.size() - 1;
	else
		end = index + staffnum / 2;
	//�����ģ����֪����*10��9�ף���
	MatrixXd MatT(end-start+1,10);
	
	for (auto i = start; i !=end+1; i++)
	{
		double tao = 2 * (data[i].time - data[start].time) / (data[end].time-data[start].time) - 1;
		for (auto j = 0; j !=10; j++)
		{
			if (j == 0)
				MatT(i-start,j) = 1;
			else if (j == 1)
				MatT(i-start,j) = tao;
			else
				MatT(i - start, j) = 2 * tao*MatT(i - start, j-1) - MatT(i - start, j-2);
		}
	}
	MatrixXd MatL(end-start+1,3);
	for (auto i = start; i !=end+1; i++)
	{
		MatL(i-start,0) = data[i].x;
		MatL(i - start, 1) = data[i].y;
		MatL(i - start, 2) = data[i].z;
	}
	//test
	cout << MatT;
	cout << endl << endl;
	cout << MatL << endl << endl;


	MatrixXd parameter=(MatT.transpose()*MatT).inverse()*MatT.transpose()*MatL;
	//test
	cout << parameter << endl;

	Point result;
	double tao = 2 * (t - data[start].time) / (data[end].time - data[start].time) - 1;
	VectorXd T(10);
	for (auto j = 0; j != 10; j++)
	{
		if (j == 0)
			T(j) = 1;
		else if (j == 1)
			T(j) = tao;
		else
			T(j) = 2 * tao*T(j-1) - T(j-2);
	}
	//test
	cout << T << endl;

	MatrixXd coor = T.transpose()*parameter;

	//test
	cout << coor << endl;
	result.x = coor(0);
	result.y = coor(1);
	result.z = coor(2);
	result.time = t;
	return result;
}

