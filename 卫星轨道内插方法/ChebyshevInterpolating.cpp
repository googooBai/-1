#include<vector>
#include<iostream>
#include"DataStructure.h"
#include"Function.h"
#include"Eigen\Dense"

using namespace Eigen;
using namespace std;

Point ChebyshevInterpolating(vector<Point> data,Time t,int staff)
{
	//查找内插点位置
	int posfirst = 0, posend = data.size() - 1;
	int index=FindIndex(posfirst, posend, t, data);
	
	
	//内插起始点与终结点
	if (index > (data.size() - staff - 1))
		index = data.size() - staff - 1;
	//矩阵规模（（staff+1）*（staff*1））
	MatrixXd MatT(staff+1,staff+1);
	
	for (auto i = index; i <=index+staff; i++)
	{
		double tao = 2 * (data[i].time - data[index].time) / (data[index+staff].time-data[index].time) - 1;
		for (auto j = 0; j <=staff; j++)
		{
			if (j == 0)
				MatT(i-index,j) = 1;
			else if (j == 1)
				MatT(i-index,j) = tao;
			else
				MatT(i - index, j) = 2 * tao*MatT(i - index, j-1) - MatT(i - index, j-2);
		}
	}
	MatrixXd MatL(staff+1,3);
	for (auto i = index; i <=index+staff; i++)
	{
		MatL(i- index,0) = data[i].x;
		MatL(i - index, 1) = data[i].y;
		MatL(i - index, 2) = data[i].z;
	}
	


	MatrixXd parameter=(MatT.transpose()*MatT).inverse()*MatT.transpose()*MatL;
	

	Point result;
	double tao = 2 * (t - data[index].time) / (data[index+staff].time - data[index].time) - 1;
	VectorXd T(staff+1);
	for (auto j = 0; j <=staff; j++)
	{
		if (j == 0)
			T(j) = 1;
		else if (j == 1)
			T(j) = tao;
		else
			T(j) = 2 * tao*T(j-1) - T(j-2);
	}
	

	MatrixXd coor = T.transpose()*parameter;

	result.x = coor(0);
	result.y = coor(1);
	result.z = coor(2);
	result.time = t;
	return result;
}

