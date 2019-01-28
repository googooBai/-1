#include<vector>
#include"DataStructure.h"
using namespace std;




Point LagrangeInterpolating(vector<Point> data, Time t)
{
	Point result{ 0.0,0.0,0.0,0.0 };
	//查找内插点的位置
	int posfirst=0, posend=data.size()-1;
	int index=(posfirst+posend)/2;
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


	/*
	while (t > data[index].time)
	{
		++index;
		if (index == data.size())
			break;
	}
	*/
	//阶数
	int staffnum = 9;
	//内插起始点与终结点
	int start, end;
	if (index - staffnum / 2 < 0)
		start = 0;
	else
		start = index - staffnum / 2;
	if (index + staffnum/2> data.size() - 1)
		end = data.size() - 1;
	else
		end = index +staffnum/2;

	for (auto k=start;k!=end+1;k++)
	{
		double tempx = 1,tempy=1,tempz=1;
		for (auto i = start; i != end+1; i++)
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
		result.time = t;
	}
	return result;
}