#include<vector>
#include"DataStructure.h"
#include"Function.h"
using namespace std;




Point LagrangeInterpolating(vector<Point> data, Time t, int staff)
{
	Point result{ 0.0,0.0,0.0,0.0 };
	//查找内插点的位置
	int posfirst=0, posend=data.size()-1;
	int index = FindIndex(posfirst, posend, t, data);
	//
	if (index > (data.size() - staff - 1))
		index = data.size() - staff - 1;

	for (auto k=index;k<=index+staff;k++)
	{
		double tempx = 1,tempy=1,tempz=1;
		for (auto i = index; i <=index+staff; i++)
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