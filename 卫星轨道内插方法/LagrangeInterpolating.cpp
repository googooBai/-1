#include<vector>
#include"DataStructure.h"
using namespace std;




Point LagrangeInterpolating(vector<Point> known, Time t)
{
	Point result{ 0.0,0.0,0.0,0.0 };
	for (auto k=known.begin();k!=known.end();k++)
	{
		double tempx = 1,tempy=1,tempz=1;
		for (auto i = known.begin(); i != known.end(); i++)
		{
			if (k != i)
			{
				tempx *= (t - i->time) / (k->time - i->time);
				tempy *= (t - i->time) / (k->time - i->time);
				tempz *= (t - i->time) / (k->time - i->time);
			}			
		}
		result.x += tempx*k->x;
		result.y += tempy*k->y;
		result.z += tempz*k->z;
	}
	return result;
}