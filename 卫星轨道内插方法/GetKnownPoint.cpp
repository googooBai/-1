#include"Function.h"

using namespace std;

namespace ClassInter {

	vector<Point> GetKnownPoint(vector<SourceData>& data,int interval)
	{
		for (auto i = data.begin(); i != data.end(); i++)
		{
			i->is_knownp = false;
		}
		Point temp;
		vector<Point> knownp;
		int index = 0;
		while(index <= data.size() - 1)
		{
			temp.time = data[index].time_d;
			temp.x = data[index].x;
			temp.y = data[index].y;
			temp.z = data[index].z;
			temp.xx = data[index].xx;
			temp.yy = data[index].yy;
			temp.zz = data[index].zz;
			knownp.push_back(temp);
			data[index].is_knownp = true;
			index += interval;
		}

		return knownp;
	}
}
