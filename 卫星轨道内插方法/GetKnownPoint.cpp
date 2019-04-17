#include"Function.h"

using namespace std;

namespace ClassInter {

	vector<Point> GetKnownPoint(vector<SourceData> data,int interval)
	{
		Point temp;
		vector<Point> knownp;
		int index = 0;
		while(index <= data.size() - 1)
		{
			temp.time = data[index].time_d;
			temp.x = data[index].x;
			temp.y = data[index].y;
			temp.z = data[index].z;
			knownp.push_back(temp);
			index += interval;
		}

		return knownp;
	}
}
