#include"Function.h"

using namespace std;

vector<Point> GetKnownPoint(vector<SourceData> data)
{
	Point temp;
	vector<Point> knownp;
	int index = 0;
	while(index <= data.size() - 1)
	{
		temp.time = data[index].time;
		temp.x = data[index].x;
		temp.y = data[index].y;
		temp.z = data[index].z;
		knownp.push_back(temp);
		index += 10;
	}
	if (data.size() != 1&&index!=data.size()-1)
	{
		temp.time = data[data.size() - 1].time;
		temp.x = data[data.size() - 1].x;
		temp.y = data[data.size() - 1].y;
		temp.z = data[data.size() - 1].z;
		knownp.push_back(temp);
	}
	return knownp;
}