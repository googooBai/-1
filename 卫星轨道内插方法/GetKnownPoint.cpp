#include"DataStructure.h"
#include"Function.h"

using namespace std;

vector<Point> GetKnownPoint(vector<SourceData> data)
{
	Point temp;
	vector<Point> knownp;
	Time start= data.front().time;
	Time end = data.back().time;
	for (; start <= end; start += min * 15)
	{
		int pos_start, pos_end;
		pos_start = 0;
		pos_end = data.size() - 1;
		int index = (pos_start+pos_end)/2;
		while (pos_start < pos_end)
		{
			if (start >= data[index].time&&start < data[index + 1].time)
				break;
			else if (start < data[index].time)
			{
				pos_end = index;
				index = (pos_start + pos_end) / 2;
			}
			else
			{
				pos_start = index + 1;
				index = (pos_start + pos_end) / 2;
			}
		}
		temp.time = data[index].time;
		temp.x = data[index].x;
		temp.y = data[index].y;
		temp.z = data[index].z;
		knownp.push_back(temp);
	}
	return knownp;
}