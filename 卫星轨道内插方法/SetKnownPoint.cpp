#include<iostream>
#include"DataStructure.h"
#include<vector>

using namespace std;
vector<Point> SetKnownPoint(vector<SourceData> data)
{
	int count = data.size();
	int step = count / 9;
	vector<Point> knownpoint;
	for (int i = 0; i < 9; i++)
	{
		Point temp;
		temp.x = data[i*step].x;
		temp.y = data[i*step].y;
		temp.z = data[i*step].z;
		temp.time = data[i*step].time;
		knownpoint.push_back(temp);
	}
	return knownpoint;
}