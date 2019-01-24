#include<iostream>
#include"DataStructure.h"
#include"Function.h"
#include<vector>

using namespace std;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> interpolating;
	for (auto i = data.begin(); i != data.end(); i++)
	{
		interpolating.push_back(LagrangeInterpolating(data, i->time));		
	}
	system("pause");
	return 0;
}