#include<iostream>
#include"DataStructure.h"
#include"Function.h"
#include<vector>

using namespace std;

int main()
{
	string filename{ "D:\\��ϰ\\��ҵ���\\Data\\���ι��\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> interpolating;
	for (auto i = data.begin(); i != data.end(); i++)
	{
		interpolating.push_back(LagrangeInterpolating(data, i->time));		
	}
	system("pause");
	return 0;
}