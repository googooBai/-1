#include<iostream>
#include<fstream>
#include"DataStructure.h"
#include"Function.h"
#include<vector>

using namespace std;

int main()
{
	string filename{ "D:\\��ϰ\\��ҵ���\\Data\\���ι��\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	double start, end;
	//start = data.front().time;
	//end = data.back().time;
	//const double min = 5 / (60 * 24);
	vector<Point> interpolating;
	for (auto i = data.begin(); i !=data.end(); i++)
	{
		interpolating.push_back(LagrangeInterpolating(data, i->time));
	}
	ofstream fout;
	fout.open("output.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);
	for (auto i = interpolating.begin(); i != interpolating.end(); i++)
	{
		fout << i->time << " " << i->x << " " << i->y << " " << i->z << endl;
	}

	system("pause");
	return 0;
}