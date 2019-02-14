#include<iostream>
#include<fstream>
#include"DataStructure.h"
#include"Function.h"
#include<vector>

using namespace std;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> knownp = GetKnownPoint(data);
	double start, end;
	start = knownp.front().time;
	end = knownp.back().time;
	vector<Point> interpolating1,interpolating2;
	while (start <= end)
	{
		interpolating1.push_back(LagrangeInterpolating(knownp, start,11));
		interpolating2.push_back(ChebyshevInterpolating(knownp, start,11));
		start += min/6; //10s
	}
	ofstream fout;
	fout.open("output.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);
	
	fout.setf(ios::fixed);
	/*for (auto i = interpolating1.begin(); i != interpolating1.end(); i++)
	{
		fout.precision(15);
		fout << i->time;
		fout.precision(4);
		fout<< " " << i->x << " " << i->y << " " << i->z << endl;
	}*/
	for (auto i = 0; i < interpolating1.size(); i++)
	{
		fout.precision(15);
		fout << interpolating1[i].time;
		fout.precision(4);
		fout << " " << interpolating1[i].x-interpolating2[i].x << " "
			<< interpolating1[i].y - interpolating2[i].y << " "
			<< interpolating1[i].z - interpolating2[i].z
			<< endl;
	}
	fout.close();
	

	Point temp;
	fout.open("different.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);
	fout.setf(ios::fixed);
	fout.precision(15);
	fout.precision(6);
	for (auto i = data.begin(); i !=data.end(); i++)
	{
		temp=LagrangeInterpolating(knownp, i->time,5);
		fout << temp.time;

		fout << " " << temp.x - i->x << " "
			<< temp.y - i->y << " "
			<< temp.z - i->z
			<< endl;
	}
	fout.close();
	system("pause");
	return 0;
}

