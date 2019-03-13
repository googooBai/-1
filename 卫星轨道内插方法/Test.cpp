#include"ClassInterpolate.h"
#include"Function.h"



int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-02.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> knownp = GetKnownPoint(data);
	double start, end;
	start = knownp.front().time;
	end = knownp.back().time;
	LarrangeInterpolate L;
	ChebyshevInterpolate C;
	Spline3 S;
	/*vector<Point> x(101); 
	vector<SourceData> u(100001);
	for (int i = 0; i <= 5; i++)
	{
		x[i].time =-5+i;
		x[i].x = 1.0/(x[i].time*x[i].time+1);
	}
	for (int i = 0; i <=25; i++)
	{
		u[i].time = -5+0.2*i;
		u[i].x= 1.0 / (u[i].time*u[i].time + 1);
		Point nn =S.interpolate(x, u[i].time, 5);
		cout << nn.time << " " << nn.x << endl;
	}*/

	//while (start <= end)
	//{
	//	//interpolating2.push_back(C.interpolate(knownp, start,9));
	//	interpolating1.push_back(Spline(knownp, start,11));
	//	start += minute; //10s
	//}

	
	ofstream fout;
	
	fout.open("output.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);


	CmpError(S, 15, data, knownp, fout);
	
	fout.close();
	/*fout.setf(ios::fixed);
	for (auto i = interpolating2.begin(); i != interpolating2.end(); i++)
	{
		fout.precision(15);
		fout << i->time;
		fout.precision(4);
		fout<< " " << i->x << " " << i->y << " " << i->z << endl;
	}
	for (auto i = 0; i < interpolating1.size(); i++)
	{
		fout.precision(15);
		fout << interpolating1[i].time;
		fout.precision(4);
		fout << " " << interpolating1[i].x-interpolating2[i].x << " "
			<< interpolating1[i].y - interpolating2[i].y << " "
			<< interpolating1[i].z - interpolating2[i].z
			<< endl;
	}*/
	
	

	/*Point temp;
	fout.open("output.txt");
	if (!fout.is_open())
		exit(EXIT_FAILURE);
	fout.setf(ios::fixed);
	for (auto i = data.begin(); i !=data.end(); i++)
	{
		temp=Spline(knownp, i->time,11);
		fout.precision(15);
		fout << temp.time;

		fout.precision(6);
		fout << "	 " << temp.x-i->x  << "	 "
			<< temp.y -i->y << "	 "
			<< temp.z -i->z
			<< endl;
	}
	fout.close();*/
	system("pause");
	return 0;
}

