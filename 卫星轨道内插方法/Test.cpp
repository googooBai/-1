#include"ClassInterpolate.h"
#include"Function.h"

using namespace ClassInter;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-06.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> oa; vector<Point> ob;
	filename = "D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-06-02.DAT";
	GetOrbitPosition(oa, ob, filename);
	
	
	int type,slope,interval;
	//		ofstream fout, fout2, fout3,fout4;
	//		fout4.open("knownp.txt");
	//		if (!fout4.is_open())
	//			exit(EXIT_FAILURE);
	//		type = 1; interval = 2;
	///*for (interval = 3; interval >= 1; interval--)
	//{*/
	//	for (slope = 7; slope <= 25; slope++)
	//	{
	//		vector<Point> knownp = GetKnownPoint(data, interval*6);

	//		fout.open("g_error.txt");
	//		if (!fout.is_open())
	//			exit(EXIT_FAILURE);

	//		fout2.open("d_error.txt");
	//		if (!fout2.is_open())
	//			exit(EXIT_FAILURE);

	//		fout3.open("processed.txt");
	//		if (!fout3.is_open())
	//			exit(EXIT_FAILURE);


	//		if (type == 0)
	//		{
	//			fout4 << "Funciton:Larrange   Slope: " << slope << "    interval:   " << interval<<"min"<<endl;
	//			LarrangeInterpolate L(knownp, slope);
	//			L.CmpError(data, oa, fout, fout2, fout3, fout4);
	//		}
	//		else if (type == 1)
	//		{
	//			fout4 << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//			cout << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//			ChebyshevInterpolate C(knownp, slope);
	//			C.CmpError(data, oa, fout, fout2, fout3, fout4);

	//		}
	//		else if (type == 2)
	//		{
	//			fout4 << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//			cout << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//			Polynomial P(knownp, slope);
	//			P.CmpError(data, oa, fout, fout2, fout3, fout4);
	//		}
	//		else if (type == 3)
	//		{
	//			Spline3v2 S(knownp, slope);
	//			CmpError(S, slope, data, knownp, oa, fout, fout2, fout3);
	//		}
	//		else
	//			cout << "Wrong Input!\n";

	//		fout.close();
	//		fout2.close();
	//		fout3.close();
	//		cout << "complete\n";
	//	}
	////}
	//		fout4.close();
	while (1)
	{	
		cout << "0. Lar || 1. Che ||2. Poly ||3. Spline \n";
		cin >> type;
		cin >> slope;
		cin >> interval;
		vector<Point> knownp = GetKnownPoint(data,interval);	
		ofstream fout,fout2,fout3,fout4;
	
		fout.open("g_error.txt");
		if (!fout.is_open())
			exit(EXIT_FAILURE);
		
		fout2.open("d_error.txt");
		if (!fout2.is_open())
			exit(EXIT_FAILURE);

		fout3.open("processed.txt");
		if (!fout3.is_open())
			exit(EXIT_FAILURE);

		fout4.open("knownp.txt");
		if (!fout3.is_open())
			exit(EXIT_FAILURE);

		if (type == 0)
		{
			LarrangeInterpolate L(knownp, slope);
			L.CmpError(data, oa, fout, fout2, fout3,fout4);
			//CmpError(L, slope, data, knownp, oa, fout, fout2, fout3);
		}
		else if (type == 1)
		{
			ChebyshevInterpolate C(knownp, slope);
			C.CmpError(data, oa, fout, fout2, fout3, fout4);

		}
		else if (type == 2)
		{
			Polynomial P(knownp, slope);
			P.CmpError(data,oa, fout, fout2, fout3, fout4);
		}
		else if (type == 3)
		{
			Spline3v2 S(knownp, slope);
			//CmpError(S, slope, data, knownp, oa, fout, fout2, fout3, fout4);
		}
		else
			cout << "Wrong Input!\n";
			
		fout.close();
		fout2.close();
		cout << "complete\n";
	}	

	system("pause");
	return 0;
}

