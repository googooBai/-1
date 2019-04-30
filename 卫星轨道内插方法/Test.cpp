#include"ClassInterpolate.h"
#include"Function.h"

using namespace ClassInter;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-05.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> oa; vector<Point> ob;
	filename = "D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-05-02.DAT";
	GetOrbitPosition(oa, ob, filename);
	
	int type,slope,interval;
	ofstream fout, fout2, fout3,fout4;

	//fout4.open("knownp.txt");
	//if (!fout4.is_open())
	//	exit(EXIT_FAILURE);

	//vector<Point> knownp = GetKnownPoint(data, 18);

	////ChebyshevInterpolate C(knownp,26);
	////Polynomial P(knownp, 9);
	//Spline3v2 S(knownp, 18);
	//fout4 << "original point: " << endl;
	//fout4
	//	<< setw(22) << "Time"
	//	<< setw(15) << "x"
	//	<< setw(15) << "y"
	//	<< setw(15) << "z"
	//	<< endl;
	//fout4.setf(ios_base::fixed);
	//fout4.precision(15);
	//fout4 << setw(22) << data[2000].time_d;
	//fout4.precision(3);
	//fout4
	//	<< setw(15) << data[2000].x
	//	<< setw(15) << data[2000].y
	//	<< setw(15) << data[2000].z
	//	<< endl;
	//S.interpolate(data[2000].time_d,fout4);
	//fout4.close();
	
	/*{

		interval = 1;
		type = 3;
		for (interval = 3; interval >= 1; interval--)
		{
			for (slope = 18; slope <= 35; slope++)
			{
				vector<Point> knownp = GetKnownPoint(data, interval*6);

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
				if (!fout4.is_open())
					exit(EXIT_FAILURE);

				if (type == 0)
				{
					fout4 << "Funciton:Larrange   Slope: " << slope << "    interval:   " << interval<<"min"<<endl;
					LarrangeInterpolate L(knownp, slope);
					L.CmpError(data, oa, fout, fout2, fout3, fout4);
				}
				else if (type == 1)
				{
					fout4 << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
					cout << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
					ChebyshevInterpolate C(knownp, slope);
					C.CmpError(data, oa, fout, fout2, fout3, fout4);

				}
				else if (type == 2)
				{
					fout4 << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
					cout << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
					Polynomial P(knownp, slope);
					P.CmpError(data, oa, fout, fout2, fout3, fout4);
				}
				else if (type == 3)
				{
					cout << "Funciton:Spline   Slope: " << slope << "    interval:   " << interval << "min" << endl;
					Spline3v2 S(knownp, slope);
					S.CmpError(data, oa, fout, fout2, fout3, fout4);
				}
				else
					cout << "Wrong Input!\n";

				fout.close();
				fout2.close();
				fout3.close();
				fout4.close();
				cout << "complete\n";
			}
		}
	}*/

	{

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
				S.CmpError(data, oa, fout, fout2, fout3, fout4);
			}
			else
				cout << "Wrong Input!\n";
			
			fout.close();
			fout2.close();
			cout << "complete\n";
		}	
	}

	system("pause");
	return 0;
}

