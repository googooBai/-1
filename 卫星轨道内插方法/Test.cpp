#include"ClassInterpolate.h"
#include"Function.h"

using namespace ClassInter;

int main()
{
	string fileg[6]{
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-01.txt",
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-02.txt",
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-03.txt", 
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-04.txt", 
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-05.txt", 
		"D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-06.txt" };
	string filed[6]{
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-01-02.DAT",
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-02-02.DAT",
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-03-02.DAT",
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-04-02.DAT",
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-05-02.DAT",
		"D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-06-02.DAT"
	};
	string outkp_20s[6]{
		"kp_20s_1.txt",
		"kp_20s_2.txt",
		"kp_20s_3.txt",
		"kp_20s_4.txt",
		"kp_20s_5.txt",
		"kp_20s_6.txt"
	};

	string outip_20s[6]{
		"ip_20s_1.txt",
		"ip_20s_2.txt",
		"ip_20s_3.txt",
		"ip_20s_4.txt",
		"ip_20s_5.txt",
		"ip_20s_6.txt"
	};

	vector<Point> oa; vector<Point> ob;
	ofstream fout, fout2, fout3,fout4;
	int type,slope,interval;
	
	//for (int j = 0; j < 6; j++)
	//{
	//	vector<SourceData> data = GetData(fileg[j]);
	//	GetOrbitPosition(oa, ob, filed[j]);
	//	//fout4.open("knownp.txt");
	//	//if (!fout4.is_open())
	//	//	exit(EXIT_FAILURE);
	//	vector<Point> knownp = GetKnownPoint(data, 2);

	//	fout.open(outkp_20s[j]);
	//	fout2.open(outip_20s[j]);
	//	for (int i = 0; i < data.size(); i++)
	//	{
	//		int ob_index;
	//		{
	//			double t = data[i].time_d;
	//			int fp = 0, ep = oa.size() - 1;
	//			ob_index = (fp + ep) >> 1;
	//			while (fp < ep)
	//			{
	//				if (abs(t - oa[ob_index].time) <= 1e-6)
	//					break;
	//				else if (t - oa[ob_index].time > 1e-6)
	//				{
	//					fp = ob_index;
	//					ob_index = (fp + ep) >> 1;
	//				}
	//				else if (t - oa[ob_index].time < -1e-6)
	//				{
	//					ep = ob_index;
	//					ob_index = (fp + ep) >> 1;
	//				}
	//			}
	//		}
	//		data[i].de_x = oa[ob_index].x - data[i].x;
	//		data[i].de_y = oa[ob_index].y - data[i].y;
	//		data[i].de_z = oa[ob_index].z - data[i].z;

	//		if (data[i].is_knownp)
	//		{
	//			fout.setf(ios_base::fixed);
	//			fout.precision(15);
	//			fout << setw(22) << data[i].time_d;
	//			fout.precision(3);
	//			fout << setw(15) << data[i].x
	//				 << setw(15) << data[i].y
	//				 << setw(15) << data[i].z
	//				 << setw(15) << data[i].de_x
	//				 << setw(15) << data[i].de_y
	//				 << setw(15) << data[i].de_z
	//				 << endl;
	//		}
	//		else
	//		{
	//			fout2.setf(ios_base::fixed);
	//			fout2.precision(15);
	//			fout2 << setw(22) << data[i].time_d;
	//			fout2.precision(3);
	//			fout2 << setw(15) << data[i].x
	//				  << setw(15) << data[i].y
	//				  << setw(15) << data[i].z
	//				  << setw(15) << data[i].de_x
	//				  << setw(15) << data[i].de_y
	//				  << setw(15) << data[i].de_z
	//				  << endl;
	//		}
	//	}
	//	fout.close();
	//	fout2.close();
	//}

	//////////ChebyshevInterpolate C(knownp,26);
	//////////Polynomial P(knownp, 9);
	//Spline3v2 S(knownp, 10);
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
	////
	//{

	//	interval =1 ;
	//	type = 2;
	//	/*for (interval = 3; interval >= 1; interval--)
	//	{*/
	//		for (slope = 13; slope <= 13; slope++)
	//		{
	//			vector<Point> knownp = GetKnownPoint(data, 3);

	//			fout.open("g_error.txt");
	//			if (!fout.is_open())
	//				exit(EXIT_FAILURE);

	//			fout2.open("d_error.txt");
	//			if (!fout2.is_open())
	//				exit(EXIT_FAILURE);

	//			fout3.open("processed.txt");
	//			if (!fout3.is_open())
	//				exit(EXIT_FAILURE);

	//			fout4.open("knownp.txt");
	//			if (!fout4.is_open())
	//				exit(EXIT_FAILURE);

	//			if (type == 0)
	//			{
	//				cout << "Funciton:Larrange   Slope: " << slope << "    interval:   " << interval<<"min"<<endl;
	//				LarrangeInterpolate L(knownp, slope);
	//				L.CmpError(data, oa, fout, fout2, fout3, fout4);
	//			}
	//			else if (type == 1)
	//			{
	//				//fout4 << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//				cout << "Funciton:Chebyshev   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//				ChebyshevInterpolate C(knownp, slope);
	//				C.CmpError(data, oa, fout, fout2, fout3, fout4);

	//			}
	//			else if (type == 2)
	//			{
	//				//fout4 << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//				cout << "Funciton:Polynomial   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//				Polynomial P(knownp, slope);
	//				P.CmpError(data, oa, fout, fout2, fout3, fout4);
	//			}
	//			else if (type == 3)
	//			{
	//				cout << "Funciton:Spline   Slope: " << slope << "    interval:   " << interval << "min" << endl;
	//				Spline3v2 S(knownp, slope);
	//				S.CmpError(data, oa, fout, fout2, fout3, fout4);
	//			}
	//			else
	//				cout << "Wrong Input!\n";

	//			fout.close();
	//			fout2.close();
	//			fout3.close();
	//			fout4.close();
	//			cout << "complete\n";
	//		}
	//	//}
	//}

	{

		while (1)
		{
			cout << "file number:\n";
			int filenum =0;
			cin >> filenum;
			vector<SourceData> data = GetData(fileg[filenum]);
			GetOrbitPosition(oa, ob, filed[filenum]);
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
			fout3.close();
			fout4.close();
			cout << "complete\n";
		}	
	}

	system("pause");
	return 0;
}

