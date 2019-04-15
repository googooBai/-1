#include"ClassInterpolate.h"
#include"Function.h"

using namespace ClassInter;

int main()
{
	string filename{ "D:\\练习\\毕业设计\\Data\\几何轨道\\graceA-kinOrb-2008-01-01.txt" };
	vector<SourceData> data = GetData(filename);
	vector<Point> oa; vector<Point> ob;
	filename = "D:\\练习\\毕业设计\\Data\\动力学轨道\\GRACE-2008-01-01-02.DAT";
	GetOrbitPosition(oa, ob, filename);
	
	
	int type,slope,interval;
	while (1)
	{	
		cout << "0. Lar || 1. Che ||2. Poly \n";
		cin >> type;
		cin >> slope;
		cin >> interval;
		vector<Point> knownp = GetKnownPoint(data,interval);	
		ofstream fout,fout2;
	
		fout.open("output.txt");
		if (!fout.is_open())
			exit(EXIT_FAILURE);
		
		fout2.open("different.txt");
		if (!fout2.is_open())
			exit(EXIT_FAILURE);

		if (type == 0)
		{
			LarrangeInterpolate L(knownp, slope);
			CmpError(L, slope, data, knownp, oa, fout, fout2);
		}
		else if (type == 1)
		{
			ChebyshevInterpolate C(knownp, slope);
			CmpError(C, slope, data, knownp, oa, fout, fout2);

		}
		else if (type == 2)
		{
			Polynomial P(knownp, slope);
			CmpError(P, slope, data, knownp, oa, fout, fout2);
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

