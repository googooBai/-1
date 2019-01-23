#include<vector>
#include"DataStructure.h"
#include"/¡∑œ∞/æÿ’Û/Class_Matrix/Class_Matrix/Class_Matrix.h"
using namespace std;

Matrix ChebyshevInterpolating(vector<Point> known,Time dt,Time t);
Matrix ChebyshevInterpolating(vector<Point> known,Time dt,Time t)
{
	Matrix MatT(known.size(), known.size());
	//vector<double> *Ti = new vector<double>[known.size()];
	for (auto i = 0; i <known.size(); i++)
	{
		double tao = 2 * (known[i].time - known[0].time) / dt - 1;
		for (auto j = 0; j < known.size(); j++)
		{
			if (j == 0)
				MatT.data[i*known.size() + j] = 0;
			else if (j == 1)
				MatT.data[i*known.size() + j] = tao;
			else
				MatT.data[i*known.size() + j] = 2 * tao*MatT.data[i*known.size() + j - 1] - MatT.data[i*known.size() + j - 2];			
		}
	}
	Matrix MatL(known.size(), 3);
	for (auto i = 0; i < known.size(); i++)
	{
		MatL.data[i*known.size() + 0] = known[i].position.x;
		MatL.data[i*known.size() + 1] = known[i].position.y;
		MatL.data[i*known.size() + 2] = known[i].position.z;
	}
	return (MatT.transpose()*MatT).inverse()*(MatT.transpose()*MatL).inverse();
}

