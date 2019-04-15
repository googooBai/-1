#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	Polynomial::Polynomial(const vector<Point>& data, int slope) :knownp(data),Slope(slope)
	{
		//Caculate V=TA-L
		T = new MatrixXd*[data.size() - Slope - 1];
		A = new MatrixXd*[data.size() - Slope - 1];
		L = new MatrixXd*[data.size() - Slope - 1];
		V = new MatrixXd*[data.size() - Slope - 1];
		error[0] = new double[data.size() - Slope - 1];
		error[1] = new double[data.size() - Slope - 1];
		error[2] = new double[data.size() - Slope - 1];
		for (int i = 0; i <= knownp.size() - Slope - 1; i++)
		{
			T[i] = new MatrixXd(Slope + 1, Slope + 1);
			A[i] = new MatrixXd(Slope + 1, 3);
			L[i] = new MatrixXd(Slope + 1, 3);
			V[i] = new MatrixXd(Slope + 1, 3);

			for (int j = 0; j <= Slope; j++)
			{
				(*T[i])(j, 0) = 1;
				for (int k = 1; k <= Slope; k++)
					(*T[i])(j, k) = (*T[i])(j, k - 1) *
					(2 * (knownp[i + j].time - knownp[i].time) / (knownp[i + Slope].time - knownp[i].time) - 1);
				(*L[i])(j, 0) = this->knownp[i + j].x;
				(*L[i])(j, 1) = this->knownp[i + j].y;
				(*L[i])(j, 2) = this->knownp[i + j].z;
			}
			(*A[i]) = ((*T[i]).transpose()*(*T[i])).inverse()*(*T[i]).transpose()*(*L[i]);
			(*V[i]) = (*T[i])*(*A[i]) - (*L[i]);

			MatrixXd e = (*V[i]).transpose()*(*V[i]);
			error[0][i] = sqrt(e(0, 0) / Slope);
			error[1][i] = sqrt(e(1, 1) / Slope);
			error[2][i] = sqrt(e(2, 2) / Slope);
		}
	}

	Point Polynomial::interpolate(Time t)
	{
		int posfirst = 0, posend = knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, knownp);
		//¥¶¿Ì±ﬂΩÁ
		if (index < Slope / 2)
			index = 0;
		else if (index >(knownp.size() - (Slope - Slope / 2) - 1))
			index = knownp.size() - Slope - 1;
		else
			index = index - Slope / 2;
		Point result{ 0.0,0.0,0.0,0.0 };
		result.time = t;
		t = 2 * (t - knownp[index].time) / (knownp[index + Slope].time - knownp[index].time) - 1;
		for (int i = Slope; i >= 0; i--)
		{
			result.x *= t;
			result.x += (*A[index])(i, 0);
			result.y *= t;
			result.y += (*A[index])(i, 1);
			result.z *= t;
			result.z += (*A[index])(i, 2);
		}
		return result;
	}
}