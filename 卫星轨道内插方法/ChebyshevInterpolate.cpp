#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {
	ChebyshevInterpolate::ChebyshevInterpolate(const vector<Point>& data, int slope):Slope(slope),knownp(data)
	{
		T = new MatrixXd*[data.size() - Slope];
		L = new MatrixXd*[data.size() - Slope];
		A = new MatrixXd*[data.size() - Slope];
		V = new MatrixXd*[data.size() - Slope];
		error[0] = new double[data.size() - Slope];
		error[1] = new double[data.size() - Slope];
		error[2] = new double[data.size() - Slope ];
		for (int i = 0; i <= knownp.size() - Slope-1; i++)
		{
			T[i] = new MatrixXd(Slope + 1, Slope + 1);
			A[i] = new MatrixXd(Slope + 1, 3);
			L[i] = new MatrixXd(Slope + 1, 3);
			V[i] = new MatrixXd(Slope + 1, 3);

			//Caculate T
			for (int j = i; j <= i + Slope; j++)
			{
				double tao =
					2 * (knownp[j].time - knownp[i].time) /
					(knownp[i + Slope].time - knownp[i].time)
					- 1;
				for (int k = 0; k <= Slope; k++)
				{
					if (k == 0)
						(*T[i])(j - i, k) = 1;
					else if (k == 1)
						(*T[i])(j - i, k) = tao;
					else
						(*T[i])(j - i, k) = 2 * tao*(*T[i])(j - i, k - 1) - (*T[i])(j - i, k - 2);
				}
				//Set L
				(*L[i])(j - i, 0) = knownp[j].x;
				(*L[i])(j - i, 1) = knownp[j].y;
				(*L[i])(j - i, 2) = knownp[j].z;
			}

			(*A[i]) = ((*T[i]).transpose()*(*T[i])).inverse()*(*T[i]).transpose()*(*L[i]);
			(*V[i]) = (*T[i])*(*A[i]) - (*L[i]);

			
			MatrixXd e=(*V[i]).transpose()*(*V[i]);
			error[0][i] = sqrt(e(0, 0)/Slope);
			error[1][i] = sqrt(e(1, 1)/Slope);
			error[2][i] = sqrt(e(2, 2)/Slope);			
		}
	}

	Point ChebyshevInterpolate::interpolate(Time t)
	{
		Point res;
		res.time = t;

		int posfirst = 0, posend = this->knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, this->knownp);

		//内插起始点与终结点
		if (index < Slope / 2)
			index = 0;
		else if (index >(this->knownp.size() - (Slope - Slope / 2) - 1))
			index = this->knownp.size() - Slope - 1;
		else
			index = index - Slope / 2;

		double tao =
			2 * (t - this->knownp[index].time) /
			(this->knownp[index + Slope].time - this->knownp[index].time)
			- 1;

		VectorXd MatT(Slope + 1);
		for (auto j = 0; j <= Slope; j++)
		{
			if (j == 0)
				MatT(j) = 1;
			else if (j == 1)
				MatT(j) = tao;
			else
				MatT(j) = 2 * tao*MatT(j - 1) - MatT(j - 2);
		}

		MatrixXd coor = MatT.transpose()*(*A[index]);
		res.x = coor(0);
		res.y = coor(1);
		res.z = coor(2);
		
		return res;
	}

}