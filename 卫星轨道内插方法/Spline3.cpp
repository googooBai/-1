#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {
	
	Spline3v2::Spline3v2(const vector<Point>& data, int slope):knownp(data),Slope(slope)
	{
		B = new MatrixXd*[knownp.size() - Slope];
		M = new MatrixXd*[knownp.size() - Slope];
		d = new MatrixXd*[knownp.size() - Slope];
		h = new vector<double>*[knownp.size() - Slope];
		for (int i = 0; i <= knownp.size() - Slope-1; i++)
		{
			B[i] = new MatrixXd(Slope + 1, Slope + 1);
			*B[i] = MatrixXd::Zero(Slope + 1, Slope + 1);
			M[i] = new MatrixXd(Slope + 1, 3);
			d[i] = new MatrixXd(Slope + 1, 3);
			h[i]=new vector<double>(Slope);

			for (int j = 0; j < Slope; j++)
			{
				(*h[i])[j] = knownp[i + j + 1].time - knownp[i + j].time;
			}

			for (int j = 0; j <= Slope; j++)
			{
				if (j == 0 || j == Slope)
				{
					(*B[i])(j, j) = 2;
					(*d[i])(j, 0) = (*d[i])(j, 1) = (*d[i])(j, 2) = 0;
				}
				else
				{
					(*B[i])(j, j - 1) = (*h[i])[j - 1] / ((*h[i])[j - 1] + (*h[i])[j]);
					(*B[i])(j, j) = 2;
					(*B[i])(j, j + 1) = (*h[i])[j] / ((*h[i])[j - 1] + (*h[i])[j]);

					(*d[i])(j, 0) =
						6 * (
						(knownp[i+j+1].x - knownp[i+j].x) / (*h[i])[j] - (knownp[i+j].x - knownp[i+j-1].x) / (*h[i])[j - 1]
							)
						/ ((*h[i])[j - 1] + (*h[i])[j]);
					(*d[i])(j, 1) =
						6 * (
						(knownp[i+j+1].y - knownp[i+j].y) / (*h[i])[j] - (knownp[i+j].y - knownp[i+j-1].y) / (*h[i])[j - 1]
							)
						/ ((*h[i])[j - 1] + (*h[i])[j]);
					(*d[i])(j, 2) =
						6 * (
						(knownp[i+j+1].z - knownp[i+j].z) / (*h[i])[j] - (knownp[i+j].z - knownp[i+j-1].z) / (*h[i])[j - 1]
							)
						/ ((*h[i])[j - 1] + (*h[i])[j]);
				}
			}

			(*M[i]) = (*B[i]).inverse()*(*d[i]);
		}
	}

	Point Spline3v2::interpolate(Time t)
	{
		Point res;
		res.time = t;
		int index = FindIndex(0, knownp.size() - 1, t, knownp);
		int index2 = index;
		if (index < Slope / 2)
			index = 0;
		else if (index >(knownp.size() - (Slope - Slope / 2) - 1))
			index = knownp.size() - Slope - 1;
		else
			index = index - Slope / 2;

		res.x
			= (*M[index])(index2 - index, 0)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 0)*(t - knownp[index2].time)*(t - knownp[index2].time)*(t - knownp[index2].time) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].x - (*M[index])(index2 - index, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(knownp[index2 + 1].time - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].x - (*M[index])(index2 - index + 1, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - knownp[index2].time) / (*h[index])[index2 - index];
		res.y
			= (*M[index])(index2 - index, 1)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 1)*(t - knownp[index2].time)*(t - knownp[index2].time)*(t - knownp[index2].time) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].y - (*M[index])(index2 - index, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(knownp[index2 + 1].time - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].y - (*M[index])(index2 - index + 1, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - knownp[index2].time) / (*h[index])[index2 - index];
		res.z
			= (*M[index])(index2 - index, 2)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t)*(knownp[index2 + 1].time - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 2)*(t - knownp[index2].time)*(t - knownp[index2].time)*(t - knownp[index2].time) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].z - (*M[index])(index2 - index, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(knownp[index2 + 1].time - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].z - (*M[index])(index2 - index + 1, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - knownp[index2].time) / (*h[index])[index2 - index];

		return res;
	}

	//Point Spline3v2::interpolate(vector<Point> data, Time t, int staff)
	//{
	//	int posfirst = 0, posend = data.size() - 1;
	//	int index = FindIndex(posfirst, posend, t, data);

	//	Point result;
	//	result.time = t;

	//	//处理边界
	//	if (index < staff / 2)
	//		index = 0;
	//	else if (index >(data.size() - (staff - staff / 2) - 1))
	//		index = data.size() - staff - 1;
	//	else
	//		index = index - staff / 2;

	//	//复制已知点
	//	vector<Point> knownp(staff + 1);
	//	copy(data.begin() + index, data.begin() + index + staff + 1, knownp.begin());

	//	//计算步长
	//	vector<double> h(staff);
	//	for (int i = 1; i <= staff; i++)
	//	{
	//		h[i - 1] = knownp[i].time - knownp[i - 1].time;
	//	}

	//	//BM = d
	//	MatrixXd B = Eigen::MatrixXd::Zero(staff + 1, staff + 1);
	//	MatrixXd M(staff + 1, 3);
	//	MatrixXd d(staff + 1, 3);

	//	for (int i = 0; i < staff + 1; i++)
	//	{
	//		if (i == 0 || i == staff)
	//		{
	//			B(i, i) = 2;
	//			d(i, 0) = d(i, 1) = d(i, 2) = 0;
	//		}
	//		else
	//		{
	//			B(i, i - 1) = h[i - 1] / (h[i - 1] + h[i]);
	//			B(i, i) = 2;
	//			B(i, i + 1) = h[i] / (h[i - 1] + h[i]);

	//			d(i, 0) =
	//				6 * (
	//				(knownp[i + 1].x - knownp[i].x) / h[i] - (knownp[i].x - knownp[i - 1].x) / h[i - 1]
	//					)
	//				/ (h[i - 1] + h[i]);
	//			d(i, 1) =
	//				6 * (
	//				(knownp[i + 1].y - knownp[i].y) / h[i] - (knownp[i].y - knownp[i - 1].y) / h[i - 1]
	//					)
	//				/ (h[i - 1] + h[i]);
	//			d(i, 2) =
	//				6 * (
	//				(knownp[i + 1].z - knownp[i].z) / h[i] - (knownp[i].z - knownp[i - 1].z) / h[i - 1]
	//					)
	//				/ (h[i - 1] + h[i]);
	//		}
	//	}

	//	M = B.inverse()*d;
	//	index = FindIndex(0, knownp.size() - 1, t, knownp);
	//	Point res;
	//	res.time = t;
	//	res.x
	//		= M(index, 0)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
	//		+ M(index + 1, 0)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
	//		+ (knownp[index].x - M(index, 0)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
	//		+ (knownp[index + 1].x - M(index + 1, 0)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
	//	res.y
	//		= M(index, 1)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
	//		+ M(index + 1, 1)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
	//		+ (knownp[index].y - M(index, 1)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
	//		+ (knownp[index + 1].y - M(index + 1, 1)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];
	//	res.z
	//		= M(index, 2)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t)*(knownp[index + 1].time - t) / 6 / h[index]
	//		+ M(index + 1, 2)*(t - knownp[index].time)*(t - knownp[index].time)*(t - knownp[index].time) / 6 / h[index]
	//		+ (knownp[index].z - M(index, 2)*h[index] * h[index] / 6)*(knownp[index + 1].time - t) / h[index]
	//		+ (knownp[index + 1].z - M(index + 1, 2)*h[index] * h[index] / 6)*(t - knownp[index].time) / h[index];

	//	return res;
	//}
}