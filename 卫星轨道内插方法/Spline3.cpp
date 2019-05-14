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
				(*h[i])[j] = Normalization(knownp, i, Slope, knownp[i + j + 1].time) - Normalization(knownp, i, Slope, knownp[i + j].time);
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

		t = Normalization(knownp, index, Slope, t);
		res.x
			= (*M[index])(index2 - index, 0)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 0)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].x - (*M[index])(index2 - index, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].x - (*M[index])(index2 - index + 1, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		res.y
			= (*M[index])(index2 - index, 1)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 1)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].y - (*M[index])(index2 - index, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].y - (*M[index])(index2 - index + 1, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		res.z
			= (*M[index])(index2 - index, 2)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 2)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].z - (*M[index])(index2 - index, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].z - (*M[index])(index2 - index + 1, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		/*res.x
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
			+ (knownp[index2 + 1].z - (*M[index])(index2 - index + 1, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - knownp[index2].time) / (*h[index])[index2 - index];*/

		return res;
	}

	Point Spline3v2::interpolate(Time t,ostream& os)
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

		t = Normalization(knownp, index, Slope, t);
		res.x
			= (*M[index])(index2 - index, 0)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 0)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].x - (*M[index])(index2 - index, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].x - (*M[index])(index2 - index + 1, 0)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		res.y
			= (*M[index])(index2 - index, 1)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 1)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].y - (*M[index])(index2 - index, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].y - (*M[index])(index2 - index + 1, 1)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		res.z
			= (*M[index])(index2 - index, 2)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / 6 / (*h[index])[index2 - index]
			+ (*M[index])(index2 - index + 1, 2)*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time))*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / 6 / (*h[index])[index2 - index]
			+ (knownp[index2].z - (*M[index])(index2 - index, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(Normalization(knownp, index, Slope, knownp[index2 + 1].time) - t) / (*h[index])[index2 - index]
			+ (knownp[index2 + 1].z - (*M[index])(index2 - index + 1, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - Normalization(knownp, index, Slope, knownp[index2].time)) / (*h[index])[index2 - index];
		/*res.x
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
		+ (knownp[index2 + 1].z - (*M[index])(index2 - index + 1, 2)*(*h[index])[index2 - index] * (*h[index])[index2 - index] / 6)*(t - knownp[index2].time) / (*h[index])[index2 - index];*/
		//test
		os << "inter point: " << endl;
		os
			<< setw(22) << "Time"
			<< setw(15) << "x"
			<< setw(15) << "y"
			<< setw(15) << "z"
			<< endl;
		os.setf(ios_base::fixed);
		os.precision(15);
		os << setw(22) << res.time;
		os.precision(3);
		os
			<< setw(15) << res.x
			<< setw(15) << res.y
			<< setw(15) << res.z
			<< endl;

		os << "known point: " << endl;
		os
			<< setw(22) << "Time"
			<< setw(15) << "x"
			<< setw(15) << "y"
			<< setw(15) << "z"
			<< setw(25) << "xx"
			<< setw(25) << "yy"
			<< setw(25) << "zz"
			<< endl;
		for (int j = index; j <= index + Slope; j++)
		{
			os.setf(ios_base::fixed);
			os.precision(15);
			os << setw(22) << knownp[j].time;
			os.precision(3);
			os
				<< setw(15) << knownp[j].x
				<< setw(15) << knownp[j].y
				<< setw(15) << knownp[j].z
				;
			os.precision(15);
			os
				<< setw(25) << knownp[j].xx
				<< setw(25) << knownp[j].yy
				<< setw(25) << knownp[j].zz
				<< endl;
		}
		return res;
	}

	void Spline3v2::CmpError(
		vector<SourceData>& data,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2,
		ostream& os3,
		ostream& os4
	)
	{
		double error_x = 0.0, error_y = 0.0, error_z = 0.0;
		double RMSE_X = 0.0, RMSE_Y = 0.0, RMSE_Z = 0.0;
		int not_knownp = 0;
		os
			<< setw(22) << "Time"
			<< setw(15) << "error x"
			<< setw(15) << "error y"
			<< setw(15) << "error z"
			<< endl;
		os2
			<< setw(22) << "Time"
			<< setw(15) << "error x"
			<< setw(15) << "error y"
			<< setw(15) << "error z"
			<< endl;

		os3
			<< setw(22) << "Time"
			<< setw(15) << "g_error x"
			<< setw(15) << "g_error y"
			<< setw(15) << "g_error z"
			<< setw(15) << "d_error x"
			<< setw(15) << "d_error y"
			<< setw(15) << "d_error z"
			<< endl;

		for (int i = 0; i < data.size(); i++)
		{
			data[i].is_effect_by_knownp = false;
			data[i].is_g = false;
			data[i].is_d = false;
			data[i].is_big_error = false;
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= this->Slope / 2 && index <= (knownp.size() - (this->Slope - this->Slope / 2) - 1))
			{
				if (!data[i].is_knownp)
				{
					Point temp = this->interpolate(data[i].time_d);
					data[i].ix = temp.x;
					data[i].iy = temp.y;
					data[i].iz = temp.z;
				}
				else
				{
					data[i].ix = data[i].x;
					data[i].iy = data[i].y;
					data[i].iz = data[i].z;
				}

				int ob_index;
				{
					double t = data[i].time_d;
					int fp = 0, ep = ob.size() - 1;
					ob_index = (fp + ep) >> 1;
					while (fp < ep)
					{
						if (abs(t - ob[ob_index].time) <= 1e-6)
							break;
						else if (t - ob[ob_index].time > 1e-6)
						{
							fp = ob_index;
							ob_index = (fp + ep) >> 1;
						}
						else if (t - ob[ob_index].time < -1e-6)
						{
							ep = ob_index;
							ob_index = (fp + ep) >> 1;
						}
					}
				}

				data[i].de_x = ob[ob_index].x - data[i].x;
				data[i].de_y = ob[ob_index].y - data[i].y;
				data[i].de_z = ob[ob_index].z - data[i].z;

				//x
				data[i].error_x = data[i].ix - data[i].x;
				if (!data[i].is_knownp)
				{
					RMSE_X += data[i].error_x*data[i].error_x;
					not_knownp++;
				}
				//y
				data[i].error_y = data[i].iy - data[i].y;
				if (!data[i].is_knownp)
					RMSE_Y += data[i].error_y*data[i].error_y;
				//z
				data[i].error_z = data[i].iz - data[i].z;
				if (!data[i].is_knownp)
					RMSE_Z += data[i].error_z*data[i].error_z;

				os.setf(ios_base::fixed);
				os.precision(15);
				os << setw(22) << data[i].time_d;
				os.precision(3);
				os
					<< setw(15) << data[i].error_x
					<< setw(15) << data[i].error_y
					<< setw(15) << data[i].error_z
					<< endl;

				os2.setf(ios_base::fixed);
				os2.precision(15);
				os2 << setw(22) << data[i].time_d;
				os2.precision(3);
				os2
					<< setw(15) << data[i].de_x
					<< setw(15) << data[i].de_y
					<< setw(15) << data[i].de_z
					<< endl;
			}
		}
		RMSE_X = sqrt(RMSE_X / not_knownp);
		RMSE_Y = sqrt(RMSE_Y / not_knownp);
		RMSE_Z = sqrt(RMSE_Z / not_knownp);

		//粗差处理
		int total{ 0 };
		int mistake{ 0 };
		int hit{ 0 };
		int knownp_index = 0;
		int big_error = 0;
		int hit_big_error = 0;
		int cx = 0, cy = 0, cz = 0;

		for (int i = 0; i < data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= this->Slope / 2 && index <= (knownp.size() - (this->Slope - this->Slope / 2) - 1))
			{
				if (abs(data[i].de_x) > 0.06 ||
					abs(data[i].de_y) > 0.06 ||
					abs(data[i].de_z) > 0.06)
				{
					data[i].is_d = true;
					total++;
				}
				else
					data[i].is_d = false;

				if (data[i].is_knownp)
				{
					data[i].is_g = false;
					os4.setf(ios_base::fixed);
					os4.precision(15);
					os4 << setw(22) << data[i].time_d;
					os4.precision(3);
					os4 << setw(15) << data[i].de_x
						<< setw(15) << data[i].de_y
						<< setw(15) << data[i].de_z
						<< setw(15) << data[i].is_d
						<< endl;
					if (data[i].is_d)
					{
						abs(data[i].de_x)>0.06 ? cx++ : cx;
						abs(data[i].de_y)>0.06 ? cy++ : cy;
						abs(data[i].de_z)>0.06 ? cz++ : cz;

						for (int j = 1; j <= 5; j++)
						{
							data[i - j].is_effect_by_knownp = data[i - j].is_knownp ? false : true;
							data[i + j].is_effect_by_knownp = data[i + j].is_knownp ? false : true;
						}
					}
				}
				else
				{
					if (abs(data[i].error_x) > 2 * RMSE_X ||
						abs(data[i].error_y) > 2 * RMSE_Y ||
						abs(data[i].error_z) > 2 * RMSE_Z)
						data[i].is_g = true;
					else
						data[i].is_g = false;


					if (data[i].is_g&&data[i].is_d)
						hit++;
					else if (data[i].is_g && !data[i].is_d)
						mistake++;
				}

				if (abs(data[i].de_x) > 0.1 ||
					abs(data[i].de_y) > 0.1 ||
					abs(data[i].de_z) > 0.1)
				{
					big_error++;
					data[i].is_big_error = true;
					if (data[i].is_g)
						hit_big_error++;
				}

				if (!data[i].is_g)
				{
					os3.setf(ios_base::fixed);
					os3.precision(15);
					os3 << setw(22) << data[i].time_d;
					os3.precision(3);
					os3
						<< setw(15) << data[i].error_x
						<< setw(15) << data[i].error_y
						<< setw(15) << data[i].error_z
						<< setw(15) << data[i].de_x
						<< setw(15) << data[i].de_y
						<< setw(15) << data[i].de_z
						<< endl;
				}
			}
		}

		int n_effectbyknownp = 0;
		double ef_RMSE_X = 0.0, ef_RMSE_Y = 0.0, ef_RMSE_Z = 0.0;
		int ef_error = 0, ef_hit = 0, ef_mistake = 0;
		int ef_big = 0, ef_hit_big = 0;
		for (int i = 0; i < data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= this->Slope / 2 && index <= (knownp.size() - (this->Slope - this->Slope / 2) - 1))
			{
				if (data[i].is_effect_by_knownp)
				{
					n_effectbyknownp++;
					ef_RMSE_X += data[i].error_x*data[i].error_x;
					ef_RMSE_Y += data[i].error_y*data[i].error_y;
					ef_RMSE_Z += data[i].error_z*data[i].error_z;
					if (data[i].is_d)
					{
						ef_error++;
					}
					data[i].is_g&&data[i].is_d ? ef_hit++ : ef_hit;
					data[i].is_g && !data[i].is_d ? ef_mistake : ef_mistake;
					if (data[i].is_big_error)
					{
						ef_big++;
						data[i].is_g ? ef_hit_big++ : ef_hit_big;
					}
				}
			}
		}
		ef_RMSE_X = sqrt(ef_RMSE_X / n_effectbyknownp);
		ef_RMSE_Y = sqrt(ef_RMSE_Y / n_effectbyknownp);
		ef_RMSE_Z = sqrt(ef_RMSE_Z / n_effectbyknownp);

		std::cout << "Amount of error: " << total << endl;
		std::cout << "Amount of big error(more than 1cm): " << big_error << "	  eliminate: " << hit_big_error << endl;
		std::cout.setf(ios_base::fixed);
		std::cout.precision(1);
		std::cout << "Rate of sucessful: %" << (double)hit / total * 100 << endl;
		std::cout << "Rate of mistake: %" << (double)mistake / data.size() * 100 << endl;
		std::cout.precision(3);
		std::cout
			<< setw(8) << RMSE_X
			<< setw(8) << RMSE_Y
			<< setw(8) << RMSE_Z
			<< endl;
		std::cout << "Amount of effected by known point which have gross: " << n_effectbyknownp << endl;
		std::cout << "Amount of effected by konwn is error: " << ef_error << endl;
		std::cout << "Amount of effected by konwon is big error： " << ef_big << "    eliminate: " << ef_hit_big << endl;
		std::cout << "cx: " << cx << "	cy: " << cy << "    cz: " << cz << endl;
		std::cout.precision(1);
		std::cout << "Rate of sucessful: %" << (double)ef_hit / ef_error * 100 << endl;
		std::cout << "Rate of mistake: %" << (double)ef_mistake / n_effectbyknownp * 100 << endl;
		std::cout.precision(3);
		std::cout
			<< setw(8) << ef_RMSE_X
			<< setw(8) << ef_RMSE_Y
			<< setw(8) << ef_RMSE_Z
			<< endl;
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