#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {
	
	Point LarrangeInterpolate::interpolate(Time t)
	{
		//查找内插点的位置
		int posfirst = 0, posend = this->knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, this->knownp);
		//处理边界
		if (index < this->Slope / 2)
			index = 0;
		else if (index >(this->knownp.size() - (this->Slope - this->Slope / 2) - 1))
			index = this->knownp.size() - this->Slope - 1;
		else
			index = index - this->Slope / 2;
		Point result{ 0.0,0.0,0.0,0.0 };
		result.time = t;
		//内插过程
		for (auto k = index; k <= index + this->Slope; k++)
		{
			double tempx = 1, tempy = 1, tempz = 1;
			for (auto i = index; i <= index + this->Slope; i++)
			{
				if (k != i)
				{
					tempx *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
					tempy *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
					tempz *= (t - this->knownp[i].time) / (this->knownp[k].time - this->knownp[i].time);
				}
			}
			result.x += tempx*this->knownp[k].x;
			result.y += tempy*this->knownp[k].y;
			result.z += tempz*this->knownp[k].z;
		}
		return result;
	}

	void LarrangeInterpolate::CmpError(
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

		for (int i = 0; i < data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= this->Slope / 2  && index <= (knownp.size() - (this->Slope - this->Slope / 2) - 1))
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
	}
}