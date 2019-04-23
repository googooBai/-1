#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	void CmpError(Interpolate& inter,
		int staff,
		vector<SourceData>& data,
		const vector<Point>& knownp,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2,
		ostream& os3)
	{
		double error_x = 0.0, error_y = 0.0, error_z = 0.0;
		double RMSE_X = 0.0, RMSE_Y = 0.0, RMSE_Z = 0.0;
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
			<< setw(15) << "error x"
			<< setw(15) << "error y"
			<< setw(15) << "error z"
			<< endl;
		for (int i=0;i<data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= staff / 2 && index <= (knownp.size() - (staff - staff / 2) - 1))
			{
				Point temp = inter.interpolate(data[i].time_d);
				data[i].ix = temp.x;
				data[i].iy = temp.y;
				data[i].iz = temp.z;

				int ob_index;
				{
					double t = data[i].time_d;
					int fp = 0, ep = ob.size()-1;
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
				error_x = (temp.x - data[i].x);
				data[i].error_x = error_x;
				RMSE_X += error_x*error_x;
				//y
				error_y = (temp.y - data[i].y);
				data[i].error_y = error_y;
				RMSE_Y += error_y*error_y;
				//z
				error_z = (temp.z - data[i].z);
				data[i].error_z = error_z;
				RMSE_Z += error_z*error_z;

				os.setf(ios_base::fixed);
				os.precision(15);
				os << setw(22) << data[i].time_d;
				os.precision(3);
				os
					<< setw(15) << error_x
					<< setw(15) << error_y
					<< setw(15) << error_z
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

		RMSE_X = sqrt(RMSE_X / data.size());
		RMSE_Y = sqrt(RMSE_Y / data.size());
		RMSE_Z = sqrt(RMSE_Z / data.size());
		
		//´Ö²î´¦Àí
		int total{ 0 };
		int mistake{ 0 };
		int hit{ 0 };
		int knownp_index = 0;
		for (int i = 0; i < data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= staff / 2  && index <= (knownp.size() - (staff - staff / 2) - 1))
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

				if (abs(data[i].error_x) > 2* RMSE_X ||
					abs(data[i].error_y) > 2* RMSE_Y ||
					abs(data[i].error_z) > 2* RMSE_Z)
					data[i].is_g = true;
				else
					data[i].is_g = false;


				if (data[i].is_g&&data[i].is_d)
					hit++;
				else if (data[i].is_g && !data[i].is_d)
					mistake++;

				if (!data[i].is_g&&data[i].is_d)
				{
					os3.setf(ios_base::fixed);
					os3.precision(15);
					os3 << setw(22) << data[i].time_d;
					os3.precision(3);
					os3
						<< setw(15) << data[i].error_x
						<< setw(15) << data[i].error_y
						<< setw(15) << data[i].error_z
						<< endl;
				}
			}			
			
		}
		std::cout << "Amount of error: " << total << endl;
		std::cout.setf(ios_base::fixed);
		std::cout.precision(1);
		std::cout << "Rate of sucessful: %"<<(double)hit/total*100<<endl;
		std::cout << "Rate of mistake: %" << (double)mistake / data.size() * 100 << endl;

		/*os
			<< setw(22) << "Total"
			<< setw(15) << "RMSE_X"
			<< setw(15) << "RMSE_Y"
			<< setw(15) << "RMSE_Z"
			<< setw(15) << "RMSE_3D"
			<< endl;

		
		os.precision(3);
		os
			<< setw(22) << " "
			<< setw(15) << RMSE_X
			<< setw(15) << RMSE_Y
			<< setw(15) << RMSE_Z
			<< setw(15) << RMSE_3D
			<< endl;
	*/
	}
}

