#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {

	void CmpError(Interpolate& inter,
		int staff,
		vector<SourceData>& data,
		const vector<Point>& knownp,
		const vector<Point>& ob,
		ostream& os,
		ostream& os2)
	{
		double error_x = 0.0, error_y = 0.0, error_z = 0.0;
		double RMSE_X = 0.0, RMSE_Y = 0.0, RMSE_Z = 0.0;
		double D_X = 0.0, D_Y = 0.0, D_Z = 0.0;
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
		for (int i=0;i<data.size(); i++)
		{
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= staff / 2 && index <= (knownp.size() - (staff - staff / 2) - 1))
			{
				Point temp = inter.interpolate(data[i].time_d);
				data[i].ix = temp.x;
				data[i].iy = temp.y;
				data[i].iz = temp.z;

				data[i].is_d = false;
				data[i].is_g = false;

				data[i].de_x = ob[2 * i].x - data[i].ix;
				data[i].de_y = ob[2 * i].y - data[i].iy;
				data[i].de_z = ob[2 * i].z - data[i].iz;
				//D_X += data[i].de_x*data[i].de_x;
				//D_Y += data[i].de_y*data[i].de_y;
				//D_Z += data[i].de_z*data[i].de_z;

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

			}
		}

		RMSE_X = sqrt(RMSE_X / data.size());
		RMSE_Y = sqrt(RMSE_Y / data.size());
		RMSE_Z = sqrt(RMSE_Z / data.size());
		//D_X = sqrt(D_X / data.size());
		//D_Y = sqrt(D_Y / data.size());
		//D_Z = sqrt(D_Z / data.size());


		//´Ö²î´¦Àí
		int total{ 0 };
		int mistake{ 0 };
		int hit{ 0 };
		for (int i = 0; i < data.size(); i++)
		{
			if (abs(data[i].error_x) > 2 * RMSE_X ||
				abs(data[i].error_y) > 2 * RMSE_Y ||
				abs(data[i].error_z) > 2 * RMSE_Z)
				data[i].is_g = true;
			else
				data[i].is_g = false;

			if (abs(data[i].de_x) > 2*0.03 ||
				abs(data[i].de_y) > 2*0.03 ||
				abs(data[i].de_z) > 2*0.03)
			{
				data[i].is_d = true;
				total++;
			}
			else
				data[i].is_d = false;

			if (data[i].is_g&&data[i].is_d)
				hit++;
			else if (data[i].is_g && !data[i].is_d)
				mistake++;
			
			if (!data[i].is_g)
			{
				os2.setf(ios_base::fixed);
				os2.precision(15);
				os2 << setw(22) << data[i].time_d;
				os2.precision(3);
				os2
					<< setw(15) << data[i].error_x
					<< setw(15) << data[i].error_y
					<< setw(15) << data[i].error_z
					<< endl;
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

