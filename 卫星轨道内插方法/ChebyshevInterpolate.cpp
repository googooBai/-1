#include"ClassInterpolate.h"
#include"Function.h"

namespace ClassInter {
	ChebyshevInterpolate::ChebyshevInterpolate(const vector<Point>& data, int slope) :Slope(slope), knownp(data)
	{
		//Caculate V=TA-L
		A = new MatrixXd*[knownp.size() - m_C + 1];
		error[0] = new double[knownp.size() - m_C + 1];
		error[1] = new double[knownp.size() - m_C + 1];
		error[2] = new double[knownp.size() - m_C + 1];
		is_knownp_error = new bool[knownp.size()]{ false };
		start = new Time[knownp.size() - m_C + 1];
		end = new Time[knownp.size() - m_C + 1];
		for (int i = 0; i <= knownp.size() - m_C; i++)
		{
			
			A[i] = new MatrixXd(Slope + 1, 3);
			MatrixXd newT(m_C, Slope + 1);
			MatrixXd newPx = MatrixXd::Zero(m_C, m_C);
			MatrixXd newPy = MatrixXd::Zero(m_C, m_C);
			MatrixXd newPz = MatrixXd::Zero(m_C, m_C);
			MatrixXd newQx = MatrixXd::Zero(m_C, m_C);
			MatrixXd newQy = MatrixXd::Zero(m_C, m_C);
			MatrixXd newQz = MatrixXd::Zero(m_C, m_C);
			MatrixXd newL(m_C, 3);
			MatrixXd V;

			vector<int> knownp_index(m_C);
			for (int j = 0; j < m_C; j++)
				knownp_index[j] = i + j;
			start[i] = knownp[i].time;
			end[i] = knownp[i + m_C - 1].time;
			//Caculate T
			for (int j = i; j < i + m_C; j++)
			{
				double tao =
					2 * (knownp[j].time - start[i]) /
					(end[i] - start[i])
					- 1;
				for (int k = 0; k <= Slope; k++)
				{
					if (k == 0)
						newT(j - i, k) = 1;
					else if (k == 1)
						newT(j - i, k) = tao;
					else
						newT(j - i, k) = 2 * tao*newT(j - i, k - 1) - newT(j - i, k - 2);
				}
				//Set L
				newL(j - i, 0) = knownp[j].x;
				newL(j - i, 1) = knownp[j].y;
				newL(j - i, 2) = knownp[j].z;

				//Set P
				newPx(j - i, j - i) = 0.0001 / knownp[j].xx;
				newPy(j - i, j - i) = 0.0001 / knownp[j].yy;
				newPz(j - i, j - i) = 0.0001 / knownp[j].zz;
				//Set P.inverse()
				newQx(j - i, j - i) = knownp[j].xx;
				newQy(j - i, j - i) = knownp[j].yy;
				newQz(j - i, j - i) = knownp[j].zz;
			}

			//Test
			(*A[i]).col(0) = (newT.transpose()*newPx*newT).inverse()*newT.transpose()*newPx*newL.col(0);			
			(*A[i]).col(1) = (newT.transpose()*newPy*newT).inverse()*newT.transpose()*newPy*newL.col(1);
			(*A[i]).col(2) = (newT.transpose()*newPz*newT).inverse()*newT.transpose()*newPz*newL.col(2);

			/*ofstream fout;
			fout.open("log.txt");
			fout.setf(ios_base::fixed);
			fout.precision(15);
			fout << newT << endl << endl;
			fout << newL << endl << endl;
			fout << newPx << endl << endl;
			fout << newPy << endl << endl;
			fout << newPz << endl << endl;
			fout.close();
			cout<<i;*/

			//bool check;	//粗差探测
			//do
			//{
			//	ofstream fout;
			//	fout.open("log.txt");
			//	check = false;
			//	MatrixXd T(newT), L(newL),
			//		Px(newPx), Py(newPy), Pz(newPz),
			//		Qx(newQx), Qy(newQy), Qz(newQz);

			//	/*(*A[i]).col(0) = (T.transpose()*Px*T).inverse()*T.transpose()*Px*L.col(0);
			//	(*A[i]).col(1) = (T.transpose()*Py*T).inverse()*T.transpose()*Py*L.col(1);
			//	(*A[i]).col(2) = (T.transpose()*Pz*T).inverse()*T.transpose()*Pz*L.col(2);*/
			//	(*A[i]).col(0) = (T.transpose()*T).inverse()*T.transpose()*L.col(0);
			//	(*A[i]).col(1) = (T.transpose()*T).inverse()*T.transpose()*L.col(1);
			//	(*A[i]).col(2) = (T.transpose()*T).inverse()*T.transpose()*L.col(2);
			//	//
			//	if (T.rows() == Slope + 1)
			//		break;

			//	V = T*(*A[i]) - L;

			//	fout.setf(ios_base::fixed);
			//	fout.precision(15);
			//	fout << newT << endl << endl;
			//	fout << newL << endl << endl;
			//	fout << newPx << endl << endl;
			//	fout << newPy << endl << endl;
			//	fout << newPz << endl << endl;
			//	fout << (*A[i]) << endl << endl;
			//	fout << V << endl << endl;
			//	fout.close();
			//	//Caculate Mean Error
			//	/*MatrixXd e;
			//	e = (V.col(0)).transpose()*Px*V.col(0);
			//	error[0][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
			//	e = (V.col(1)).transpose()*Py*V.col(1);
			//	error[1][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
			//	e = (V.col(2)).transpose()*Pz*V.col(2);
			//	error[2][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));*/
			//	MatrixXd e;
			//	e = (V.col(0)).transpose()*V.col(0);
			//	error[0][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
			//	e = (V.col(1)).transpose()*V.col(1);
			//	error[1][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
			//	e = (V.col(2)).transpose()*V.col(2);
			//	error[2][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
			//	//Caculate Qvv
			//	MatrixXd Qvv_x, Qvv_y, Qvv_z;
			//	/*Qvv_x = Qx - T*(T.transpose()*Px*T).inverse()*T.transpose();
			//	Qvv_y = Qy - T*(T.transpose()*Py*T).inverse()*T.transpose();
			//	Qvv_z = Qz - T*(T.transpose()*Pz*T).inverse()*T.transpose();*/
			//	Qvv_x = T*(T.transpose()*T).inverse()*T.transpose();
			//	Qvv_y = T*(T.transpose()*T).inverse()*T.transpose();
			//	Qvv_z = T*(T.transpose()*T).inverse()*T.transpose();

			//	//Check Error
			//	int index_of_error = -1;
			//	double max_error = 0;
			//	double ex, ey, ez;
			//	for (int j = 0; j < T.rows(); j++)
			//	{
			//		ex = abs(V(j, 0) / error[0][i] / sqrt(1-Qvv_x(j, j)));
			//		ey = abs(V(j, 1) / error[1][i] / sqrt(1-Qvv_y(j, j)));
			//		ez = abs(V(j, 2) / error[2][i] / sqrt(1-Qvv_z(j, j)));
			//		//cout << ex <<" "<< ey << " " << ez << endl;
			//		if (ex > 3 || ey > 3 || ez > 3) {
			//			/*if (ex > max_error)
			//			{
			//			max_error = ex;
			//			index_of_error = j;
			//			}
			//			if (ey > max_error)
			//			{
			//			max_error = ey;
			//			index_of_error = j;
			//			}
			//			if (ez > max_error)
			//			{
			//			max_error = ez;
			//			index_of_error = j;
			//			}*/
			//			if (ex*ex + ey*ey + ez*ez > max_error)
			//			{
			//				max_error = ex*ex + ey*ey + ez*ez;
			//				index_of_error = j;
			//			}
			//			check = true;
			//		}
			//	}
			//	if (check)
			//	{
			//		//vector<Point> test;
			//		//test.clear();
			//		//for (int j = 0; j < knownp_index.size(); j++)
			//		//	test.push_back(knownp[knownp_index[j]]);
			//		//if (test.size() == 42)
			//		//{
			//		//	for (int j = 0; j < test.size(); j++)
			//		//	{
			//		//		fout.precision(15);
			//		//		fout << test[j].time << " " << test[j].x << endl;
			//		//	}
			//		//	fout << endl;
			//		//	MatrixXd MATT(42, Slope + 1);
			//		//	MatrixXd MATL(42, 3);
			//		//	MatrixXd MATPx = MatrixXd::Zero(42, 42);
			//		//	MatrixXd MATPy = MatrixXd::Zero(42, 42);
			//		//	MatrixXd MATPz = MatrixXd::Zero(42, 42);
			//		//	MatrixXd MATQx = MatrixXd::Zero(42, 42);
			//		//	MatrixXd MATQy = MatrixXd::Zero(42, 42);
			//		//	MatrixXd MATQz = MatrixXd::Zero(42, 42);
			//		//	for (int j = 0; j < 42; j++)
			//		//	{
			//		//		double tao1 =
			//		//			2 * (test[j].time - start) /
			//		//			(end - start)
			//		//			- 1;
			//		//		for (int k = 0; k <= Slope; k++)
			//		//		{
			//		//			if (k == 0)
			//		//				MATT(j, k) = 1;
			//		//			else if (k == 1)
			//		//				MATT(j, k) = tao1;
			//		//			else
			//		//				MATT(j, k) = 2 * tao1*MATT(j, k - 1) - MATT(j, k - 2);
			//		//		}
			//		//		//Set L
			//		//		MATL(j, 0) = test[j].x;
			//		//		MATL(j, 1) = test[j].y;
			//		//		MATL(j, 2) = test[j].z;

			//		//		//Set P
			//		//		MATPx(j, j) = 1 / test[j].xx;
			//		//		MATPy(j, j) = 1 / test[j].yy;
			//		//		MATPz(j, j) = 1 / test[j].zz;
			//		//		//Set P.inverse()
			//		//		MATQx(j, j) = test[j].xx;
			//		//		MATQy(j, j) = test[j].yy;
			//		//		MATQz(j, j) = test[j].zz;
			//		//	}
			//		//	MatrixXd MATA;
			//		//	MATA = (MATT.transpose()*MATPx*MATT).inverse()*MATT.transpose()*MATPx*MATL.col(0);
			//		//	fout.precision(15);
			//		//	fout << MATA;
			//		//}
			//		is_knownp_error[knownp_index[index_of_error]] = true;
			//		for (int j = index_of_error; j != knownp_index.size() - 1; j++)
			//		{
			//			knownp_index[j] = knownp_index[j + 1];
			//		}
			//		knownp_index.pop_back();
			//		start[i] = knownp[knownp_index.front()].time;
			//		end[i] = knownp[knownp_index.back()].time;					
			//		
			//		newT = MatrixXd(T.rows() - 1, Slope + 1);
			//		newL = MatrixXd(T.rows() - 1, 3);
			//		newPx = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		newPy = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		newPz = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		newQx = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		newQy = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		newQz = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
			//		for (int j = 0; j != index_of_error; j++)
			//		{
			//			//newT.row(j) = T.row(j);
			//			newL.row(j) = L.row(j);
			//			newPx(j, j) = Px(j, j);
			//			newPy(j, j) = Py(j, j);
			//			newPz(j, j) = Pz(j, j);
			//			newQx(j, j) = Qx(j, j);
			//			newQy(j, j) = Qy(j, j);
			//			newQz(j, j) = Qz(j, j);
			//		}
			//		for (int j = index_of_error + 1; j < T.rows(); j++)
			//		{
			//			//newT.row(j - 1) = T.row(j);
			//			newL.row(j - 1) = L.row(j);
			//			newPx(j - 1, j - 1) = Px(j, j);
			//			newPy(j - 1, j - 1) = Py(j, j);
			//			newPz(j - 1, j - 1) = Pz(j, j);
			//			newQx(j-1, j-1) = Qx(j, j);
			//			newQy(j-1, j-1) = Qy(j, j);
			//			newQz(j-1, j-1) = Qz(j, j);
			//		}
			//		for (int j = 0; j < knownp_index.size(); j++)
			//		{
			//			double tao =
			//				2 * (knownp[knownp_index[j]].time - start[i]) /
			//				(end[i] - start[i])
			//				- 1;
			//			for (int k = 0; k <= Slope; k++)
			//			{
			//				if (k == 0)
			//					newT(j, k) = 1;
			//				else if (k == 1)
			//					newT(j, k) = tao;
			//				else
			//					newT(j, k) = 2 * tao*newT(j, k - 1) - newT(j, k - 2);
			//			}
			//		}
			//	}
			//	
			//} while (check);
		}
	}
	//ChebyshevInterpolate::ChebyshevInterpolate(const vector<Point>& data, int slope) :Slope(slope), knownp(data)
	//{
	//	//Caculate V=TA-L
	//	A = new MatrixXd*[knownp.size() - m_C + 1];
	//	error[0] = new double[knownp.size() - m_C + 1];
	//	error[1] = new double[knownp.size() - m_C + 1];
	//	error[2] = new double[knownp.size() - m_C + 1];
	//	is_knownp_error = new bool[knownp.size()]{ false };
	//	for (int i = 0; i <= knownp.size() - m_C; i++)
	//	{
	//		int obvation = m_C;
	//		for (int j = i; j < i + m_C; j++)
	//		{
	//			obvation = is_knownp_error[j] ? obvation-1 : obvation;
	//		}
	//		A[i] = new MatrixXd(Slope + 1, 3);
	//		MatrixXd newT(obvation, Slope + 1);
	//		MatrixXd newPx = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newPy = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newPz = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newQx = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newQy = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newQz = MatrixXd::Zero(obvation, obvation);
	//		MatrixXd newL(obvation, 3);
	//		MatrixXd V;

	//		vector<int> knownp_index(obvation);

	//		//Caculate T
	//		int cur = 0;
	//		for (int j = i; j < i + m_C; j++)
	//		{
	//			if (!is_knownp_error[j])
	//			{
	//				double tao =
	//					2 * (knownp[j].time - knownp[i].time) /
	//					(knownp[i + m_C - 1].time - knownp[i].time)
	//					- 1;
	//				for (int k = 0; k <= Slope; k++)
	//				{
	//					if (k == 0)
	//						newT(cur, k) = 1;
	//					else if (k == 1)
	//						newT(cur, k) = tao;
	//					else
	//						newT(cur, k) = 2 * tao*newT(cur, k - 1) -newT(cur, k - 2);
	//				}
	//				//Set L
	//				newL(cur, 0) = knownp[j].x;
	//				newL(cur, 1) = knownp[j].y;
	//				newL(cur, 2) = knownp[j].z;

	//				//Set P
	//				newPx(cur, cur) = 1 / knownp[j].xx;
	//				newPy(cur, cur) = 1 / knownp[j].yy;
	//				newPz(cur, cur) = 1 / knownp[j].zz;
	//				//Set P.inverse()
	//				newQx(cur, cur) =  knownp[j].xx;
	//				newQy(cur, cur) =  knownp[j].yy;
	//				newQz(cur, cur) =  knownp[j].zz;
	//				knownp_index[cur] = j;
	//				cur++;
	//			}					
	//		}
	//		//test if make a obvation*obvation matrix
	//		if (cur != obvation)
	//			cout << "no ways!"<<endl;

	//		//Test
	//		/*(*A[i]).col(0) = (newT.transpose()*newPx*newT).inverse()*newT.transpose()*newPx*newL.col(0);
	//		(*A[i]).col(1) = (newT.transpose()*newPy*newT).inverse()*newT.transpose()*newPy*newL.col(1);
	//		(*A[i]).col(2) = (newT.transpose()*newPz*newT).inverse()*newT.transpose()*newPz*newL.col(2);*/
	//		bool check;	//粗差探测
	//		do
	//		{
	//			check = false;
	//			MatrixXd T(newT), L(newL),
	//				Px(newPx), Py(newPy), Pz(newPz),
	//				Qx(newQx), Qy(newQy), Qz(newQz);

	//			(*A[i]).col(0) = (T.transpose()*Px*T).inverse()*T.transpose()*Px*L.col(0);
	//			(*A[i]).col(1) = (T.transpose()*Py*T).inverse()*T.transpose()*Py*L.col(1);
	//			(*A[i]).col(2) = (T.transpose()*Pz*T).inverse()*T.transpose()*Pz*L.col(2);
	//			V = T*(*A[i]) - L;
	//			cout << i << endl;
	//			cout << V << endl<<endl;
	//			//Caculate Mean Error
	//			MatrixXd e;
	//			e = (V.col(0)).transpose()*Px*V.col(0);
	//			error[0][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
	//			e = (V.col(1)).transpose()*Py*V.col(1);
	//			error[1][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
	//			e = (V.col(2)).transpose()*Pz*V.col(2);
	//			error[2][i] = sqrt(e(0, 0) / (T.rows() - Slope - 1));
	//			//Caculate Qvv
	//			MatrixXd Qvv_x, Qvv_y, Qvv_z;
	//			Qvv_x = Qx - T*(T.transpose()*Px*T).inverse()*T.transpose();
	//			Qvv_y = Qy - T*(T.transpose()*Py*T).inverse()*T.transpose();
	//			Qvv_z = Qz - T*(T.transpose()*Pz*T).inverse()*T.transpose();
	//			//Check Error
	//			int index_of_error = -1;
	//			double max_error = 0;
	//			double ex, ey, ez;
	//			for (int j = 0; j < T.rows(); j++)
	//			{
	//				ex = abs(V(j, 0) / error[0][i] / sqrt(Qvv_x(j, j)));
	//				ey = abs(V(j, 1) / error[1][i] / sqrt(Qvv_y(j, j)));
	//				ez = abs(V(j, 2) / error[2][i] / sqrt(Qvv_z(j, j)));
	//				if (ex > 3 || ey > 3 || ez > 3) {
	//					/*if (ex > max_error)
	//					{
	//					max_error = ex;
	//					index_of_error = j;
	//					}
	//					if (ey > max_error)
	//					{
	//					max_error = ey;
	//					index_of_error = j;
	//					}
	//					if (ez > max_error)
	//					{
	//					max_error = ez;
	//					index_of_error = j;
	//					}*/
	//					if (ex*ex + ey*ey + ez*ez > max_error)
	//					{
	//						max_error = ex*ex + ey*ey + ez*ez;
	//						index_of_error = j;
	//					}
	//					check = true;
	//				}
	//			}
	//			if (check)
	//			{
	//				is_knownp_error[knownp_index[index_of_error]] = true;
	//				for (int j = index_of_error; j != knownp_index.size() - 1; j++)
	//				{
	//					knownp_index[j] = knownp_index[j + 1];
	//				}
	//				knownp_index.pop_back();
	//				newT = MatrixXd(T.rows() - 1, Slope + 1);
	//				newL = MatrixXd(T.rows() - 1, 3);
	//				newPx = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				newPy = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				newPz = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				newQx = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				newQy = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				newQz = MatrixXd::Zero(T.rows() - 1, T.rows() - 1);
	//				for (int j = 0; j != index_of_error; j++)
	//				{
	//					newT.row(j) = T.row(j);
	//					newL.row(j) = L.row(j);
	//					newPx(j, j) = Px(j, j);
	//					newPy(j, j) = Py(j, j);
	//					newPz(j, j) = Pz(j, j);
	//					newQx(j, j) = Qx(j, j);
	//					newQy(j, j) = Qy(j, j);
	//					newQz(j, j) = Qz(j, j);
	//				}
	//				for (int j = index_of_error + 1; j < T.rows(); j++)
	//				{
	//					newT.row(j - 1) = T.row(j);
	//					newL.row(j - 1) = L.row(j);
	//					newPx(j - 1, j - 1) = Px(j, j);
	//					newPy(j - 1, j - 1) = Py(j, j);
	//					newPz(j - 1, j - 1) = Pz(j, j);
	//					newQx(j - 1, j - 1) = Qx(j, j);
	//					newQy(j - 1, j - 1) = Qy(j, j);
	//					newQz(j - 1, j - 1) = Qz(j, j);
	//				}
	//			}
	//		} while (check);
	//	}
	//}

	Point ChebyshevInterpolate::interpolate(Time t,ostream& os)
	{
		Point res;
		res.time = t;

		int posfirst = 0, posend = this->knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, this->knownp);

		//内插起始点与终结点
		if (index < m_C/2)
			index = 0;
		else if (index >this->knownp.size()-1 - (m_C - m_C / 2 - 1))
			index = this->knownp.size() - m_C;
		else
			index = index - m_C/2;

		double tao =
			2 * (t - this->knownp[index].time) /
			(this->knownp[index + m_C-1].time - this->knownp[index].time)
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
		for (int j = index; j <= index + m_C; j++)
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

	Point ChebyshevInterpolate::interpolate(Time t)
	{
		Point res;
		res.time = t;

		int posfirst = 0, posend = this->knownp.size() - 1;
		int index = FindIndex(posfirst, posend, t, this->knownp);

		//内插起始点与终结点
		if (index < m_C / 2)
			index = 0;
		else if (index >this->knownp.size() - 1 - (m_C - m_C / 2 - 1))
			index = this->knownp.size() - m_C;
		else
			index = index - m_C / 2;

		double tao =
			2 * (t - this->knownp[index].time) /
			(this->knownp[index + m_C - 1].time - this->knownp[index].time)
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
		ofstream fout;
		fout.open("tao.txt");
		fout.setf(ios_base::fixed);
		fout.precision(15);
		fout << MatT<<endl<<endl;
		fout << coor << endl << endl;
		fout.close();
		return res;
	}


	void ChebyshevInterpolate::CmpError(
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

		for (int i = 0; i<data.size(); i++)
		{
			data[i].is_effect_by_knownp = false;
			data[i].is_g = false;
			data[i].is_d = false;
			data[i].is_big_error = false;
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
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= m_C / 2 && index <= this->knownp.size()-1 - (m_C - m_C / 2 - 1))
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
			/*if (data[i].is_knownp)
			{
				if (this->is_knownp_error[knownp_index++])
					data[i].is_g = true;
				else
					data[i].is_g = false;
			}*/
			int index = FindIndex(0, knownp.size() - 1, data[i].time_d, knownp);
			if (index >= m_C / 2 && index <= this->knownp.size()-1 - (m_C -m_C/2-1))
			{
				index = index - m_C / 2;
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

				if (!data[i].is_knownp)
				{
					if (abs(data[i].error_x) > 2 * RMSE_X ||
						abs(data[i].error_y) > 2 * RMSE_Y ||
						abs(data[i].error_z) > 2 * RMSE_Z)
					/*if (abs(data[i].error_x) >  this->error[0][index] * 3 ||
						abs(data[i].error_y) >  this->error[1][index] * 3 ||
						abs(data[i].error_z) >  this->error[2][index] * 3)*/
						data[i].is_g = true;
					else
						data[i].is_g = false;
				}


				if (data[i].is_g&&data[i].is_d)
					hit++;
				else if (data[i].is_g && !data[i].is_d)
					mistake++;

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
	/*	os4 << "Amount of error: " << total << endl;
		os4 << "Amount of big error(more than 1cm): " << big_error << "	  eliminate: " << hit_big_error << endl;
		os4.setf(ios_base::fixed);
		os4.precision(1);
		os4 << "Rate of sucessful: %" << (double)hit / total * 100 << endl;
		os4 << "Rate of mistake: %" << (double)mistake / data.size() * 100 << endl;
		os4.precision(3);
		os4
			<< setw(8) << RMSE_X
			<< setw(8) << RMSE_Y
			<< setw(8) << RMSE_Z
			<< endl;*/
		

		/*total = 0;
		hit = 0;
		mistake = 0;
		for (int i = 0; i != data.size(); i++)
		{
			if (data[i].is_knownp)
			{
				os4.setf(ios_base::fixed);
				os4.precision(15);
				os4 << setw(22) << data[i].time_d;

				if (data[i].is_g)
					os4 << setw(15) << 1;
				else
					os4 << setw(15) << 0;

				if (data[i].is_d)
				{
					total++;
					os4 << setw(15) << 1 << endl;
				}
				else
					os4 << setw(15) << 0 << endl;

				if (data[i].is_g&&data[i].is_d)
					hit++;
				if (data[i].is_g && !data[i].is_d)
					mistake++;
			}
		}*/
		/*std::cout << "Amount of known point: " << total << endl;
		std::cout << "Amount of hit known point: " << hit << endl;
		std::cout << "Amount of mistake konwn point: " << mistake << endl;*/
	}
}