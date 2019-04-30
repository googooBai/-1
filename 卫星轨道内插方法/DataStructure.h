#ifndef DATASTUCTURE_H
#define DATASTUCTURE_H
namespace ClassInter {

	using Time = double;
	const Time minute = 1.0 / (24 * 60);	//����ת��������
	const int m_C = 30;	//�б�ѩ���������
	const int m_P = 30;	//N�׶���ʽ��������
	struct Point {	//�ڵ���������
		//��Ԫ
		Time time;
		//����
		double x;
		double y;
		double z;
		//Э����
		double xx;
		double yy;
		double zz;
	};

	struct SourceData
	{
		//����
		double x;
		double y;
		double z;
		//�ڲ�ֵ
		double ix;
		double iy;
		double iz;
		//Э����
		double xx;
		double yy;
		double zz;
		double xy;
		double xz;
		double yz;
		//��Ԫ
		Time time_g;	//����ѧ�����Ԫ����λ���죩
		Time time_d;	//����ѧ�����Ԫ����λ���룩
		//�ڲ�ֵ�뼸��ѧ����ϲ�
		double error_x;	
		double error_y;	
		double error_z;
		//����ѧ����붯��ѧ����ϲ�
		double de_x;
		double de_y;
		double de_z;
		//Flag
		bool is_g;	//����ѧ����Ƿ���
		bool is_d;	//����ѧ����Ƿ���
		bool is_knownp;	//�Ƿ���Ϊ��֪��
	};
}
#endif // !DATASTUCTURE_H

