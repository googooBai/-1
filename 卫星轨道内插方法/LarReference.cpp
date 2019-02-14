//#include<stdio.h>
//
//#include<math.h>
//
//#include<stdlib.h>
//
//#define PI 3.1415926
//
//
//
//double * getX1(int n) {//�õ���һ��ڵ��X���� 
//
//	double * x = (double *)malloc((n + 1) * sizeof(double));
//
//	int i = 0;
//
//	for (i = 0; i <= n; i++) {
//
//		x[i] = -5 + 10.0 * i / (double)n;
//
//	}
//
//	return x;
//
//}
//
//double * getX2(int n) {//�õ��ڶ���ڵ��X���� 
//
//	double * x = (double *)malloc((n + 1) * sizeof(double));
//
//	int i = 0;
//
//	for (i = 0; i <= n; i++) {
//
//		//x[i] = -5 * cos((2 * i + 1)*PI / (2 * n + 2));
//		x[i] = -5 + 10*cos(i*PI / n);
//	}
//
//	return x;
//
//}
//
//double * getY() {//�õ�Ҫ�ȽϵĲ��Խڵ��X���� 
//
//	double *y = (double *)malloc(501 * sizeof(double));
//
//	int j;
//
//	for (j = 0; j <= 500; j++) {
//
//		y[j] = -5 + (double)10 * j / 500;
//
//	}
//
//	return y;
//
//}
//
//double f(double x) {//f��x�� 
//
//	return (double)1 / (1 + x*x);
//
//}
//
//double * getFX(double * x, int n) {//��һ��ڵ����룬�õ�һ���Ӧ��f��x�� 
//
//	int i;
//
//	double * fx = (double *)malloc((n + 1) * sizeof(double));
//
//	for (i = 0; i <= n; i++) {
//
//		fx[i] = f(x[i]);
//
//	}
//
//	return fx;
//
//}
//
//double lagrange(double *x, double *fx, int n, double input) {//�������ն���ʽ������һ��ڵ㣬һ���Ӧ�ڵ��fx���õ�input����������ս�� 
//
//	int i, j;
//
//	double result = 0.0;
//
//	double temp = 1;
//
//	for (i = 0; i <= n; i++) {
//
//		temp = 1;
//
//		for (j = 0; j <= n; j++) {
//
//			if (j != i) {
//
//				temp = temp * (input - x[j]) / (x[i] - x[j]);
//
//			}
//
//		}
//
//		result = result + temp * fx[i];
//
//	}
//
//	return result;
//
//}
//
//void lagrangeProcess(int flag) {
//
//	double *x5, *x10, *x20, *x40;
//
//	double *fx5, *fx10, *fx20, *fx40;
//
//	double* y;
//
//	double error5, error10, error20, error40;
//
//	int i, j;
//
//	if (flag == 1) {
//
//		x5 = getX1(5);
//
//		x10 = getX1(10);
//
//		x20 = getX1(20);
//
//		x40 = getX1(40);
//
//	}
//	else if (flag == 2) {
//
//		x5 = getX2(5);
//
//		x10 = getX2(10);
//
//		x20 = getX2(20);
//
//		x40 = getX2(40);
//
//	}
//	else {
//
//		printf("error\n");
//
//		return;
//
//	}
//
//	fx5 = getFX(x5, 5);
//
//	fx10 = getFX(x10, 10);
//
//	fx20 = getFX(x20, 20);
//
//	fx40 = getFX(x40, 40);
//
//	y = getY();
//
//	error5 = error10 = error20 = error40 = 0.0;
//
//	double abserror5, abserror10, abserror20, abserror40;
//
//	for (i = 0; i <= 500; i++) {
//
//		/*printf("f(%0.2f) = %f, while langrange5(%0.2f) = %lf\n",y[i],f(y[i]),y[i],lagrange(x5,fx5,5,y[i]));
//
//		printf("f(%0.2f) = %f, while langrange10(%0.2f) = %lf\n",y[i],f(y[i]),y[i],lagrange(x10,fx10,10,y[i]));
//
//		printf("f(%0.2f) = %f, while langrange20(%0.2f) = %lf\n",y[i],f(y[i]),y[i],lagrange(x20,fx20,20,y[i]));
//
//		printf("f(%0.2f) = %f, while langrange40(%0.2f) = %lf\n",y[i],f(y[i]),y[i],lagrange(x40,fx40,40,y[i]));*/
//
//		abserror5 = fabs(lagrange(x5, fx5, 5, y[i]) - f(y[i]));
//
//		abserror10 = fabs(lagrange(x10, fx10, 10, y[i]) - f(y[i]));
//
//		abserror20 = fabs(lagrange(x20, fx20, 20, y[i]) - f(y[i]));
//
//		abserror40 = fabs(lagrange(x40, fx40, 40, y[i]) - f(y[i]));
//
//		if (abserror5>error5)
//
//			error5 = abserror5;
//
//		if (abserror10>error10)
//
//			error10 = abserror10;
//
//		if (abserror20>error20)
//
//			error20 = abserror20;
//
//		if (abserror40>error40)
//
//			error40 = abserror40;
//
//	}
//
//	if (flag == 1) {
//
//		printf("��һ��ڵ�����\n");
//
//	}
//
//	else {
//
//		printf("�ڶ���ڵ�����\n");
//
//	}
//
//	printf("6���ڵ㣬����������Ϊ��%0.15e\n", error5);
//
//	printf("11���ڵ㣬����������Ϊ��%0.15e\n", error10);
//
//	printf("21���ڵ㣬����������Ϊ��%0.15e\n", error20);
//
//	printf("41���ڵ㣬����������Ϊ��%0.15e\n", error40);
//
//	free(x5);
//
//}
//
//
//
//int main() {
//
//	lagrangeProcess(1);
//
//	lagrangeProcess(2);
//
//	getchar();
//
//	return 0;
//
//}