#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<math.h>
#include<conio.h>
int count_elem;			//кол-во конечных элементов
int left_usl;			//тип кр у-я на левой границе
int right_usl;			//тип кр у-я на правой границе
int N;				//размерность СЛАУ
int *ig;
double *di;
double *gg;
double *f;
double *setka;			//узлы сетки 
double *v;				//решение
double gamma = 1.;
double beta1 = -1.;
double beta2 = 2.;
double lambda(double x);
void clear_mat();			//функция обнуляющая СЛАУ
FILE *out = fopen("output.txt", "w+");
double* q_last; //массив весов с предыдущего временого слоя
double* q_new; //массив узлов, с нового временого слоя
double* q_temp; //массив узлов, для промежуточных итераций
double* q_min; //массив, используеммый для минимизации невязки
double* q_min_rp; //массив, используеммый для минимизации невязки хранящий правую часть матрицы
double* q_min_R; //массив, используемый для минимизации - вектор невзяки. Так же использутеся для определения невязки, при выходе из процесса решения СНУ

double t_now, h_t, k_t, eps, t_new; //информация о времени



void matr_M();
void matr_G();
void F_vect();
void gen_mat();
int razlo_LLt();
int gauss();
int output();


double lambda(double x) //коэффициент лямбда
{
	return 1;
}

double func(double x, double t)	//возвращает значение функции правой части в заданной точке
{
	/*if (x < 2) return 0.;
	if (x < 3) return 1.;
	if (x < 7) return 0.25;*/
	//return -4.*x+x*x;
	return x;
	//return -2;
	//return x-2.;
}

double u_toch(double x, double t)	//возвращает аналитическое значение функции
{
	/*if (x <= 2) return x + 1.;
	if (x <= 3) return (-0.5*(x - 3.)*(x-3.)+3.5);
	if (x <= 7) return(-(1./8.)*(x - 3.)*(x - 3.) + 3.5);*/
	return x*t;
	//return x*t*x;
}

double du_toch(double x, double t)	//возвращает производную аналитическое значение функции
{
	//return (-0.25*(x - 3.));
	return t;
	//return t*2.*x;
}

void usl(double t)
{
	if (left_usl == 1)
	{
		di[0] = 1e+8;
		f[0] = 1e+8 * u_toch(setka[0], t);
	}
	else if (left_usl == 2)
	{
		f[0] += -lambda(setka[0])*du_toch(setka[0], t);
	}
	else
	{
		di[0] += beta1;
		f[0] += setka[0] * setka[0] * beta1*u_toch(setka[0], t);
	}
	if (right_usl == 1)
	{
		di[N - 1] = 1e+8;
		f[N - 1] = 1e+8 * u_toch(setka[count_elem], t);
	}
	else if (right_usl == 2)
		f[N - 1] += lambda(setka[count_elem])*du_toch(setka[count_elem], t);
	else
	{
		di[N - 1] += beta2;
		f[N - 1] += beta2*u_toch(setka[count_elem], t);
	}
}


void gen_mat()			//функция выделяющая память под СЛАУ и задающая портрет
{					//и считывающая сетку из файла
	FILE *in;
	int i, j;
	in = fopen("Text.txt", "r");
	fscanf(in, "%i", &count_elem);
	N = 2 * count_elem + 1;		//вычисляем размерность СЛАУ

	ig = new int[N + 1];		//выделяем память под массив указателей на строки

								//заполняем массив указателей на строки
								ig[0] = 1;
								ig[1] = 1;
								ig[2] = 2;
								ig[3] = 4;
								ig[4] = 5;
								ig[5] = 7;
								ig[6] = 8;
								ig[7] = 10;

	//ig[0] = 1;
	//ig[1] = 1;
	//ig[2] = 2;
	//ig[3] = 4;
	//ig[4] = 7;
	//ig[5] = 8;
	//ig[6] = 10;
	//ig[7] = 13;
	//ig[8] = 14;
	//ig[9] = 16;
	//ig[10] = 19;


	di = new double[N];
	f = new double[N];
	gg = new double[ig[N] - 1];
	v = new double[N];

	fscanf(in, "%i", &left_usl);
	fscanf(in, "%i", &right_usl);
	fclose(in);

	q_last = new double[N];
	q_new = new double[N];
	q_temp = new double[N];

	in = fopen("Text1.txt", "r");
	setka = new double[count_elem + 1];

	for (i = 0; i <= count_elem; i++)
	{
		fscanf(in, "%lf ", &setka[i]);
	}

	fclose(in);

	for (i = 0; i < N; i++)
	{
		q_last[i] = 0;
	}
	in = fopen("time.txt", "r");
	fscanf(in, "%lf %lf %lf %d %lf", &t_now, &h_t, &k_t);

}


void clear_mat()		//функция очищающая СЛАУ
{
	int i;

	//очистка вектора правой части и диагонали матрицы
	for (i = 0; i<N; i++)
	{
		f[i] = 0;
		di[i] = 0;
	}

	//очистка нижнего треугольника
	for (i = 0; i<ig[N]; i++)
		gg[i] = 0;

	return;
}


int LLT_decomposition()
{
	int i0, i1, j, ki, kj, d;
	double sum_l, sum_d;

	for (int i = 0; i < N; i++)
	{
		i0 = ig[i]; //индекс первого элемента строки
		i1 = ig[i + 1]; //индекс первого элемента следующей строки
		j = i - (i1 - i0); //номер столбца
		sum_d = 0; //обнуление суммы для элементов d(i)
		for (int k = i0; k < i1; k++, j++) //пока не пройдены все элементы строки
		{
			sum_l = 0; //обнуление суммs для элементов l(ij)
			ki = ig[i]; //индекс первого элемента в i-той строке
			kj = ig[j]; //индекс первого элемента в j-том столбце
			d = k - ki - (ig[j + 1] - ig[j]);

			if (d < 0) //если количество элементов в столбце больше числа просмотренных элементов в i-той строке
				kj += abs(d); //смещаем индекс элемента в j-том столбце на разницу количества элементов в строке и столбце
			else ki += d; //смещаем индекс элемента в i-той строке на разницу
			for (; ki < k; ki++, kj++)
			{
				sum_l += gg[ki] * gg[kj];	//суммирование произведения элементов l(ik)*l(kj) для элементов l(ij)
			}
			gg[k] = (gg[k] - sum_l) / di[j]; //вычисление l(ij) по формуле
			sum_d += gg[k] * gg[k]; //суммирование произведения элементов l(ik)*l(ik) с целью вычислить диагональные элементы
		}
		//условие разложимости матрицы в LU(sq) разложение
		if (di[i] - sum_d <0 || i>0 && di[i] == 0) //если диагональные элементы матрицы L отрицательные (корень) или они нулевые (знаменатель) 
		{
			printf("Ошибка: разложить матрицу нельзя");
			return 0;
		}
		else di[i] = sqrt(di[i] - sum_d); //вычисление диагональных элементов по формуле
	}
	return 1;
}


int output()		//вывод результата в файл
{
	int i;
	fprintf(out, "\n");
	double sum1 = 0, sum2 = 0;
	for (i = 0; i <= count_elem; i++) { sum1 += pow(f[2 * i] - u_toch(setka[i], t_new), 2); sum2 += pow(u_toch(setka[i], t_new), 2); }

	for (i = 0; i <= count_elem; i++)
		fprintf(out, "%2.4lf\t%2.7le\t%2.7le\t%2.7le\n", setka[i], f[2 * i], u_toch(setka[i], t_new), abs(f[2 * i] - u_toch(setka[i], t_new)));
	fprintf(out, "%.2le", sqrt(sum1) / sqrt(sum2));
	return 0;
}

void Ly_is_b()
{
	int i, j, k, i0, i1;
	double sum = 0;
	//проход по всем строкам
	for (i = 0; i < N; i++)
	{
		sum = 0;
		//индекс первого элемента строки
		i0 = ig[i];
		//индекс первого элемента следующей строки
		i1 = ig[i + 1];
		//номер столбца
		j = i - (i1 - i0);
		for (k = i0; k < i1; k++, j++)
			//суммирование произведения элементов l(ij)*y(i)
			sum += gg[k] * f[j];
		f[i] = f[i] - sum;
		f[i] /= di[i];
	}
}

//Oбратный ход
void Ux_is_y()
{
	int i, j, k;
	int i_last, num;
	for (i = N - 1; i >= 0; i--)
	{
		//деление на диагональный элемент
		f[i] /= di[i];
		//индекс последнего элемента следующей строки
		i_last = ig[i + 1] - 1;
		//количество элементов в i-той строке
		num = ig[i + 1] - ig[i];
		//номер столбца в котором находится первый элемент строки
		j = i - num - 1;
		//пока не пройдены все элементы столбца
		for (k = i - 1; k > j; k--, i_last--)
			f[k] -= gg[i_last] * f[i];
	}
}

void main()
{
	gen_mat();
	for (int i = 0; i < 10; i++) {
		for (int j = 0; j < N; j++)
			q_new[i] = q_last[i];
		t_new = t_now + h_t;
		clear_mat();
		matr_M();
		matr_G();
		F_vect();
		usl(t_new);
		//LU_decomposition();
		//Ly_is_b();
		//Ux_is_y();
		if (!LLT_decomposition()) {
			printf("Error razlo_LLt\n");
			_getch();
			return;
		}
		Ly_is_b();
		Ux_is_y();
		//if (gauss())	return;	
		output();
		for (int k = 0; k < N; k++)
			q_new[k] = f[k];
		double* ch_pt;
		ch_pt = q_new;
		q_new = q_last;
		q_last = ch_pt;

		t_now = t_new;
		h_t *= k_t;
	}
	return;
}

//void Mq() {
//	umnozh M na q_last;
//	int i;
//		for (i = 0; i<count_elem; i++)
//		{
//
//			di[2 * i] += h*(4. / 30.)*gamma / h_t;
//			di[2 * i + 1] += h*(16. / 30.)*gamma / h_t;
//			di[2 * i + 2] += h*(4. / 30.)*gamma / h_t;
//
//			gg[ig[2 * i + 1]] += h*(2. / 30.)*gamma / h_t;
//			gg[ig[2 * i + 2]] += h*(-1. / 30.)*gamma / h_t;
//			gg[ig[2 * i + 2] + 1] += h*(2. / 30.)*gamma / h_t;
//		}
//}

void matr_M()		//добавляет в глобальную матрицу матрицу масс
{
	int i;
	double h;
	double M11, M22, M33, M12, M13, M23;
	for (i = 0; i<count_elem; i++)
	{
		//вычисляем шаг
		h = setka[i + 1] - setka[i];

		//вычисляем матрицу массы
		M11 = h*(4. / 30.)*gamma / h_t;
		M22 = h*(16. / 30.)*gamma / h_t;
		M33 = h*(4. / 30.)*gamma / h_t;
			
		di[2 * i] += M11;
		di[2 * i + 1] += M22;
		di[2 * i + 2] += M33;

		gg[ig[2 * i + 1]] += h*(2. / 30.)*gamma / h_t;
		gg[ig[2 * i + 2]] += h*(-1. / 30.)*gamma / h_t;
		gg[ig[2 * i + 2] + 1] += h*(2. / 30.)*gamma / h_t;
	}

	return;
}

void F_vect()//генерирует вектор правой части, массив должен быть чистым
{
	int i;
	double h, a, b, c, d;
	double c11, c22, c33, c44, c21, c31, c32, c41, c42, c43;
	double M11, M22, M33, M12, M13, M23;
	for (i = 0; i<count_elem; i++)
	{
		//вычисляем шаг
		h = setka[i + 1] - setka[i];
		M11 = h*(4. / 30.)*gamma / h_t;
		M22 = h*(16. / 30.)*gamma / h_t;
		M33 = h*(4. / 30.)*gamma / h_t;
		M13 = h*(-1. / 30.)*gamma / h_t;
		M12 = h*(2. / 30.)*gamma / h_t;
		M23 = h*(2. / 30.)*gamma / h_t;
		//вычисляем вектор правой части
		a = func(setka[i], t_new);
		b = func(setka[i] + h / 2., t_new);
		c = func(setka[i + 1] - 0.1, t_new);

		f[2 * i] += h*(4.*a + 2.*b + -1.*c) / 30. + (M11*q_last[0]+M12*q_last[1]+M13*q_last[2]);
		f[2 * i + 1] += h*(2.*a + 16.*b + 2.*c) / 30. + (M12*q_last[0] + M22*q_last[1] + M23*q_last[2]);
		f[2 * i + 2] += h*(-1.*a + 2.*b + 4.*c) / 30. + (M13*q_last[0] + M23*q_last[1] + M33*q_last[2]);
	}
	return;
}

void matr_G()	//добавляет в глобальную матрицу матрицу жесткости
{
	int i;
	double h;
	double lambda1, lambda2;

	for (i = 0; i<count_elem; i++)
	{
		//вычисляем шаг
		h = setka[i + 1] - setka[i];
		lambda1 = lambda(setka[i]);
		lambda2 = lambda(setka[i + 1]);

		//di[2 * i] += (lambda1*(7.0)) /( h* 3.0);
		//di[2 * i + 1] += (lambda1*(16.0)) / (h* 3.0);
		//di[2 * i + 2] += (lambda1*(7.0)) / (h* 3.0);

		di[2 * i] += (lambda1*(11.) + lambda2*(3.)) / (h* 6.0);
		di[2 * i + 1] += (lambda1*(16.) + lambda2*(16.)) / (h* 6.0);
		di[2 * i + 2] += (lambda1*(3.) + lambda2*(11.)) / (h* 6.0);

		gg[ig[2 * i + 1]] += (lambda1*(-12.) + lambda2*(-4.)) / (h* 6.0);
		gg[ig[2 * i + 2]] += (lambda1*(1.) + lambda2*(1.)) / (h* 6.0);
		gg[ig[2 * i + 2] + 1] += (lambda1*(-4.) + lambda2*(-12.)) / (h* 6.0);

		//gg[ig[2 * i + 1] - 2] += (lambda1*(-8.0) )/ (h* 3.0);
		//gg[ig[2 * i + 1] - 1] +=( lambda1*(1.0)) / (h* 3.0);
		////		gg[ig[2 * i + 2] - 1] += h*(38.0*setka[i] * h + 38.0*setka[i] * setka[i] + 15.0*h*h) / 3360.0*gamma;
		//gg[ig[2 * i + 1]] += (lambda1*(-8.0)) / (h* 3.0);
	}

	return;
}

//void matr_G()	//добавляет в глобальную матрицу матрицу жесткости
//{
//	int i;
//	double h;
//	double lambda1, lambda2;
//
//	for (i = 0; i<count_elem; i++)
//	{
//		//вычисляем шаг
//		h = setka[i + 1] - setka[i];
//		//lambda1 = lambda(setka[i]);
//		lambda1 = (lambda(setka[i]) + lambda(setka[i + 1])) / 2.;
//		//lambda2 = lambda(setka[i + 1]);
//		//вычисляем коэффициенты для разложения коэффициента диффузии
//
//		//вычисляем матрицу жескости
//
//		//di[3 * i] += (lambda1*(131. / 40.) + lambda2*(17. / 10.)) / h;
//		//di[3 * i + 1] += (lambda1*(999. / 20.) + lambda2*(675. / 4.)) / h;
//		//di[3 * i + 2] += (lambda1*(27. / 8.) + lambda2*(297. / 40.)) / h;
//		//di[3 * i + 3] += (lambda1*(17. / 40.) + lambda2*(131. / 40.)) / h;
//
//		//gg[ig[3 * i + 2] - 2] += (lambda1*(-109. / 2160.) + lambda2*(-17. / 2160.)) / h;
//		//gg[ig[3 * i + 2] - 1] += (lambda1*(39. / 40.) + lambda2*(3. / 8.)) / h;
//		//gg[ig[3 * i + 2]] += (lambda1*(1. / 4.) + lambda2*(-1. / 80.)) / h;
//		//gg[ig[3 * i + 3] - 1] += (lambda1*(-297. / 80.) + lambda2*(-297. / 80.)) / h;
//		//gg[ig[3 * i + 3]] += (lambda1*(3. / 80.) + lambda2*(39. / 40.)) / h;
//		//gg[ig[3 * i + 3] + 1] += (lambda1*(-153. / 80.) + lambda2*(-981. / 80.)) / h;
//
//
//		di[3 * i] += lambda1*(148.) / (40.*h);
//		di[3 * i + 1] += lambda1*(432.) / (40.*h);
//		di[3 * i + 2] += lambda1*(432.) / (40.*h);
//		di[3 * i + 3] += lambda1*(148.) / (40.*h);
//
//		gg[ig[3 * i + 1]] += lambda1*(-189.) / (40.*h);
//		gg[ig[3 * i + 2]] += lambda1*(54.) / (40.*h);
//		gg[ig[3 * i + 2] + 1] += lambda1*(-297.) / (40.*h);
//		gg[ig[3 * i + 3]] += lambda1*(-13.) / (40.*h);
//		gg[ig[3 * i + 3] + 1] += lambda1*(54.) / (40.*h);
//		gg[ig[3 * i + 3] + 2] += lambda1*(-189.) / (40.*h);
//	}
//
//
//	return;
//}
//
//
//void matr_M()		//добавляет в глобальную матрицу матрицу масс
//{
//	int i;
//	double h;
//
//	for (i = 0; i<count_elem; i++)
//	{
//		//вычисляем шаг
//		h = setka[i + 1] - setka[i];
//
//		//вычисляем матрицу массы
//
//
//		di[3 * i] += (h*gamma*128.) / 1680.0;
//		di[3 * i + 1] += (h*gamma*648.) / 1680.0;
//		di[3 * i + 2] += (h*gamma*648.) / 1680.0;
//		di[3 * i + 3] += (h*gamma*128.) / 1680.0;
//
//		gg[ig[3 * i + 1]] += (h*gamma*99.) / 1680.0;
//		gg[ig[3 * i + 2]] += (-h*gamma*36.) / 1680.0;
//		gg[ig[3 * i + 2] + 1] += (-h*gamma*81.) / 1680.0;
//		gg[ig[3 * i + 3]] += (h*gamma*19.) / 1680.0;
//		gg[ig[3 * i + 3] + 1] += (-h*gamma*36.) / 1680.0;
//		gg[ig[3 * i + 3] + 2] += (h*gamma*99.) / 1680.0;
//
//
//	}
//
//	return;
//}
//
//void F_vect()//генерирует вектор правой части, массив должен быть чистым
//{
//	int i;
//	double h, a, b, c, d;
//	double c11, c22, c33, c44, c21, c31, c32, c41, c42, c43;
//
//	for (i = 0; i<count_elem; i++)
//	{
//		//вычисляем шаг
//		h = setka[i + 1] - setka[i];
//
//		//вычисляем вектор правой части
//		a = func(setka[i]);
//		b = func(setka[i] + h / 3.);
//		c = func(setka[i] + 2 * h / 3.);
//		d = func(setka[i + 1]);
//
//		f[3 * i] += h*(128.*a + 99.*b - 36.*c + 19.*d) / 1680.;
//		f[3 * i + 1] += h*(99.*a + 648.*b - 81.*c - 36. * d) / 1680.;
//		f[3 * i + 2] += h*(-36.*a - 81.*b + 648.*c + 99. * d) / 1680.;
//		f[3 * i + 3] += h*(19.*a - 36.*b + 99.*c + 128. * d) / 1680.;
//
//	}
//	return;
//}
//
//
//
//
//
//
//
//
//
