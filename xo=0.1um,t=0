#include<iostream>
#include<complex>
#include <fstream>
using namespace std;

double Delta_b;
double Delta_r;
double hbar = 1.056 * pow(10, -34);
double const pi = 3.14159;
double const k = 1.38 * pow(10, -23);
double const T = 1 * pow(10, -6);
double m_at = 1.44316 * pow(10, -25);
complex<double> l;
complex<double> b, n, r;
void Define_E(double Int);
double E[3];
double matrix_bn, matrix_nr;
void Deltas_calc()
{
	Delta_b = 0;
	Delta_r = 0;
}

//int const size_mass = 10;
//complex<double> mass[3][size_mass];

double coor = pow(10, -7);
double p_step = 0.02 * hbar / coor;

double define_p_edge()
{
	complex<double> sum = 0;
	double px_edge = p_step*0.5;
	while (abs(sum) < (0.5 - pow(10, -8)))
	{
		sum += pow(pi * m_at / 0.5 * k * T, -0.5) * exp(-px_edge * px_edge * 0.5 / m_at / k / T) * p_step;
		px_edge += p_step;
	}
	cout << "p_edge equals" << " " << px_edge << endl;
	cout << "size is" << " " << 2 * px_edge / p_step << endl;
	return px_edge;
}
int size_mass;

void Init(complex<double>* cb_wf, double cb, double cr)
{
	double sum = 0;
	int a = 0, b = 0;
	double px;
	for (int i = 0; i < 2 * size_mass; i++)
	{
		px = (i - size_mass) * p_step;
		cb_wf[i] = cb * pow(pi * m_at / 0.5 * k * T, -0.5) * exp(-(px * px) * 0.5 / m_at / k / T);
		//cout << cb_wf[i]<<" "<<i << endl;
		sum += abs(cb_wf[i]) * p_step;// pow(pi * m_at / 0.5 * k * T, -0.5)* exp(-px * px * 0.5 / m_at / k / T)* (cb * cb + cr * cr)* p_step;
	}

	if (abs(sum - 1) > 0.01)
		cout << "problems with normalization" << endl;
}
ofstream fout("2.txt");
complex<double> a_;
complex<double> b_;

complex<double> exp_div_delta(double delta, double px_, double px, double time)
{
	double d = -delta - (px_ * px_ * 0.5 / m_at - px * px * 0.5 / m_at) / hbar;
	return exp(l * d * time) / d;
	//return double(1) / d;
}

complex<double> exp_delta(double delta, double px_, double px, double time)
{
	double d = -delta - (px_ * px_ * 0.5 / m_at - px * px * 0.5 / m_at) / hbar;
	return exp(l * d * time);
	//return double(1) / d;
}

/*void Runge_Cutte(double p,double time_step)
{
	double k1 = 2 * pi / (780 * pow(10, -9)), k2 = 2 * pi / (480 * pow(10, -9));
	double delta_p = 0, delta_1 = 1, delta_2 = 1, delta_3;
	complex<double> element1, element2;
	double time_of_op = 0;
	complex<double> Omega_eff, Omega_b, Omega_r;
	//complex<double> Omega_eff_conj;
	double Omega_nb = 14.1 * pow(10, 6), Omega_rn = 14.1 * pow(10, 6);
	complex<double> sum;
	double a[4];
	complex<double> mass_k[2][4];
	double px;
	double koef;
	a[0] = p;
	a[1] = p + p_step * 0.5;
	a[2] = p + p_step * 0.5;
	a[3] = p + p_step;
	complex<double> k_1, k_2;
	double wo = pow(10, -6);
	double x = 0.1 * wo;
	double delta = 50 * pow(10, 6);
	koef = double(1);// / double(6);
	for (int i = 0; i < 4; i++)
	{
		px = a[i];
		if (i == 0)
		{
			k_1 = 0;
			k_2 = 0;

		}
		else if (i == 3)
		{
			k_1 = time_step * mass_k[0][i - 1];
			k_2 = time_step * mass_k[1][i - 1];
		}
		else
		{
			k_1 = time_step * mass_k[0][i - 1] * 0.5;
			k_2 = time_step * mass_k[1][i - 1] * 0.5;
		}
		//delta_1 = 50 * pow(10, 6) + (hbar * hbar * k2 * k2 - 2 * (px + hbar * k1) * hbar * k2) * 0.5 / m_at / hbar;
		//delta_2 = 50 * pow(10, 6) + (hbar * hbar * k1 * k1 - 2 * (px + hbar * k1) * hbar * k1) * 0.5 / m_at / hbar;
		//delta_3 = 50 * pow(10, 6) + (hbar * hbar * k2 * k2 - 2 * (px + hbar * k2) * hbar * k2) * 0.5 / m_at / hbar;
		//delta_p = -1 / hbar * pow((hbar * k2 - hbar * k1), 2) * 0.5 / m_at - px * (k2 - k1) / m_at;
		//Omega_eff = double(Omega_1 * Omega_2) / delta_1 * l;
		//Omega_b = double(Omega_1 * Omega_1) / delta_2 * l;
		//Omega_r = double(Omega_2 * Omega_2) / delta_3 * l;
		//cout << delta_p << " " << delta_3 << endl;
		//element1 = pow(double(2) * l , -2) * Omega_eff* exp(l * delta_p * time);
		//element2 = pow(double(2) * l, -2) * mass[0][int(px / p_step + size_mass * 0.5)] * Omega * exp(l * delta_3 * time) * hbar * l;
		//if (abs(element2)> pow(10,-10))
			//cout << 1 << endl;
		//cout << int(px / p_step + size_mass * 0.5) << endl;
		//mass[0][int(px / p_step + size_mass * 0.5)];
		mass_k[0][i] = 0;
		mass_k[1][i] = 0;
		mass_k[0][i]= l * 0.5 * (c_n[int((a[i] + hbar * k1) / p_step + size_mass * 0.5)]+k_1) * Omega_nb * exp_delta(delta, a[i] + hbar * k1, a[i], time_step);
		mass_k[1][i] = l * 0.5 * (c_n[int((a[i] + hbar * k2) / p_step + size_mass * 0.5)] + k_2) * Omega_nb * exp_delta(delta, a[i] + hbar * k2, a[i], time_step);
		//cout << koef * pow(double(2) * l, -2) * Omega_eff * exp(l * delta_p * time_step) << " " << koef  * Omega_b << endl;
		//mass_k[0][i] = koef * pow(double(2) * l, -2) * Omega_eff * exp(l * delta_p * time_step) * (mass[1][int((px + hbar * (k1 - k2)) / p_step + size_mass * 0.5)] + k_1);
		//mass_k[0][i] -= koef * (mass[0][int(px / p_step + size_mass * 0.5)] + k_2) * Omega_b;
		//sum += mass[0][int(px / p_step + size_mass * 0.5)]*p_step;
		//mass_k[1][i] = koef * pow(double(2) * l, -2) * (-Omega_eff) * exp(-l * delta_p * time_step) * (mass[0][int(px / p_step + size_mass * 0.5)] + k_2);//  +
		//mass_k[1][i] += koef * (mass[1][int((px + hbar * (k1 - k2)) / p_step + size_mass * 0.5)] + k_1) * Omega_r;
		//sum+= pow(double(2) * l, -2) * Omega_eff * exp(l * delta_p * time) * mass[1][int((px + hbar * (k1 - k2)) / p_step + size_mass * 0.5)] * time_step +
			//mass[0][int(px / p_step + size_mass * 0.5)] * Omega_b * l * hbar*time_step;
		//mass[1][int(px / p_step) + 5000] += element1 * time_step + element2 * time_step;
	}
	a_ = (mass_k[0][0] + double(2) * mass_k[0][1] + double(2) * mass_k[0][2] + mass_k[0][3]) / double(6);
	b_ = (mass_k[1][0] + double(2) * mass_k[1][1] + double(2) * mass_k[1][2] + mass_k[1][3]) / double(6);
}
*/

/*complex<double>** help_massive_11;// [size_help_mass] [size_help_mass] ;
complex<double>** help_massive_12;// [size_help_mass] [size_help_mass] ;
complex<double>** help_massive_21;// [size_help_mass] [size_help_mass] ;
complex<double>** help_massive_22;// [size_help_mass] [size_help_mass] ;

void help_massive_count1()
{
	double k_;
	double wo = pow(10, -6);
	double x = 0.1 * wo;
	double Omega_nb = 14.1 * pow(10, 6), Omega_rn = 14.1 * pow(10, 6);
	double k1 = 2 * pi / (780 * pow(10, -9)), k2 = 2 * pi / (480 * pow(10, -9));


	help_massive_11 = new complex<double> * [size_help_mass];
	help_massive_12 = new complex<double> * [size_help_mass];
	help_massive_21 = new complex<double> * [size_help_mass];
	help_massive_22 = new complex<double> * [size_help_mass];
	for (int i = 0; i < size_help_mass; i++)
	{
		help_massive_11[i] = new complex<double>[size_help_mass];
		help_massive_12[i] = new complex<double>[size_help_mass];
		help_massive_21[i] = new complex<double>[size_help_mass];
		help_massive_22[i] = new complex<double>[size_help_mass];
	}

	for (int i = 0; i < size_help_mass; i++)
		for (int j = 0; j < size_help_mass; j++)
		{
			help_massive_11[i][j] = 0;
			help_massive_12[i][j] = 0;
			help_massive_21[i][j] = 0;
			help_massive_21[i][j] = 0;
		}
	for (double px_ = -px_edge; px_ < px_edge; px_ += p_step)
		for (double px = -px_edge; px < px_edge; px += p_step)
		{
			k_ = (-px_ / hbar - k1 + px / hbar);
			help_massive_11[int(px / p_step + size_help_mass * 0.5)][int(px_ / p_step + size_help_mass * 0.5)] = Omega_nb * (exp(l * k_ * x) * (-l * k_ * k_ * x * x + 2 * k_ * x + double(2) * l) - exp(-l * k_ * x) * (-l * k_ * k_ * x * x - 2 * k_ * x + double(2) * l)) / k_ / k_ / k_;
			k_ = (-px_ / hbar - k2 + px / hbar);
			help_massive_12[int(px / p_step + size_help_mass * 0.5)][int(px_ / p_step + size_help_mass * 0.5)] = Omega_rn * (exp(l * k_ * x) * (-l * k_ * k_ * x * x + 2 * k_ * x + double(2) * l) - exp(-l * k_ * x) * (-l * k_ * k_ * x * x - 2 * k_ * x + double(2) * l)) / k_ / k_ / k_;
			k_ = (-px_ / hbar + k1 + px / hbar);
			help_massive_21[int(px / p_step + size_help_mass * 0.5)][int(px_ / p_step + size_help_mass * 0.5)] = Omega_nb * (exp(l * k_ * x) * (-l * k_ * k_ * x * x + 2 * k_ * x + double(2) * l) - exp(-l * k_ * x) * (-l * k_ * k_ * x * x - 2 * k_ * x + double(2) * l)) / k_ / k_ / k_;
			k_ = (-px_ / hbar + k2 + px / hbar);
			help_massive_22[int(px / p_step + size_help_mass * 0.5)][int(px_ / p_step + size_help_mass * 0.5)] = Omega_rn * (exp(l * k_ * x) * (-l * k_ * k_ * x * x + 2 * k_ * x + double(2) * l) - exp(-l * k_ * x) * (-l * k_ * k_ * x * x - 2 * k_ * x + double(2) * l)) / k_ / k_ / k_;
		}
}*/

void Euler(complex<double>* cb_wf, complex<double>* cr_wf, int j, double time)
{
	double k1 = 2 * pi / (780 * pow(10, -9)), k2 = 2 * pi / (480 * pow(10, -9));
	double q = k1 - k2;
	double delta_p = 0, delta_1 = 1, delta_2 = 1, delta_3;
	complex<double> element1, element2;
	double time_of_op = 0;
	complex<double> Omega_eff, Omega_b, Omega_r;
	//complex<double> Omega_eff_conj;
	double Omega_nb = 14.1 * pow(10, 6), Omega_rn = 14.1 * pow(10, 6);
	complex<double> sum;
	complex<double> Int1 = 0, Int2 = 0;
	a_ = 0;
	b_ = 0;
	double wo = pow(10, -6);
	double x1 = 2 * wo;
	double x2 = 2 * wo;
	double x_gen = 2 * x1 * x2 / (x1 + x2);
	double delta = 50 * pow(10, 6);
	double ep1, ep2, ep;
	int i = j - hbar * q / p_step;
	double px_ = (i - size_mass) * p_step;
	ep1 = pow(px_ + hbar * q, 2) / 2 / m_at;
	ep2 = pow(px_ - hbar * q, 2) / 2 / m_at;
	ep = pow(px_, 2) / 2 / m_at;
	a_ = 0;
	b_ = 0;
	//a_ += l * 0.5 * c_n[int((px_ + hbar * k1) / p_step + size_mass * 0.5)] * Omega_nb * exp_delta(delta, px_ + hbar * k1, px_, time);
	//b_ += l * 0.5 * c_n[int((px_ + hbar * k2) / p_step + size_mass * 0.5)] * Omega_rn * exp_delta(delta, px_ + hbar * k2, px_, time);
	a_ = -l / double(4) / delta * pow(Omega_nb * hbar, 2) * (double(0) * cb_wf[i] / hbar / hbar + double(2) / p_step / p_step * (cb_wf[i - 1] - double(2) * cb_wf[i] + cb_wf[i + 1]) / x1 / x1);
	//a_ += -l / double(4) / delta * Omega_nb * hbar * Omega_rn * hbar * (double(0) * cr_wf[j] / hbar / hbar - double(0) / p_step / p_step * (cr_wf[j - 1] - double(2) * cr_wf[j] + cr_wf[j + 1]) / x_gen / x_gen) * exp(-l / hbar * (ep1 - ep) * double(0));
	if
		(a_ != a_)
	{
		cout << 1 << endl;
		cout << cr_wf[j - 1] - double(2) * cr_wf[j] + cr_wf[j + 1] << endl;
	}
	//b_ = -l / double(4) / delta * pow(Omega_nb * hbar, 2) * (double(0) * cr_wf[j] / hbar / hbar + double(0) / p_step / p_step * (cr_wf[j - 1] - double(2) * cr_wf[j + 1] + cr_wf[j + 1]) / x1 / x1);
	//b_ += -l / double(4) / delta * Omega_nb * hbar * Omega_rn * hbar * (double(0) * cb_wf[i] / hbar / hbar + double(0) / p_step / p_step * (cb_wf[i - 1] - double(0) * cb_wf[i] + cb_wf[i + 1]) / x_gen / x_gen) * exp(-l / hbar * (ep2 - ep) * double(0));
	if (j == 2)
	{
		cout << (cb_wf[i - 1] - double(2) * cb_wf[i] + cb_wf[i + 1]) << endl;
		cout << cb_wf[i - 1] << " " << cb_wf[i] << " " << cb_wf[i + 1] << endl;
	}
	//cout << -l / double(4) / delta * Omega_nb * hbar * Omega_rn * hbar * cb_wf[i] / hbar / hbar << " " << -l / double(4) / delta * Omega_nb * hbar * Omega_rn * hbar * ( double(2) / p_step / p_step * (cb_wf[i - 1] - double(2) * cb_wf[i] + cb_wf[i + 1]) / x_gen / x_gen)<< endl;
	if
		(b_ != b_)
	{
		cout << 2 << endl;
		cout << cb_wf[i - 1] << " " << cb_wf[i] << " " << cb_wf[i + 1] << endl;
		cout << (double(2) / p_step / p_step * (cb_wf[i - 1] - double(2) * cb_wf[i] + cb_wf[i + 1]) / x_gen / x_gen) * exp(-l / hbar * (ep2 - ep) * time) << endl;
	}
}


void Define_E(double Int)
{
	E[1] = 0.482174 * pow(Int / 3.8 * 1000, 0.5) * pow(10, 6);
	E[2] = 0.482174 * pow(Int / 3.8 * 1000, 0.5) * pow(10, 6);
}
complex<double> cb[9306];
complex<double> cr[9306];

void Rabi_osc(complex<double>* cb_wf, complex<double>* cr_wf, double time)
{
	double k1 = 2 * pi / (780 * pow(10, -9)), k2 = 2 * pi / (480 * pow(10, -9));
	double q = k1 - k2;
	double time_step = pow(10, -8);
	double time_of_op = 0;
	complex<double> Omega_eff, Omega_b, Omega_r;
	//complex<double> Omega_eff_conj;
	double Omega_1 = 50 * pow(10, 6), Omega_2 = 50 * pow(10, 6);
	complex<double> sum;
	double a, b;
	complex<double> mass_k[2][4];
	double px;
	double koef;
	int i;
	complex<double> b__[2] = { 0,0 }, r__[2] = { 0,0 };
	while (abs(time - time_of_op) > 0.5 * time_step)
	{
		sum = 0;
		for (int i = 0; i < 2 * size_mass; i++)
		{
			cb[i] = 0;
			cr[i] = 0;
		}
		for (int j = 1; j < 2 * size_mass - 1; j++)
		{
			i = j - hbar * q / p_step;
			if (abs(i - size_mass) > size_mass - 1) continue;
			Euler(cb_wf, cr_wf, j, time_step);
			cb[i] = -a_ * time_step;
			//cout << abs(a_) << endl;
			//cr[j] = -b_ * time_step;
		}
		for (int j = 1; j < 2 * size_mass - 1; j++)
		{
			i = j - hbar * q / p_step;
			if (abs(i - size_mass) > size_mass - 1) continue;
			cb_wf[i] += cb[i];
			//fout << abs(cb[j]) << " " << j << endl;
			//cr_wf[j] += cr[j];
		}

		//cout << sum << endl;
		time_of_op += time_step;
	}
}
void proverka()
{
	double time_step = 0.0001;
	complex<double> a = 1, b = 0;
	for (double time = 0; time < 1; time += time_step)
	{
		a -= l * double(10) * b * time_step;
		b -= l * double(10) * a * time_step;
		cout << pow(abs(a), 2) << endl;
	}
	cout << 1;
}

void trace_out_of_impulse(complex<double>* cb_wf, complex<double>* cr_wf, bool fl)
{
	b = r = 0;
	for (int i = 1; i < 2 * size_mass - 1; i++)
	{
		b += cb_wf[i] * p_step;
		//cout << b <<" "<<i<< endl;
		r += cr_wf[i] * p_step;
		//cout << abs(cb_wf[i])<<" "<<abs(cr_wf[i])<<" "<<abs(b)<<" "<<abs(r) <<" "<< i << endl;

		if (fl == 1)
			fout << abs(cr_wf[i - 1] - double(2) * cr_wf[i] + cr_wf[i + 1]) << " " << i << endl;
		fout << abs(cr_wf[i])  << endl;

	}
}

int main()
{
	//imaginary unit 
	l = -1;
	l = sqrt(l);
	//proverka();
	double k1 = 2 * pi / (780 * pow(10, -9)), k2 = 2 * pi / (480 * pow(10, -9));
	double q = k1 - k2;
	//Defining of p_edge
	static const double p_edge = define_p_edge();

	//Creating arrays of wave_functions
	size_mass = p_edge / p_step;
	complex<double>* c_b = new complex<double>[int(2 * size_mass)];
	//complex<double>* c = new complex<double>[int(2 * size_mass)];
	Init(c_b, 1, 0);
	double x1 = pow(10, -6);
	double time_step = pow(10, -10);
	double Omega_nb = 14.1 * pow(10, 6), Omega_rn = 14.1 * pow(10, 6);
	double delta = 50 * pow(10, 6);
	double ep, ep1, ep2;
	double px_;
	double w=60000;
	complex<double> c[10000];
	double po = pow(hbar*m_at*w,0.5);
	double xo = hbar / po;
	//1
	fout.precision(15);
	for (double time = 100; time < 10; time += time_step)
	{
		for (int i = 0; i < 2 * size_mass; i++)
		{
			c[i] = 0;
		}
		for (double p_x = -3 * po; p_x < 3 * po; p_x += p_step)
		{
			complex<double>sum = 0;
			ep1 = p_x * p_x / 2 / m_at;
			for (double p__x = -3 * po; p__x < 3 * po; p__x += p_step)
			{
				ep2 = p__x * p__x / 2 / m_at;
				complex<double>sum1 = 0;
				for (double x = -x1; x < x1; x += x1 * 0.0005)
				{
					sum1 += exp(l / hbar * (p__x - p_x)*x) * exp(-x * x / x1 / x1)* x1 * 0.0005;
				}
				sum += l*sum1*exp(-p__x* p__x/po/po)*p_step/pi/hbar/double(8)*exp(-l*time/hbar*(ep2-ep1));
			}
			fout << abs(sum) << " " << p_x / po << endl;

		}
		cout << "heh" << endl;
	}
	for (double time = pow(10, -4); time < 1000*pow(10, -4); time += time_step)
	{
		for (int i = 0; i < 2 * size_mass; i++)
		{
			c[i] = 0;
		}
		for (double p_x = -3 * po; p_x < 3 * po; p_x += p_step)
		{
			complex<double>sum = l * exp(-p_x * p_x / po / po) / double(4);
			//sum += hbar*hbar/x1/x1*l * pow(Omega_nb, 2) / delta * exp(-p_x * p_x / po / po) * ((-2*p_x/po/po)*(double(1)+l * p_x / m_at / hbar * time) + (l*time/m_at/hbar-pow(p_x * time / m_at / hbar, 2))) / double(4);
			sum += hbar * hbar / x1 / x1 * l   * exp(-p_x * p_x / po / po) * ((-2 / po / po + l * time / m_at / hbar) + pow(-2 * p_x / po / po + l * p_x * time / m_at  / hbar, 2)) / double(4);
			fout << abs(sum) << " " << p_x /po << endl;
		}
		cout << "heh" << endl;
	}
	for (int i = 0; i < 2 * size_mass; i++)
	{
		fout << abs(c_b[i]) << " " << i << endl;
		cout << c[i] << endl;
	}
	//fout << pow(abs(b), 2) << " " << pow(abs(n), 2) << " " << pow(abs(r), 2) << endl;
	fout.close();
	system("pause");
	return 0;
}
