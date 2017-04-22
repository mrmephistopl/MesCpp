#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;
const double eps = 1e-12; // sta³a przybli¿enia zera

struct Wezel {
	double t;
	double x;
};

struct Element {
	int id1;
	int id2;
	double l;
	double t[2];
	double localH[2][2];
	double localP[2];
};

struct s {
	int iw;	//ilosc wezlow
	int ie; //ilosc elementow
	double rMax, rMin;
	double time;
	double dTau, alfa;
	double tempBegin, tempEnv;
	double C, Ro, k, dR; 
	double a;
	double W[2], E[2], N1[2], N2[2];
	Element *dane;
	Wezel *dane2;
	void GenerujSiatkeZPliku();
	void DrukujSiatke();
};


struct Rozw_uk_rownan{
	double **GH;
	double *GP;
	double *GT;
};

bool gauss(s s, int n, double ** AB, double * X);


int main(){

	fstream plik1;
	s s;

	int np = 2;
	double tpTau, Rp;
	s.GenerujSiatkeZPliku();
	s.DrukujSiatke();
	double t0 = s.tempBegin;
	double t1 = s.tempBegin;

	for (int x = 0; x < s.time; x += s.dTau){
		for (int i = 0; i < s.ie; i++){

			if (x == 0)
			{
				s.dane[i].t[0] = t0;
				s.dane[i].t[1] = t1;
			}
			double x1 = s.dane2[i].x;
			double x2 = s.dane2[i + 1].x;
			//Zerowanie macierzy lokalnych

			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					s.dane[i].localH[j][k] = 0;
				}
				s.dane[i].localP[j] = 0;
			}

			for (int j = 0; j < np; j++) {
				tpTau = s.N1[j] * s.dane[i].t[0] + s.N2[j] * s.dane[i].t[1];
				Rp = (s.N1[j] * x1) + (s.N2[j] * x2);

				s.dane[i].localH[0][0] += s.k*Rp * 1 / s.dR + s.C*s.Ro*s.dR*Rp * 1 * s.N1[j] * s.N1[j] / s.dTau;
				s.dane[i].localH[0][1] += -s.k*Rp * 1 / s.dR + s.C*s.Ro*s.dR*Rp * 1 * s.N1[j] * s.N2[j] / s.dTau;
				s.dane[i].localH[1][0] = s.dane[i].localH[0][1];
				s.dane[i].localH[1][1] += s.k*Rp * 1 / s.dR + s.C*s.Ro*s.dR*Rp * 1 * s.N2[j] * s.N2[j] / s.dTau;
				s.dane[i].localP[0] += -s.C*s.Ro*s.dR / s.dTau * tpTau*s.N1[j] * Rp * 1;
				s.dane[i].localP[1] += -s.C*s.Ro*s.dR / s.dTau * tpTau*s.N2[j] * Rp * 1;

			}
			if (i == s.ie - 1)
			{
				s.dane[i].localH[1][1] += (2 * s.rMax * s.alfa);
				s.dane[i].localP[1] += -2 * s.alfa*s.rMax*s.tempEnv;
			}
		}

		//Wyswietlanie macierzy lokalnych


		for (int i = 0; i < s.ie; i++)
		{
			cout << "///////ELEMENT - " << i << "///////" << endl;
			cout << "id1 - " << s.dane[i].id1 << "\t";
			cout << "id2 - " << s.dane[i].id2 << "\t";
			cout << "l - " << s.dane[i].l << "\t";
			cout << "LH" << endl;
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					cout << s.dane[i].localH[j][k] << "\t";
				}
				cout << endl;
			}
			cout << endl;
			cout << "LP" << endl;
			for (int j = 0; j < 2; j++)
			{
				cout << s.dane[i].localP[j] << "\t";
			}
			cout << endl;
		}


	//
		Rozw_uk_rownan rozw;
		rozw.GH = new double*[s.iw];
		rozw.GP = new double[s.iw];
		rozw.GT = new double[s.iw];
		{
			for (int j = 0; j < s.iw; j++)
			{
				rozw.GH[j] = new double[s.iw];
			}
		}
		for (int i = 0; i < s.iw; i++)
		{
			rozw.GP[i] = 0;
			rozw.GT[i] = 0;
			for (int j = 0; j < s.iw; j++)
			{
				rozw.GH[i][j] = 0;
			}
		}
		//Uzupelnienie macierzy globalnych K
		for (int i = 0; i < s.ie; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				rozw.GP[j + i] += s.dane[i].localP[j];
				for (int k = 0; k < 2; k++)
				{
					rozw.GH[j + i][k + i] += s.dane[i].localH[j][k];
				}
			}
		}

		//Tworzenie macierzy F
		double ** HP = new double*[s.iw];
		double ** HP1 = new double*[s.iw];

		for (int i = 0; i < s.iw; i++)
		{
			HP[i] = new double[s.ie + 1];
			HP1[i] = new double[s.ie + 1];
		}
		for (int i = 0; i < s.iw; i++)
		{
			for (int j = 0; j < s.iw; j++)
			{
				HP[i][j] = rozw.GH[i][j];
				HP1[i][j] = rozw.GH[i][j];
			}
			HP[i][s.iw] = rozw.GP[i] * (-1);
			HP1[i][s.iw] = rozw.GP[i] * (-1);
		}
		//Wyswietlenie K
		cout << endl << "Macierz globalna K" << endl;
		for (int j = 0; j < s.iw; j++)
		{
			for (int k = 0; k < s.iw + 1; k++)
			{
				cout << HP[j][k] << "\t";
			}
			cout << endl;
		} //Wyœwietplenie F
		cout << endl << "Wektor F" << endl;
		for (int j = 0; j < s.iw; j++)
		{
			cout << rozw.GP[j] << "\t";

		}
		cout << endl;

		//Oblicznie Gaussa

		if (gauss(s, s.iw, HP, rozw.GT))
		{
			for (int i = 0; i < s.iw; i++)
				cout << "t" << i << " = " << setw(9) << s.dane2[i].t << endl;
		}
		else
		{
			cout << "DZIELNIK ZERO\n";
		}

		HP = HP1;

		//Zapis do pliku
		plik1.open("wyniki.txt", ios::app);
		if (plik1.is_open())
		{
			plik1 << endl << "////K/////" << endl;
			for (int j = 0; j < s.iw; j++)
			{
				for (int k = 0; k < s.iw + 1; k++)
				{
					plik1 << HP[j][k] <<  "\t";
				}
				plik1 << endl;
			}
			plik1 << endl << "////F/////" << endl;
			for (int j = 0; j < s.iw; j++)
			{
				plik1 <<  rozw.GP[j] <<  "\t";

			}plik1 << endl;
				plik1 << endl << "Temperatura po " << s.dTau + x << " sekundach: " << endl<< endl;
			for (int i = 0; i < s.iw; i++)
				plik1 << "t" << i << " = " << setw(9) << s.dane2[i].t << endl;
		}

		plik1.close();
	}



	cout << endl;
	system("PAUSE");
	return 0;
}

void s::GenerujSiatkeZPliku(){

	fstream plik;
	plik.open("dane.txt", std::ios::in);
	if (plik.good())
	{
		plik >> this->iw;
		this->ie = this->iw - 1;
		plik >> this->rMax;
		plik >> this->rMin;
		plik >> this->time;
		plik >> this->dTau;
		plik >> this->alfa;
		plik >> this->tempBegin;
		plik >> this->tempEnv;
		plik >> this->C;
		plik >> this->Ro;
		plik >> this->k;
		this->dR = (this->rMax - this->rMin) / this->ie;
		this->W[0] = this->W[1] = 1; // wektor wag dla punktow calkowania
		this->E[0] = -0.5773502692; // wspolrzedne punktow calkowania
		this->E[1] = -this->E[0];
		this->N1[0] = 0.5 * (1 - this->E[0]); // funkcje ksztaltu dla naszego przypadku
		this->N1[1] = 0.5 * (1 - this->E[1]);
		this->N2[0] = 0.5 * (1 + this->E[0]);
		this->N2[1] = 0.5 * (1 + this->E[1]);
		this->dane = new Element[ie];
		this->dane2 = new Wezel[iw];
		for (int i = 0; i < this->iw; i++){
			this->dane2[i].x = i*this->dR;
			this->dane2[i].t = this->tempBegin;
		}

		for (int i = 0; i < this->ie; i++){
			this->dane[i].id1 = i;
			this->dane[i].id2 = i + 1;
			this->dane[i].l = this->dane2[i + 1].x - this->dane2[i].x;
			for (int j = 0; j < 2; j++)
			{
				for (int k = 0; k < 2; k++)
				{
					this->dane[i].localH[j][k] = 0;
				}
				this->dane[i].localP[j] = 0;
			}
		}
	}
	else cout << "ERROR ";
}

void s::DrukujSiatke(){
	cout << "Liczba wezlow:" << this->iw << endl;
	cout << "Liczba elementow:" << this->ie << endl;
	cout << "Promien maksymalny:" << this->rMax << endl;
	cout << "Promien minimalny:" << this->rMin << endl;
	cout << "Skok promienia:" << this->dR << endl;
	cout << "Czas procesu:" << this->time << endl;
	cout << "Skok czasowy:" << this->dTau << endl;
	cout << "Temperatura poczatkowa:" << this->tempBegin << endl;
	cout << "Temperatura otoczenia:" << this->tempEnv << endl;
	cout << "Wspolczynnik wymiany ciepla:" << this->alfa << endl;
	cout << "Wspolczcynnik przewodzenia ciepla:" << this->k << endl;
	cout << "Gestosc:" << this->Ro << endl;
	cout << "Efektywne cieplo wlasciwe:" << this->C << endl;
}

bool gauss(s s, int n, double ** AB, double * X)
{

	double m, a;
	// eliminacja wspó³czynników
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
		{
			if (abs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (int k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}
	// wyliczanie niewiadomych

	for (int i = n - 1; i >= 0; i--)
	{
		a = AB[i][n];
		for (int j = n - 1; j >= i + 1; j--)
		{
			a -= AB[i][j] * s.dane2[j].t;
		}
		if (abs(AB[i][i]) < eps) return false;

		s.dane2[i].t = a / AB[i][i];
	}
	for (int i = 0; i < s.ie; i++) {
		s.dane[i].t[0] = s.dane2[i].t;
		s.dane[i].t[1] = s.dane2[i + 1].t;
	}
	return true;
}