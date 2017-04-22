#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;
const double eps = 1e-12; // sta³a przybli¿enia zera

struct Element {
	int id1;
	int id2;
	double S;
	double L;
	double k;
	double localH[2][2];
	double localP[2];
};

struct Wezel {
	double temp;
	double warunek; //warunki brzegowe
	double x; //polozenie
};

struct Siatka {
	int iw;	//ilosc wezlow
	int ie; //ilosc elementow
	double q;
	double alfa;
	double temp_ot;
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

bool gauss(int n, double ** AB, double * X);

int main(){

    Siatka siatka;
	Rozw_uk_rownan rozw;
	siatka.GenerujSiatkeZPliku();
	siatka.DrukujSiatke();
	//cout << siatka.dane->id1;

	////////////////////////////////////////////////

	for (int i = 0; i < siatka.ie; i++){
		double c = siatka.dane[i].S*siatka.dane[i].k / siatka.dane[i].L;
		for (int j = 0; j < 2; j++){
			for (int k = 0; k < 2; k++){
				if ((j + k) % 2 == 0){
					if (i == siatka.ie - 1 && j + k == 2) siatka.dane[i].localH[j][k] = c + siatka.alfa*siatka.dane[i].S;
					else	siatka.dane[i].localH[j][k] = c;
				}
				else siatka.dane[i].localH[j][k] = -c;
			}
		}
		siatka.dane[i].localP[0] = siatka.dane2[i].warunek*-1;
		siatka.dane[i].localP[1] = siatka.dane2[i + 1].warunek*-1;
	}

	//Wyswietlanie
	/*
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			cout << siatka.dane[i].localH[i][j] <<" ";
		}
		cout << endl;
	}
	for (int i = 0; i < 2; i++){
		cout << siatka.dane[i].localP[i] <<endl;
	}
	*/
/////////////////////////////////////////////
	rozw.GH = new double*[siatka.iw];
	rozw.GP = new double[siatka.iw];
	rozw.GT = new double[siatka.iw];
	
		for (int i = 0; i < siatka.iw; i++){
			rozw.GH[i] = new double[siatka.iw];
		}
	
		for (int i = 0; i < siatka.iw; i++){
			rozw.GP[i] = 0; rozw.GT[i] = 0;
			for (int j = 0; j < siatka.iw; j++){
			rozw.GH[i][j] = 0;
			}
		}

	for (int i = 0; i < siatka.ie; i++){
		for (int j = 0; j < 2; j++){
			rozw.GP[i+j] += siatka.dane[i].localP[j];
			for (int k = 0; k < 2; k++){
				rozw.GH[i+j][i+k] += siatka.dane[i].localH[j][k];
			}
		}
	}
	//Wyœwietlanie
	/*
	for (int i = 0; i < siatka.ie; i++){
		cout << rozw.GP[i] << endl;
		}
	for (int i = 0; i < siatka.iw; i++){
		for (int j = 0; j < siatka.iw; j++){
			cout << rozw.GH[i][j] << " ";
		}
		cout << endl;
	}
	*/
	///////////////////////////////////////////////////////////////////////////////
	cout << "METODA ELIMINACJI GAUSSA\n"<<endl;
	
	//scalanie
	double **AB;
	AB = new double*[siatka.iw];
	for (int i = 0; i < siatka.iw; i++){
		AB[i] = new double[siatka.iw + 1];
	}
	for (int i = 0; i < siatka.iw; i++){
		for (int j = 0; j < siatka.iw; j++) {
			AB[i][j] = rozw.GH[i][j];
		}
		AB[i][siatka.iw] = rozw.GP[i];
	}

	//Wyœwietlanie
	cout << "Macierz po scaleniu: " << endl;
	for (int i = 0; i < siatka.iw; i++){
		for (int j = 0; j <= siatka.iw; j++) {
			cout << AB[i][j]<<" ";
		}
		cout << endl;
	}
	cout << endl;
	
	double *X;
	X = new double[siatka.iw];
	if(gauss(siatka.iw, AB,rozw.GT))
	{
		for (int i = 0; i < siatka.iw; i++)
			cout << "t"<<i  <<" = " << rozw.GT[i] << endl;
			
	}
	
	
	system("PAUSE");
	return 0;
}


void Siatka::GenerujSiatkeZPliku(){

	fstream plik;
	plik.open("dane.txt", std::ios::in);
	if (plik.good())
	{
		plik >> this->iw;
		this->ie = this->iw - 1;
		plik >> this->q;
		plik >> this->alfa;
		plik >> this->temp_ot;
		this->dane = new Element[ie];
		this->dane2 = new Wezel[iw];

		for (int i = 0; i < this->iw; i++){
			plik >> this->dane2[i].warunek;
			plik >> this->dane2[i].x;
		}

		for (int i = 0; i < this->ie; i++){
			plik >> this->dane[i].id1;
			plik >> this->dane[i].id2;
			plik >> this->dane[i].S;
			plik >> this->dane[i].k;
		}

		if (this->dane2[0].warunek == 1){
			for (int i = 0; i < this->iw; i++)
			{
				if (this->dane2[i].warunek == 1)
				{
					this->dane2[i].warunek = this->q*this->dane[i].S;
				}
				else if (this->dane2[i].warunek == 2)
				{
					this->dane2[i].warunek = -(this->temp_ot * this->dane[i - 1].S)*this->alfa;
				}
				else this->dane2[i].warunek = 0;
			}
		}
		
		for (int i = 0; i < this->ie; i++){
			this->dane[i].L = this->dane2[i+1].x - this->dane2[i].x;
		}
		plik.close();
	}
	else cout << "ERROR ";
}
void Siatka::DrukujSiatke(){
	cout << "Liczba wezlow:" << this->iw << endl;
	cout << "Liczba elementow:" << this->ie << endl;
	cout << "Strumien ciepla:" << this->q << endl;
	cout << "Alfa:" << this->alfa << endl;
	cout << "Temepratura otoczenia:" << this->temp_ot << endl;
	for (int i = 0; i < this->ie; i++){
		cout << "Nr el: " << i << "id1: " << this->dane[i].id1 << "id2: " << this->dane[i].id2 << "Powierzchia el:" << this->dane[i].S << "Wspolczynnik k:" << this->dane[i].k << endl;
	}
	for (int i = 0; i < this->iw; i++){
		cout << "Nr wezla: " << i << "Warunek brzegowy: " << this->dane2[i].warunek << endl;;
	}
	for (int i = 0; i < this->ie; i++){
		cout << "Dlugosc wezla: " << i << "wynosi: " << this->dane[i].L << endl;
	}
}

bool gauss(int n, double ** AB, double * X)
{

	int i, j, k;
	double m, s;

	// eliminacja wspó³czynników

	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (fabs(AB[i][i]) < eps) return false;
			m = -AB[j][i] / AB[i][i];
			for (k = i + 1; k <= n; k++)
				AB[j][k] += m * AB[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = AB[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= AB[i][j] * X[j];
		if (fabs(AB[i][i]) < eps) return false;
		X[i] = s / AB[i][i];
	}
	return true;
}