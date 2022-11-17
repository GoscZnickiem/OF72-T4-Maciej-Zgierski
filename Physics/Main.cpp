// Kod programu stanowi�cy cz�� rozwi�zania zadania T4 I stopnia 72. Olimpiady Fizycznej
// Autor programu: Maciej Zgierski
// Listopad 2022

// Program wyznacza i zapisuje wyniki kilku symulacji do plik�w o nazwie T4_s<nr>.txt oraz T4_wsp_tarcia.txt
// Wypisuje tak�e u�yteczne informacje o symulacjach w konsoli

// W celu wykonania wykres�w nale�y pos�u�y� si� oddzielnym programem (np. arkuszem kalkulacyjnym)

#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

//Sta�e

const double g = 9.81;
const double L = 0.1;

const double CZAS_SYMULACJI = 10;

//Pocz�tkowe warto�ci symulacji

const double ALPHA = -M_PI / 2;
const double ALPHA_DOT = 0.0;

//Funkcja obliczaj�ca przyspieszenie k�towe wywo�ane momentem si�y zsuwaj�cej dzia�aj�ce na cia�o w danej chwili
double moment_sily_zsuwajacej(double alpha) {
    return - g / L * sin(alpha);
}

//Funkcja obliczaj�ca przyspieszenie k�towe wywo�ane momentem si�y tarcia dzia�aj�cej na poruszaj�cy si� klocek w danej chwili
double moment_sily_tarcia(double alpha, double omega, double f) {
    if (omega == 0) //kiedy klocek si� nie porusza tarcie kinetyczne nie jest obecne
        return 0;
    if (omega > 0) //sprawdzamy w kt�r� stron� porusza si� klocek w celu wyznaczenia kierunku dzia�ania tarcia
        return - f * (omega * omega + g / L * cos(alpha));
    return f * (omega * omega + g / L * cos(alpha));
}

//R�wnanie r�niczkowe dla alfy (f_alfa - patrz r�wnanie (7))
double f_alpha(double alpha, double omega) {
    return omega;
}

//R�wnanie r�niczkowe dla omegi (f_omega - patrz r�wnanie (7))
double f_omega(double alpha, double omega, double f) {
    return moment_sily_zsuwajacej(alpha) + moment_sily_tarcia(alpha, omega, f);
}

void symulacja(int numer_symulacji, double delta_t, double f, const char* nazwa_pliku) {

    //Przygotowanie symulacji
    double alpha = ALPHA;
    double omega = ALPHA_DOT;
    double czas = 0;

    double czas_zatrzymania = -1;
    bool stop = false; //stop jest identyfikatorem zatrzymania si� klocka

    //Otwarcie pliku tekstowego

    std::fstream plik;
    plik.open(nazwa_pliku, std::ios::out);
    bool otwarte = plik.is_open();
    if (!otwarte) {
        std::cout << "\n\n\n======================ERROR======================\n\nPlik o nazwie \"" << nazwa_pliku << "\" nie zostal otwarty\n\n=================================================\n\n\n";
    }

    //P�tla

    while (czas < CZAS_SYMULACJI + delta_t) {

        //Zapisywanie wynik�w w pliku tekstowym

        if (otwarte) {
            plik << czas << " " << alpha << " " << omega << " " << "\n";
        }

        //Metoda RK4 (wykonywana tylko je�eli klocek jeszcze si� nie zatrzyma�)

        if (!stop) {
            double k1_alpha = f_alpha(alpha, omega);
            double k1_omega = f_omega(alpha, omega, f);

            double k2_alpha = f_alpha(alpha + 0.5 * k1_alpha * delta_t, omega + 0.5 * k1_omega * delta_t);
            double k2_omega = f_omega(alpha + 0.5 * k1_alpha * delta_t, omega + 0.5 * k1_omega * delta_t, f);

            double k3_alpha = f_alpha(alpha + 0.5 * k2_alpha * delta_t, omega + 0.5 * k2_omega * delta_t);
            double k3_omega = f_omega(alpha + 0.5 * k2_alpha * delta_t, omega + 0.5 * k2_omega * delta_t, f);

            double k4_alpha = f_alpha(alpha + k3_alpha * delta_t, omega + k3_omega * delta_t);
            double k4_omega = f_omega(alpha + k3_alpha * delta_t, omega + k3_omega * delta_t, f);

            alpha += delta_t / 6 * (k1_alpha + 2 * k2_alpha + 2 * k3_alpha + k4_alpha);
            omega += delta_t / 6 * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);
        }

        czas += delta_t;

        //Sprawdzamy czy klocek si� zatrzyma�
        //Klocek si� zatrzyma je�li:
        //1. si�y tarcia i zsuwania s� sobie przeciwne
        //2. si�a tarcia jest wi�ksza od si�y zsuwania
        //3. wykonanie nast�pnej iteracji symulacji odwr�ci�oby pr�dko�� klocka
        //(tu wykorzystujemy metod� Eulera, gdy� nie jest wymagana dok�adno�� metody RK4 w pojedy�czej iteracji)

        double tarcie = moment_sily_tarcia(alpha, omega, f);
        double zsuwanie = moment_sily_zsuwajacej(alpha);

        if (!stop &&
            tarcie * zsuwanie < 0 &&
            abs(tarcie) > abs(zsuwanie) &&
            omega * (omega + delta_t * (zsuwanie + tarcie)) < 0) {

            //klocek si� zatrzyma�
            stop = true;
            omega = 0;
            czas_zatrzymania = czas;

        }

    }

    plik.close();

    std::cout << "Symulacja nr: " << numer_symulacji << "\nPozycja koncowa: alfa = " << alpha;
    if (stop)
        std::cout << "\nKlocek zatrzymal sie po " << czas_zatrzymania << "s\n";
    else
        std::cout << "\nKlocek nie zatrzymal sie w trakcie symulacji\n";
    std::cout << "Wyniki symulacji zapisane w pliku: " << nazwa_pliku << "\n\n";

}

//Prze�adowanie funkcji 'symulacja' w wypadku, gdy nie chcemy zapisywa� wynik�w
//R�ni si� jedynie brakiem zapisu wynik�w do pliku tekstowego oraz brakiem informacji wy�wietlanych w konsoli
double symulacja(double delta_t, double f) {

    double alpha = ALPHA;
    double omega = ALPHA_DOT;
    double czas = 0;

    double czas_zatrzymania = -1;
    bool stop = false;

    while (czas < CZAS_SYMULACJI + delta_t) {

        if (!stop) {
            double k1_alpha = f_alpha(alpha, omega);
            double k1_omega = f_omega(alpha, omega, f);

            double k2_alpha = f_alpha(alpha + 0.5 * k1_alpha * delta_t, omega + 0.5 * k1_omega * delta_t);
            double k2_omega = f_omega(alpha + 0.5 * k1_alpha * delta_t, omega + 0.5 * k1_omega * delta_t, f);

            double k3_alpha = f_alpha(alpha + 0.5 * k2_alpha * delta_t, omega + 0.5 * k2_omega * delta_t);
            double k3_omega = f_omega(alpha + 0.5 * k2_alpha * delta_t, omega + 0.5 * k2_omega * delta_t, f);

            double k4_alpha = f_alpha(alpha + k3_alpha * delta_t, omega + k3_omega * delta_t);
            double k4_omega = f_omega(alpha + k3_alpha * delta_t, omega + k3_omega * delta_t, f);

            alpha += delta_t / 6 * (k1_alpha + 2 * k2_alpha + 2 * k3_alpha + k4_alpha);
            omega += delta_t / 6 * (k1_omega + 2 * k2_omega + 2 * k3_omega + k4_omega);
        }

        czas += delta_t;

        double tarcie = moment_sily_tarcia(alpha, omega, f);
        double zsuwanie = moment_sily_zsuwajacej(alpha);

        if (!stop &&
            tarcie * zsuwanie < 0 &&
            abs(tarcie) > abs(zsuwanie) &&
            omega * (omega + delta_t * (zsuwanie + tarcie)) < 0) {

            //stop
            stop = true;
            omega = 0;
            czas_zatrzymania = czas;

        }

    }

    return alpha;
}

//Funkcja wyznaczaj�ca f
//Zwraca i wypisuje warto�� f w konsoli
double wyznacz_f(double a, double b, double dokladnosc, double dt) {

    //Metoda r�wnego podzia�u (ang. bisection method)

    double x = (a + b) / 2;

    do {

        double wartosc = symulacja(dt, x);

        if (wartosc == 0)
            break;

        if (wartosc * symulacja(dt, a) < 0) {
            b = x;
        } else {
            a = x;
        }

        x = (a + b) / 2;

    } while (abs(a - x) > dokladnosc);

    std::cout << "\n\nSzukana wartosc wspolczynnika tarcia wynosi:" << x << " +- " << dokladnosc << "\n\n";
    return x;

}

int main() {

    //Krok czasowy. Symulacje nale�y powt�rzy� dla kroku zmniejszonego pi�ciokrotnie.
    double delta_t = 0.001;
    //double delta_t = 0.0005;

    symulacja(1, delta_t, 0, "T4_s1.txt"); //przypadek a)
    symulacja(2, delta_t, 0.02, "T4_s2.txt"); //przypadek b)
    symulacja(3, delta_t, 0.1, "T4_s3.txt"); //przypadek c)
    symulacja(4, delta_t, 0.5, "T4_s4.txt"); //przypadek d)

    double dokladnosc = 0.01; //Delta - oznacza dok�adno�� wyznaczanej warto�ci wsp�czynnika tarcia
    double a = 0.5;
    double b = 1;

    //Symulacja przeprowadzona dla wyznaczonego wsp�czynnika tarcia w celu upewnienia si�, �e klocek faktycznie zatrzymuje si� dok�adnie na �rodku
    symulacja(3, delta_t, wyznacz_f(a, b, dokladnosc, delta_t), "T_4_tarcie.txt");

}