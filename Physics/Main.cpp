// Kod programu stanowi¹cy czêœæ rozwi¹zania zadania T4 I stopnia 72. Olimpiady Fizycznej
// Autor programu: Maciej Zgierski
// Listopad 2022

// Program wyznacza i zapisuje wyniki kilku symulacji do plików o nazwie T4_s<nr>.txt oraz T4_wsp_tarcia.txt
// Wypisuje tak¿e u¿yteczne informacje o symulacjach w konsoli

// W celu wykonania wykresów nale¿y pos³u¿yæ siê oddzielnym programem (np. arkuszem kalkulacyjnym)

#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

//Sta³e

const double g = 9.81;
const double L = 0.1;

const double CZAS_SYMULACJI = 10;

//Pocz¹tkowe wartoœci symulacji

const double ALPHA = -M_PI / 2;
const double ALPHA_DOT = 0.0;

//Funkcja obliczaj¹ca przyspieszenie k¹towe wywo³ane momentem si³y zsuwaj¹cej dzia³aj¹ce na cia³o w danej chwili
double moment_sily_zsuwajacej(double alpha) {
    return - g / L * sin(alpha);
}

//Funkcja obliczaj¹ca przyspieszenie k¹towe wywo³ane momentem si³y tarcia dzia³aj¹cej na poruszaj¹cy siê klocek w danej chwili
double moment_sily_tarcia(double alpha, double omega, double f) {
    if (omega == 0) //kiedy klocek siê nie porusza tarcie kinetyczne nie jest obecne
        return 0;
    if (omega > 0) //sprawdzamy w któr¹ stronê porusza siê klocek w celu wyznaczenia kierunku dzia³ania tarcia
        return - f * (omega * omega + g / L * cos(alpha));
    return f * (omega * omega + g / L * cos(alpha));
}

//Równanie ró¿niczkowe dla alfy (f_alfa - patrz równanie (7))
double f_alpha(double alpha, double omega) {
    return omega;
}

//Równanie ró¿niczkowe dla omegi (f_omega - patrz równanie (7))
double f_omega(double alpha, double omega, double f) {
    return moment_sily_zsuwajacej(alpha) + moment_sily_tarcia(alpha, omega, f);
}

void symulacja(int numer_symulacji, double delta_t, double f, const char* nazwa_pliku) {

    //Przygotowanie symulacji
    double alpha = ALPHA;
    double omega = ALPHA_DOT;
    double czas = 0;

    double czas_zatrzymania = -1;
    bool stop = false; //stop jest identyfikatorem zatrzymania siê klocka

    //Otwarcie pliku tekstowego

    std::fstream plik;
    plik.open(nazwa_pliku, std::ios::out);
    bool otwarte = plik.is_open();
    if (!otwarte) {
        std::cout << "\n\n\n======================ERROR======================\n\nPlik o nazwie \"" << nazwa_pliku << "\" nie zostal otwarty\n\n=================================================\n\n\n";
    }

    //Pêtla

    while (czas < CZAS_SYMULACJI + delta_t) {

        //Zapisywanie wyników w pliku tekstowym

        if (otwarte) {
            plik << czas << " " << alpha << " " << omega << " " << "\n";
        }

        //Metoda RK4 (wykonywana tylko je¿eli klocek jeszcze siê nie zatrzyma³)

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

        //Sprawdzamy czy klocek siê zatrzyma³
        //Klocek siê zatrzyma jeœli:
        //1. si³y tarcia i zsuwania s¹ sobie przeciwne
        //2. si³a tarcia jest wiêksza od si³y zsuwania
        //3. wykonanie nastêpnej iteracji symulacji odwróci³oby prêdkoœæ klocka
        //(tu wykorzystujemy metodê Eulera, gdy¿ nie jest wymagana dok³adnoœæ metody RK4 w pojedyñczej iteracji)

        double tarcie = moment_sily_tarcia(alpha, omega, f);
        double zsuwanie = moment_sily_zsuwajacej(alpha);

        if (!stop &&
            tarcie * zsuwanie < 0 &&
            abs(tarcie) > abs(zsuwanie) &&
            omega * (omega + delta_t * (zsuwanie + tarcie)) < 0) {

            //klocek siê zatrzyma³
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

//Prze³adowanie funkcji 'symulacja' w wypadku, gdy nie chcemy zapisywaæ wyników
//Ró¿ni siê jedynie brakiem zapisu wyników do pliku tekstowego oraz brakiem informacji wyœwietlanych w konsoli
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

//Funkcja wyznaczaj¹ca f
//Zwraca i wypisuje wartoœæ f w konsoli
double wyznacz_f(double a, double b, double dokladnosc, double dt) {

    //Metoda równego podzia³u (ang. bisection method)

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

    //Krok czasowy. Symulacje nale¿y powtórzyæ dla kroku zmniejszonego piêciokrotnie.
    double delta_t = 0.001;
    //double delta_t = 0.0005;

    symulacja(1, delta_t, 0, "T4_s1.txt"); //przypadek a)
    symulacja(2, delta_t, 0.02, "T4_s2.txt"); //przypadek b)
    symulacja(3, delta_t, 0.1, "T4_s3.txt"); //przypadek c)
    symulacja(4, delta_t, 0.5, "T4_s4.txt"); //przypadek d)

    double dokladnosc = 0.01; //Delta - oznacza dok³adnoœæ wyznaczanej wartoœci wspó³czynnika tarcia
    double a = 0.5;
    double b = 1;

    //Symulacja przeprowadzona dla wyznaczonego wspó³czynnika tarcia w celu upewnienia siê, ¿e klocek faktycznie zatrzymuje siê dok³adnie na œrodku
    symulacja(3, delta_t, wyznacz_f(a, b, dokladnosc, delta_t), "T_4_tarcie.txt");

}