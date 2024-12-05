// ------------------------------------------------------------------------------------------------
// Programm "Morphologie"
//   Das Programm dient zur Korrektur inhomogen ausgeleuchteter Bilder
//   - Der Name des Eingangsbildes wird (ohne ".bmp"-Endung)
//     dem Programm als Argument uebergeben.
//     (Bei Code::Blocks werden die Argumente unter
//      "Project"->"Set programs' arguments..." angegeben)
//
// B. Lang, HS Osnabrueck
// Version Dezember 2016 (RH)
// ------------------------------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <string>
#include <limits>
#include <cstring>    // strerror
#include "Img.h"
#include "BmpRead.h"
#include "ConvImg.h"
#include "BmpWrite.h"

// ---------------------------------------------------------------------------
// Groesse des Strukturierenden Elementes festlegen
// ---------------------------------------------------------------------------
#define SE_GROESSE 3

// ---------------------------------------------------------------------------
// Vektor mit relativen Positionen der Pixel eines quadratischen SEs fuellen
// ---------------------------------------------------------------------------
// Parameter:
// [in]  Diameter : Seitenlaenge des quadratischen SEs in Pixel
// Return:
// Vektor mit relativen Positionen des SEs
// ---------------------------------------------------------------------------
vector<Position> create_square_SE(int Diameter)
{
	vector<Position> ImageWindow;
	int Radius = Diameter / 2;
	for (int dy = -Radius; dy <= Radius; dy++) {
		for (int dx = -Radius; dx <= Radius; dx++) {
			ImageWindow.push_back(Position(dx, dy));
		}
	}
	return ImageWindow;
}

// ---------------------------------------------------------------------------
// Aufgabe 3: Code erstellen
// Vektor mit relativen Positionen der Pixel eines nahezu runden SEs fuellen
// ---------------------------------------------------------------------------
// Parameter:
// [in]  Diameter : Durchmesser des nahezu runden SEs in Pixel
// Return:
// Strukturierendes Element (Vektor mit Positionen)
// ---------------------------------------------------------------------------
vector<Position> create_round_SE(int Diameter)
{
	vector<Position> ImageWindow;


    int Radius = Diameter / 2;
    int RadiusSquared = Radius * Radius;


    for (int y = -Radius; y <= Radius; ++y) {
        for (int x = -Radius; x <= Radius; ++x) {
            if (x * x + y * y <= RadiusSquared) {
                ImageWindow.push_back(Position{x, y});
            }
        }
    }


    return ImageWindow;\
}

// ---------------------------------------------------------------------------
// Aufgabe 4: Code erstellen
// SE spiegeln
// ---------------------------------------------------------------------------
// Parameter:
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// gespiegeltes Strukturierendes Element (Vektor mit Positionen)
// ---------------------------------------------------------------------------
vector<Position> mirror_SE(const vector<Position> &ImageWindow)
{
	vector<Position> MirroredImageWindow(ImageWindow.size());

    for (const auto &pos : ImageWindow) {
        MirroredImageWindow.push_back(Position{-pos.get_x(), -pos.get_y()});
    }

	return MirroredImageWindow;
}

// ---------------------------------------------------------------------------
// Aufgabe 2: Code erstellen
// Erosion
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// erodiertes Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> erode(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
	Img<Pixel> eroded(src.Width(), src.Height());
    const unsigned int Width = src.Width();
    const unsigned int Height = src.Height();


    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            Pixel min_value = std::numeric_limits<Pixel>::max();
            for (const auto &pos : ImageWindow) {
                int nx = x + pos.get_x();
                int ny = y + pos.get_y();

                if (nx >= 0 && nx < Width && ny >= 0 && ny < Height) {
                    min_value = std::min(min_value, src[ny][nx]);
                }
            }
            eroded[y][x] = min_value;
        }
    }

	return eroded;
}

// ---------------------------------------------------------------------------
// Aufgabe 2: Code erstellen
// Dilation
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// dilatiertes Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> dilate(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
    const unsigned int Width = src.Width();
    const unsigned int Height = src.Height();
	Img<Pixel> dilated(src.Width(), src.Height());


    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            Pixel max_value = std::numeric_limits<Pixel>::min();
            for (const auto &pos : ImageWindow) {
                int nx = x + pos.get_x();
                int ny = y + pos.get_y();

                if (nx >= 0 && nx < Width && ny >= 0 && ny < Height) {
                    max_value = std::max(max_value, src[ny][nx]);
                }
            }
            dilated[y][x] = max_value;
        }
    }

	return dilated;
}

// ---------------------------------------------------------------------------
// Aufgabe 4: Code erstellen
// Opening
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// geoeffnetes Bild Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> opening(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
	Img<Pixel> opened;

    Img<Pixel> eroded = erode(src, ImageWindow);


    std::vector<Position> mirrored_SE = mirror_SE(ImageWindow);
     opened = dilate(eroded, mirrored_SE);

	return opened;
}

// ---------------------------------------------------------------------------
// Aufgabe 4: Code erstellen
// Closing
// ---------------------------------------------------------------------------
// Parameter:
// [in]  src         : Eingangsbild
// [in]  ImageWindow : Strukturierendes Element (Vektor mit Positionen)
// Return:
// geschlossenes Bild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> closing(const Img<Pixel> &src, const vector<Position> &ImageWindow)
{
	Img<Pixel> closed;

    Img<Pixel> dilated = dilate(src, ImageWindow);


    std::vector<Position> mirrored_SE = mirror_SE(ImageWindow);
    closed = erode(dilated, mirrored_SE);

	return closed;
}

// ---------------------------------------------------------------------------
// Aufgabe 5: Code erstellen
// Bilddifferenz
// ---------------------------------------------------------------------------
// Parameter:
// [in]  l : Bild links vom Operator '-'
// [in]  r : Bild rechts vom Operator '-'
// Return:
// Differenzbild
// ---------------------------------------------------------------------------
template<typename Pixel>
Img<Pixel> operator-(const Img<Pixel> &l, const Img<Pixel> &r)
{
	Img<Pixel> d(l.Width(), l.Height());


    for (unsigned int y = 0; y < l.Height(); ++y) {
        for (unsigned int x = 0; x < l.Width(); ++x) {
            // Berechnung der Differenz der Pixelwerte (Clamping bei Bedarf)
            d[y][x] = std::max(0, l[y][x] - r[y][x]);
        }
    }

	return d;
}


// ---------------------------------------------------------------------------
// Aufgabe 1: Code erstellen
// Optimale Schwelle
// ---------------------------------------------------------------------------
// Parameter:
// [in]  gray_image : Grauwertbild
// Return:
// Binaerbild
// ---------------------------------------------------------------------------
Img<bool> optimal_threshold(const Img<unsigned char> &gray_image)
{
	const unsigned int Width = gray_image.Width();
	const unsigned int Height  = gray_image.Height();
	Img<bool> binary_image(Width, Height);

    // Histogramm berechnen
    std::vector<unsigned int> histogram(256, 0);
    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            ++histogram[gray_image[y][x]];
        }
    }

    // Gesamtanzahl der Pixel
    const unsigned int total_pixels = Width * Height;

    // Otsu-Methode: Optimalen Schwellenwert berechnen
    double sum = 0.0;
    for (unsigned int i = 0; i < 256; ++i) {
        sum += i * histogram[i];
    }

    double sum_background = 0.0;
    unsigned int weight_background = 0;
    unsigned int weight_foreground = 0;

    double max_variance = 0.0;
    unsigned int optimal_threshold = 0;

    for (unsigned int t = 0; t < 256; ++t) {
        weight_background += histogram[t];
        if (weight_background == 0) continue;

        weight_foreground = total_pixels - weight_background;
        if (weight_foreground == 0) break;

        sum_background += t * histogram[t];

        double mean_background = sum_background / weight_background;
        double mean_foreground = (sum - sum_background) / weight_foreground;

        double variance_between = weight_background * weight_foreground *
                                  std::pow(mean_background - mean_foreground, 2);

        if (variance_between > max_variance) {
            max_variance = variance_between;
            optimal_threshold = t;
        }
    }

    // Binärbild erzeugen
    for (unsigned int y = 0; y < Height; ++y) {
        for (unsigned int x = 0; x < Width; ++x) {
            binary_image[y][x] = gray_image[y][x] > optimal_threshold;
        }
    }

	return binary_image;
}

// ------------------------------------------------------------------------------------------------
// Hauptprogramm
// ------------------------------------------------------------------------------------------------
// Parameter:
// [in] argv[1] : Name des Eingangsbildes
// ------------------------------------------------------------------------------------------------
int main(int argc, char*argv[])
{
	string filename;

	cout << "BV-Praktikum: Morphologie" << endl << endl;

//	if (argc < 2) {
//		cerr << "Dateiname nicht angegeben" << endl;
//		return -1;
//	}
	string Bildname("C:/Users/quint/CLionProjects/untitled/Inhomogen-1");

	// Bild einlesen und nach BYTE konvertieren
	Img<RGB_Pixel> rgb;
	try {
		filename = Bildname + ".bmp";
		BmpRead(filename.c_str()) >> rgb;
		cout << "Lese " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Lesen von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// Quellbild in Luminanzbild konvertieren
	// (Luminanzbild im Wertebereich 0.0 bis 1.0)
	Img<unsigned char> src(ConvImg<unsigned char, RGB_Pixel>(rgb));
	src.Margin_Mirror();
	try {
		filename = Bildname + "_src.bmp";
		BmpWrite(filename.c_str(), src);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// 1. Aufgabe: Optimale Schwelle
	// --------------------------------------------------------------------------------
	try {
		Img<bool> erstes_Binaerbild = optimal_threshold(src);
		filename = Bildname + "_binaer.bmp";
		BmpWrite(filename.c_str(), erstes_Binaerbild);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// 2. Aufgabe: Erosion und Dilation mit quadratischen SE
	// --------------------------------------------------------------------------------

	// Quadratisches Bildfenster aufbauen
	vector<Position> quadratisches_Bildfenster;
	quadratisches_Bildfenster = create_square_SE(SE_GROESSE);

	try {
		// Erosion durchfuehren
		Img<unsigned char> eroded = erode<unsigned char>(src, quadratisches_Bildfenster);
		filename = Bildname + "_eroded_sq.bmp";
		BmpWrite(filename.c_str(), eroded);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	try {
		// Dilation durchfuehren
		Img<unsigned char> dilated = dilate<unsigned char>(src, quadratisches_Bildfenster);
		filename = Bildname + "_dilated_sq.bmp";
		BmpWrite(filename.c_str(), dilated);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// 3. Aufgabe: Erosion und Dilation mit rundem SE
	// --------------------------------------------------------------------------------
	// Rundes Bildfenster aufbauen
	vector<Position> rundes_Bildfenster;
	rundes_Bildfenster = create_round_SE(SE_GROESSE);

	try {
		// Erosion durchfuehren
		Img<unsigned char> eroded = erode<unsigned char>(src, rundes_Bildfenster);
		filename = Bildname + "_eroded_rd.bmp";
		BmpWrite(filename.c_str(), eroded);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	try {
		// Dilation durchfuehren
		Img<unsigned char> dilated = dilate<unsigned char>(src, rundes_Bildfenster);
		filename = Bildname + "_dilated_rd.bmp";
		BmpWrite(filename.c_str(), dilated);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// 4. Aufgabe: Opening und Closing mit rundem SE
	// --------------------------------------------------------------------------------

	Img<unsigned char> opened;
	try {
		// Opening durchfuehren
		opened = opening<unsigned char>(src, rundes_Bildfenster);
		filename = Bildname + "_opened.bmp";
		BmpWrite(filename.c_str(), opened);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	Img<unsigned char> closed;
	try {
		// Closing durchfuehren
		closed = closing<unsigned char>(src, rundes_Bildfenster);
		filename = Bildname + "_closed.bmp";
		BmpWrite(filename.c_str(), closed);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	// --------------------------------------------------------------------------------
	// 5. Aufgabe: Zylinderhut-Bilder berechnen und optimale Schwelle anwenden
	// --------------------------------------------------------------------------------

	Img<unsigned char> WZH;
	try {
		// WZH berechnen
		WZH = src - opened;
		filename = Bildname + "_WZH.bmp";
		BmpWrite(filename.c_str(), WZH);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	try {
		Img<bool> WZH_bool = optimal_threshold(WZH);
		filename = Bildname + "_WZH_bool.bmp";
		BmpWrite(filename.c_str(), WZH_bool);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	Img<unsigned char> SZH;
	try {
		// SZH berechnen
		SZH = closed - src;
		filename = Bildname + "_SZH.bmp";
		BmpWrite(filename.c_str(), SZH);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	try {
		Img<bool> SZH_bool = optimal_threshold(SZH);
		filename = Bildname + "_SZH_bool.bmp";
		BmpWrite(filename.c_str(), SZH_bool);
		cout << "Schreibe " << filename << endl;
	} catch (const char * s) {
		cerr << "Fehler beim Schreiben von " << filename << ": " << strerror(errno) << endl;
		return -1;
	}

	return 0;
}
