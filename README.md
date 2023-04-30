Projekt 1: 
1. Program służy do wykonywania transformacji:
- XYZ -> BLH
Transformacja współrzędnych ortokartezjańskich na współrzędne geodezyjne. 
- BLH -> XYZ 
Transformacja współrzędnych geodezyjnych na współrzędne ortokartezjańskie
- XYZ -> NEU
Transformacja współrzędnych ortokartezjańskich na współrzędne topocentryczne
- BL -> 2000
Trzansformacja współrzednych geodezyjnych na współrzedne w układzie 2000
- BL -> 1992
Transformacja współrzędnych geodezyjnych na współrzędne w układzie 1992

* Transformacje w programie mogą zostac wykonane w oparciu na elipsoidy wgs84, grs80 albo Krasowskiego(kras)

2. Aby program działał poprawnie, należy zainstalować python w wersji co najmniej 3.9 oraz biblioteki numpy, argparse oraz scipy.

3. Program został napisany dla systemu, który wspiera Python w wersji 3.9 albo wyższej oraz ma zainstalowane wymagane biblioteki.

4. żeby skorzystać z programu należy w folderze z programem utworzyć plik tekstowy o nazwie dane.txt.
Program jest zbudowany w ten sposób żeby wybierał potrzebne dane do konkretnego zadania. Kolejno kolumny to: X, Y, Z, X0, Y0, Z0, fi, lam, h.
Przykładowo zadanie zostało wykonane na danych:
3937512.213,2316765.301,4345641.066,5004565.672,1241213.654,1490434.615,1.9875,0.4006,873.2
wynikami transformacji są:
XYZ2FLH= [(43.76185154872513, 30.471832446824944, -62715.849972978234)]
FLH2XYZ= (-2514997.0844853334, 5682008.7657792345, 1433995.8258535683)
XYZ2neu [array([-1067053.459,  1075551.647,  2855206.451])]
FL2PL1992= [(9483624293570.684, 17406141.295948178)]
FL2PL2000= [(9489540830568.686, 15812827.431141818)]

5. Występuje problem przy zmianie elipsoidy: program bierze tylko jedną elipsoidę
