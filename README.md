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

Projekt został rozdzielony na dwa pliku w celu przejrzystości kodu. Program składa się z pliku definicje.py oraz kod.py. Dodatkowo wymagany jest plik dane.txt z danymi, napisany w poniższy sposób:
![image](https://github.com/Em2510/Projekt_1/assets/129061319/5bf589e7-98be-48d9-b134-85f3fbc1ddda)


