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


W plikach dostępny jest powyższy plik do pobrania. 
Program wykonuje wybraną z powyższych transformacji dla dowolnej wybranej elispoidy (wgs84,grs80 lub Krasowskiego). Żeby wykorzystać możliwości programu trzeba przede wszystkich wpisać wymagane dane do pliku tekstowego. Gdy już to zrobimy, upewniamy się, że wszystkie pliki znajdują się w jednym folderze i otwieramy wiersz poleceń. Piszemy "cd"  i wklejamy ścieżkę do naszego folderu z plikami, a następnie wpisujemy python kod.py --transform  --model , po "transform" wpisujemy wybraną przez nas transformację, a po "model" wybraną elipsoidę. Poglądowo może to wygladać np. tak: 

![image](https://github.com/Em2510/Projekt_1/assets/129061319/d67fc9be-cdba-4199-9404-8db9cd056509)

Następnie powyższe wyniki zapisują się w naszym folderze jako plik tekstowy

![image](https://github.com/Em2510/Projekt_1/assets/129061319/e318d53a-1aba-47f1-b962-63c42d60c397)

Jednocześnie kod "czuwa" nad wybranej nieistniejącej w kontekście naszego kodu elipsoidy:

![image](https://github.com/Em2510/Projekt_1/assets/129061319/dceea54f-4b03-4ff9-8032-01e800511117)


Przesyłam jeszcze przykładowe użycia kodu: 




