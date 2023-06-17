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

JEDNOSTKI 
- Współrzędna X punktu wyrażona w metrach
- Współrzędna Y punktu wyrażona w metrach
- Współrzędna Z punktu wyrażona w metrach
- Współrzędna X0 odbiornika wyrażona w metrach
- Współrzędna Y0 odbiornika wyrażona w metrach
- Współrzędna Z0 odbiornika wyrażona w metrach
- Współrzędna Fi punktu wyrażona w stopniach
- Współrzędna Lambda punktu wyrażona w stopniach
- Współrzędna h punktu wyrażona w metrach
- Południk osiowy strefy odwzorowawczej wyrażony w stopniach

Wyniki będą analogicznie w tych samych jednostkach. 


* Transformacje w programie mogą zostac wykonane w oparciu na elipsoidy wgs84, grs80 albo Krasowskiego(kras)

2. Aby program działał poprawnie, należy zainstalować python w wersji co najmniej 3.9 oraz biblioteki numpy, argparse oraz scipy.

3. Program został napisany dla systemu, który wspiera Python w wersji 3.9 albo wyższej oraz ma zainstalowane wymagane biblioteki.

Projekt został rozdzielony na dwa pliku w celu przejrzystości kodu. Program składa się z pliku definicje.py oraz kod.py. Dodatkowo wymagany jest plik dane.txt z danymi, napisany w poniższy sposób:
![image](https://github.com/Em2510/Projekt_1/assets/129061319/5bf589e7-98be-48d9-b134-85f3fbc1ddda)



Program wykonuje wybraną z powyższych transformacji dla dowolnej wybranej elispoidy (wgs84,grs80 lub Krasowskiego). Żeby wykorzystać możliwości programu trzeba przede wszystkich wpisać wymagane dane do pliku tekstowego. Gdy już to zrobimy, upewniamy się, że wszystkie pliki znajdują się w jednym folderze i otwieramy wiersz poleceń. Piszemy "cd"  i wklejamy ścieżkę do naszego folderu z plikami, a następnie wpisujemy python kod.py --transform  --model , po "transform" wpisujemy wybraną przez nas transformację, a po "model" wybraną elipsoidę. Poglądowo może to wygladać np. tak: 

![image](https://github.com/Em2510/Projekt_1/assets/129061319/d67fc9be-cdba-4199-9404-8db9cd056509)

Następnie powyższe wyniki zapisują się w naszym folderze jako plik tekstowy

![image](https://github.com/Em2510/Projekt_1/assets/129061319/e318d53a-1aba-47f1-b962-63c42d60c397)

Jednocześnie kod "czuwa" nad wybranej nieistniejącej w kontekście naszego kodu elipsoidy:

![image](https://github.com/Em2510/Projekt_1/assets/129061319/dceea54f-4b03-4ff9-8032-01e800511117)


Przesyłam jeszcze przykładowe użycia kodu. 
Funkcja przedstawia transformację współrzędnych geocentrycznych odbiornika (X, Y, Z) do układu topocentrycznego aby uzyskać współrzędne topocentryczne odbiornika. Zadanie wykorzystuje także współrzędne geocentryczne anteny (X0, Y0, Z0). W przypadku tej transformacji istnieją dwa punkty wejściowe, ponieważ jest to przekształcenie dwupunktowe.

![image](https://github.com/Em2510/Projekt_1/assets/129061319/974a9b7d-6560-4146-88b8-b7fc7e75b7a7)


Kolejny przykład:

![image](https://github.com/Em2510/Projekt_1/assets/129061319/a22b218a-40f5-45b4-aead-174055c487ae)

ZNANE BŁĘDY: 
NIE ZALECAMY STOSOWANIA ELIPSOIDY KRASOWSKIEGO, PONIEWAŻ PODAJE ONA BŁĘDNE WYNIKI. 
W przypadku błędów związanych z wierszem poleceń, być może wystarczy zainstalować "python -m pip install numpy bądź "python -m pip install scipy".

