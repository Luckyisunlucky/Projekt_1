import numpy as np
from projekt1 import *

p = Transformacje(model = 'wgs84')

plik = 'wsp_inp.txt'
#odczyt z pliku
tablica = np.genfromtxt(plik, delimiter = ',', skip_header = 4)

w, r = np.shape(tablica)

wynik = np.zeros((w,10))

i= 0

for wiersz in tablica:
     f, l, h = p.xyz2plh(wiersz[0], wiersz[1], wiersz[2])
     x, y, z = p.flh2XYZ(f, l, h)
     xgk, ygk = p.fl2xy(f, l)
     x2000, y2000 = p.u2000(xgk, ygk)
     x1992, y1992 = p.u1992(xgk, ygk)
     
     wynik[i, 0] = f
     wynik[i, 1] = l
     wynik[i, 2] = h
     wynik[i, 3] = x
     wynik[i, 4] = y
     wynik[i, 5] = z
     wynik[i, 6] = x2000
     wynik[i, 7] = y2000
     wynik[i, 8] = x1992
     wynik[i, 9] = y1992           
     i+=1
print(wynik)
    
np.savetxt('wsp_out.txt', wynik, delimiter = ',', fmt = ['%10.7f', '%10.7f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f','%10.3f'], 
           header = 'Mateusz Brzeziński\n Wydział Geodezji i Kartografii\n Zamiana współrzędnych geodezyjnych')