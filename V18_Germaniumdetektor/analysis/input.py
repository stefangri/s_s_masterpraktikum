import numpy as np
import sys
from collections import OrderedDict
from pathlib import Path
import os.path

print('Dateinamen eingeben (ohne .txt): ')
dateiname = input()

pfad_zu_datei='data/{}.txt'.format(dateiname)
messdaten_txt=Path(pfad_zu_datei)

if messdaten_txt.is_file():
    print('Die Datei scheint schon erstellt zu sein')
    print('Möchtest du daran weiterschreiben (w)')
    print('oder neu beginnen [alle daten gehen verloren](n)')
    eingabe=input()
    print('\n')

    if eingabe == 'w':
        dateiname_txt=dateiname+'.txt'
        data_arra = np.genfromtxt(dateiname_txt,unpack=True,names=True)
        head=data_arra.dtype.names
        data_arra = np.genfromtxt(dateiname_txt,unpack=True)
        p = data_arra.shape[0]
        data_list=data_arra.tolist()
        print(data_list)
        print('Zeilenweise die Messwerte eingeben: ')
        u = 'a'
        while u != 'x':
            for k in range(p):
                u = input()
                if u == 'x':
                    break
                elif u == 'pop':
                    for n in range(k):
                        data_list[n].pop()
                    print(data_list[n])
                    break
                elif u == 'popl':
                    if k!=0:
                        print('Fehler, Zeile wird zurückgesetzt')
                        for n in range(int(k)):
                            data_list[n].pop()
                        print(data_list)
                        break
                    for n in range(len(data_list)):
                        data_list[n].pop()
                        print(data_list[n])
                    break
                try:
                    float(u)
                except ValueError:
                    print('Fehler')
                    print('Sicherheitsspeicherung ? Ja(j) oder Nein (n)')
                    eingabe_sicher=input()
                    if eingabe_sicher=='j':
                        for n in range(k):
                            data_list[n].pop()
                        data_array=(data_list)
                        np.savetxt(dateiname + '.txt', np.column_stack(data_array),header = str(head)) ##Hier eine anpassung
                        print('Sicherheitskopie erstellt')
                    else:
                        print('Programm beendet')

                data_list[k].append(float(u))
                print(data_list[k])
            data_array=(data_list)
        np.savetxt(f'data/{dateiname_txt}', np.column_stack(data_array),header=str(head))



            #np.savetxt(dateiname + '.txt', np.column_stack(data.values()), header = head) ##Hier eine An


    else:
        print('Anzahl der Messgrößen eingeben')
        p = int(input())
        data = OrderedDict()
        head = ""
        for i in range(0, p):
            print(i+1, 'te Messgröße eingeben')
            v = input()
            data[i] = []
            head += v
            head += " "
        print(head)
        print('Zeilenweise die Messwerte eingeben: ')
        u = 'a'
        while u != 'x':
            for k in data:
                u = input()
                if u == 'x':
                    break
                elif u=='pop':
                    for n in range(int(k)):
                        data[n].pop()
                    print(data)
                    break
                elif u=='popl':
                    if k!=0:
                        print('Fehler, Zeile wird zurückgesetzt')
                        for n in range(int(k)):
                            data[n].pop()
                        print(data)
                        break
                    else:
                        for n in range(p):
                            data[n].pop()
                        print(data)
                        break
                try:
                    float(u)
                except ValueError:
                    print('Fehler')
                    print('Sicherheitsspeicherung ? Ja(j) oder Nein (n)')
                    eingabe_sicher=input()
                    if eingabe_sicher=='j':
                        for n in range(int(k)):
                            data[n].pop()
                        np.savetxt(f'data/{dateiname_txt}', np.column_stack(data.values()), header = head)
                        print('Sicherheitskopie erstellt')
                    else:
                         print('Programm beendet')
                data[k].append(float(u))
                print(data[k])
            np.savetxt(f'data/{dateiname_txt}', np.column_stack(data.values()), header = head)

else:
    print('Anzahl der Messgrößen eingeben')
    p = int(input())
    data = OrderedDict()
    head = ""
    for i in range(0, p):
        print(i+1, 'te Messgröße eingeben')
        v = input()
        data[i] = []
        head += v
        head += " "

    print(head)
    print('Zeilenweise die Messwerte eingeben: ')
    u = 'a'
    while u != 'x':
        for k in data:
            u = input()
            if u == 'x':
                break
            elif u=='pop':
                for n in range(int(k)):
                    data[n].pop()
                print(data)
                break
            elif u=='popl':
                if k!=0:
                    print('Fehler, Zeile wird zurückgesetzt')
                    for n in range(int(k)):
                        data[n].pop()
                    print(data)
                    break
                else:
                    for n in range(p):
                        data[n].pop()
                    print(data)
                    break
            try:
                float(u)
            except ValueError:
                print('fehlerhafte Eingabe: String')
                print('Sicherheitsspeicherung ? Ja(j) oder Nein (n)')
                eingabe_sicher=input()
                if eingabe_sicher=='j':
                    for n in range(int(k)):
                        data[n].pop()
                    np.savetxt(f'data/{dateiname}.txt', np.column_stack(data.values()), header = head)
                    print('Sicherheitskopie erstellt')
                else:
                    print('Programm beendet')

            data[k].append(float(u))
            print(data[k])

    np.savetxt(f'data/{dateiname}.txt', np.column_stack(data.values()), header = head)
