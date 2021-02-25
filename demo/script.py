#test
from flask import Flask, request, redirect, url_for, render_template
from math import sqrt,ceil

import numpy as np


app = Flask(__name__)
@app.route('/plot/')
def plot():

    zadany_procent_zarazonych = int(request.args.get("zadany_procent_zarazonych", 15))
    if(zadany_procent_zarazonych>30):
        zadany_procent_zarazonych=30
    elif(zadany_procent_zarazonych<10):
        zadany_procent_zarazonych=10
    liczba_osob=int(request.args.get("liczba_osob", 10000))
    if(liczba_osob>50000):
        liczba_osob=50000
    elif(liczba_osob<1000):
        liczba_osob=1000
    liczba_zakazonych=int(request.args.get("liczba_zakazonych", 5000))
    stratowa_liczba_zakazonych=liczba_zakazonych
    if(liczba_zakazonych>liczba_osob):
        liczba_zakazonych=liczba_osob
    elif(liczba_zakazonych<1):
        liczba_zakazonych=1
    gammaA=float(request.args.get("gammaA", 0.94))  #wspolczynnik przypadkow hospitalizowanych
    if(gammaA>0.999):
        gammaA=0.999
    elif(gammaA<0.001):
        gammaA=0.001
    gammaR=float(request.args.get("gammaR", 0.5))  #wspolczynnik wyzdrowien wsrod przypadkow hospitalizowanych
    if(gammaR>0.999):
        gammaR=0.999
    elif(gammaR<0.15):
        gammaR=0.15
    gammaI=float(request.args.get("gammaI", 0.27))  #wspolczynnik wyzdrowien poza hospitalizacja
    if(gammaI>0.999):
        gammaI=0.999
    elif(gammaI<0.15):
        gammaI=0.15
    jota=float(request.args.get("jota", 1.56))   #wzgledna zarazliwosć wirusa u przypadkow hospitalizowanych
    if(jota>3):
        jota=3
    elif(jota<0.5):
        jota=0.5
    kappa=float(request.args.get("kappa", 0.25))   #wspolczynnik zarazliwosci wirusa u osob bezposrednio narazonych na zarazenie [E]
    if(kappa>1):
        kappa=1
    elif(kappa<0.1):
        kappa=0.1
    ro=float(request.args.get("ro", 0.99)) #wspolczynnik decydujacy o przydzieleniu do klasy zarazonych [I]
    if(ro>0.999):
        ro=0.999
    elif(ro<0.5):
        ro=0.5


        #Dane poczatkowe
   
    
    ilosc_krokow=1000
    

    #Wspolczyniki PID
    kp=0.001
    ki=0.0005
    kd=0.001

    #Parametry

    betaI=0  #wspolczynnik reprodukcji wirusa w przeciętnych przypadkach
    betaP=0  #wspolczynnik reprodukcji wirusa u osob wykazujacych wyzszy stopień zarazliwosci (tzw. superroznosiciele [P])
    

    #Zmienne pomocnicze
    zadana_liczba_chorych=ceil(liczba_osob*zadany_procent_zarazonych/100)
    N=liczba_osob                                   #liczba ludności
    S=liczba_osob-liczba_zakazonych                 #podatni na zarazenie
    E=0                                             #bezposrednio narazeni na efekt
    I=ceil(liczba_zakazonych*ro)                     #zarazający i wykazujący objawy:
    P=liczba_zakazonych-I                           #superroznosiciele
    H=0                                             #hospitalizowani
    R=0                                             #przypadki ozdrowiale
    e=zadana_liczba_chorych-liczba_zakazonych     #uchyb regulacji
    u=0                                             #wartosc sterujaca
    suma_e=0
    n=0
    wykres=[]
    wykres_P=[]
    wykres_S=[]
    wykres_E=[]
    wykres_I=[]
    wykres_H=[]
    wykres_R=[]
    wykres_krokow = [i for i in range(ilosc_krokow+1)]

    while(n <= ilosc_krokow):
        #poczatek e(n)

        e_poprzednie=e
        if(e_poprzednie>1000):
            e_poprzednie=1000
        if(e_poprzednie<-1000):
            e_poprzednie=-1000
        e=zadana_liczba_chorych-liczba_zakazonych
        if(e>1000):
            e=1000
        if(e<-1000):
            e=-1000
        delta_e=abs(e_poprzednie-e)
        if(delta_e>1000):
            delta_e=1000
        if(delta_e<-1000):
            delta_e=-1000
        suma_e+=e
        if(e>1000):
            suma_e=1000
        if(e<-1000):
            e=-1000
        #koniec e(n)
        #poczatek PID
        u=kp*e+ki*suma_e+kd*delta_e
        
        #koniec PID

        #poczatek sam nie wiem czego
        BetaI=u/50
        BetaP=u/150
        #koniec no tego
        # if(n<=500):
            # print(u,n,BetaI, liczba_zakazonych)
        #poczatek liczby_zarazownych(n)
        S_p=S
        E_p=E
        I_p=I
        P_p=P
        H_p=H
        liczba_zakazonych=I+P+H
        wykres_S.append(S)
        wykres_E.append(E)
        wykres_I.append(I)
        wykres_P.append(P)
        wykres_H.append(H)
        wykres_R.append(R)
        S += (R - ceil(BetaI*I_p*S_p/N) -ceil(BetaP*P_p*S_p/N) - ceil(BetaI*jota*H_p*S_p/N))
        E += (ceil(BetaI*I_p*S_p/N) + ceil(BetaP*P_p*S_p/N) + ceil(BetaI*jota*H_p*S_p/N) - ceil(kappa*E_p*ro) - ceil(kappa*E_p*(1-ro)))
        I += (ceil(kappa*E_p*ro) - ceil(gammaA*I_p) - ceil(gammaI*I_p))
        P += (ceil(kappa*E_p*(1-ro)) - ceil(gammaA*P_p) - ceil(gammaI*P_p))
        H += (ceil(gammaA*I_p) + ceil(gammaA*P_p) - ceil(gammaR*H_p))
        R = ceil(gammaI*P_p) + ceil(gammaI*I_p) + ceil(gammaR*H_p)
        if(S<0):
            E+=S
            S=0
        if(E<0):
            I+=E
            E=0
        if(I<0):
            R+=I
            I=0
        if(P<0):
            R+=P
            P=0
        if(H<0):
            R+=H
            H=0
        if(R<0):
            S+=R
            R=0
        wykres.append(liczba_zakazonych)
        #koniec h(n)
        n=n+1


    from bokeh.plotting import figure, output_file, show
    from bokeh.embed import components
    from bokeh.resources import CDN

    # # prepare some data
    # x = [number, 2, 3, 4, 5]
    # y = [6, 7, 2, 4, 5]
    # x1 = [1,1,4,1,1]
    # y1 = [10*(i**2) for i in range(5)]
    

    # # create a new plot with a title and axis labels
    # p = figure(title="simple line example", x_axis_label='x', y_axis_label='y', x_range=[0,11], y_range=(0,11),width=400,height=200, sizing_mode='scale_width')

    # # # add a line renderer with legend and line thickness
    # # p.line(x, y, legend_label="Temp.", line_width=8)
    # # p.line(x1,y1, line_color="red")
    # # p.circle(x,x, fill_color="green", size=13)
    # # p.circle(x, x, legend_label="y=x", fill_color="white", size=8)
    # # p.line(x, y1, legend_label="y=10^x^2", line_color="orange", line_dash="4 4")
    # # script1, div1 = components(p)

    w=figure(title="Liczba zarażonych", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/2)),width=400,height=200, sizing_mode='scale_width')
    w.line(wykres_krokow, wykres, line_width=4)
    wjs, wd = components(w)

    p=figure(title="Wykres liczby superroznosicieli", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/200)),width=400,height=200, sizing_mode='scale_width')
    p.line(wykres_krokow, wykres_P, line_width=4)
    pjs, pd = components(p)

    s=figure(title="Wykres liczby osób podanych na zarażenie", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N)),width=400,height=200, sizing_mode='scale_width')
    s.line(wykres_krokow, wykres_S, line_width=4)
    sjs, sd = components(s)
    
    e=figure(title="Wykres liczby osób bezpośrednio narażonych na zachorowanie", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/2)),width=400,height=200, sizing_mode='scale_width')
    e.line(wykres_krokow, wykres_E, line_width=4)
    ejs, ed = components(e)

    i=figure(title="Wykres liczby osób zarażających", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/2)),width=400,height=200, sizing_mode='scale_width')
    i.line(wykres_krokow, wykres_I, line_width=4)
    ijs, id1 = components(i)

    r=figure(title="Wykres liczby osób ozdrowiałych", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/2)),width=400,height=200, sizing_mode='scale_width')
    r.line(wykres_krokow, wykres_R, line_width=4)
    rjs, rd = components(r)

    h=figure(title="Wykres liczby osób hospitalizowanych", x_axis_label='liczba kroków', y_axis_label='liczba osób', x_range=[0,ilosc_krokow], y_range=(0,(N/2)),width=400,height=200, sizing_mode='scale_width')
    h.line(wykres_krokow, wykres_H, line_width=4)
    hjs, hd = components(h)

    cdn_js=CDN.js_files[0]
    cdn_css = CDN.css_files
    


    return render_template("plot.html",
        cdn_js=cdn_js,
        cdn_css=cdn_css,
        wjs=wjs, wd=wd,
        pjs=pjs,pd=pd,
        sjs=sjs, sd=sd,
        ejs=ejs, ed=ed,
        ijs=ijs, id1=id1,
        rjs=rjs, rd=rd,
        hjs=hjs, hd=hd,
        zadany_procent_zarazonych=zadany_procent_zarazonych,
        liczba_osob=liczba_osob,
        liczba_zakazonych=stratowa_liczba_zakazonych,
        gammaA= gammaA,
        gammaI = gammaI,
        gammaR = gammaR,
        jota= jota,
        kappa=kappa,
        ro=ro,
        zadana_liczba_chorych=zadana_liczba_chorych,
        )


@app.route("/")
def home():
	return render_template("home.html")
	
@app.route('/bravo/')
def about():
    return render_template("about.html")


if __name__== "__main__":
    # import webbrowser
    # webbrowser.open("http://127.0.0.1:5000/")
    app.run(debug=True)