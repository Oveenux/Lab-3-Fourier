import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from scipy import signal as sg
from PIL import Image as img

st.title('Fourier')
tab1, tab2, tab3 = st.tabs(["Presentación", "Serie", "Transformada"])
a = 1
p = 1
n = 1
A1 = 1
A2 = 1
A3 = 1
ch = 1
sig = 0

fig , ax = plt.subplots(2,2)
fig1 , ax1 = plt.subplots(1,2)
fig0 , ax0 = plt.subplots()
op = st.sidebar.radio( 'Seleccione la operación de Fourier',
    ('Serie', 'Transformada'))
if op == 'Serie':
    st.sidebar.header('Serie')

    # inputs para serie
    sig = st.sidebar.selectbox( 'Señal x(t) :' , [ 'Exponencial', 'Sinusoidal rectificada' , 'Cuadrada' , 'Triangular' , 'Rampa Trapezoidal' ])

    st.sidebar.write('Inserte los siguientes parámetros:')

    a = st.sidebar.number_input('Inserte un valor de amplitud:', step = 1)
    p = st.sidebar.number_input('Inserte un valor para el periodo',value = 1.0 , step = 0.1)
    n = st.sidebar.number_input('Inserte un valor para el número de armónicos', min_value = 2, step = 1)
    if a == 0:
        a = 1
    
else:
    st.sidebar.header('Tranformada')
    st.sidebar.write('Inserte los siguientes parámetros:')

    N = st.sidebar.number_input('Inserte el numero de muestras:', step = 50,min_value = 1)
    Fs = st.sidebar.number_input('Inserte una frecuencia de muestreo', step = 100, min_value = 1)
    A1 = st.sidebar.number_input('Inserte un valor de amplitud 1:', step = 1)
    A2 = st.sidebar.number_input('Inserte un valor de amplitud 2:', step = 1)
    A3 = st.sidebar.number_input('Inserte un valor de amplitud 3:', step = 1)
    f1 = st.sidebar.number_input('Inserte el valor de frecuencia 1:', step = 1)
    f2 = st.sidebar.number_input('Inserte el valor de frecuencia 2:', step = 1)
    f3 = st.sidebar.number_input('Inserte el valor de frecuencia 3:', step = 1)
    
with tab1:
    col1, col2= st.columns(2)
    with col1: 
        logo1 = img.open('Logo.PNG')
        st.image(logo1, width = 300)
    with col2:
        logo2 = img.open('UN.PNG')
        st.image(logo2, width = 300)
    st.write('Hola, mediante esta app puedes seleccionar señales para reconstruir mediante series de Fourier. '
                'Adicionalmente también puede aplicar tranformada de Fourier a una determinada señal y conocer la frecuencia de operacion de esta misma. '
                ' Además de apreciar sus graficas y la grafica resultante de las operaciones.')
with tab2:
    
    st.subheader('Serie de Fourier')

    if op == 'Serie':

        dt = 0.01
        x = 0
        xf = 0
        wo = (2*np.pi)/p
        fs = np.pi/(20*wo)
        t = np.arange(0 , p+fs , fs)
        t1 = np.arange( 0 , 2*p , fs )

        # Generacion de señales

        if sig == 'Exponencial' :
            x = a*np.exp(t)
        elif sig == 'Sinusoidal rectificada' :
            x = abs(a*np.sin (wo*t))
        elif sig == 'Cuadrada' :
            x = a*sg.square(t*wo, 0.5)
        elif sig == 'Triangular' :
            x = a*sg.sawtooth(t*wo, 0.5)
        else :
            t = np.arange(0 , p,0.01)
            p1 = 0
            p2 = p
            pe = a / (p2-p*(2/3)-p1)
            x = np.piecewise(t, [(t <= p2-p*(2/3)) , (t>p2-p*(2/3)) & (t<p2-(p/3)) , (t>=p2-(p/3))],[lambda x: (t[t<=p2-p*(2/3)]-p1)*pe , lambda x: a , lambda x: (t[t>=p2-p*(1/3)]-p2)*(-pe)])


        # Series de Fourier

        ak = np.zeros(n)
        bk = np.zeros(n)
        m = len(t)

        nk = np.zeros(n)
        fk = np.zeros(n)
        nwo = np.arange(0,n,)
        A0 = 0
        for i in range ( 1 , m ):
            A0 = A0 + (1/p)*x[i]*dt
            
        for i in range ( 1 , n ):
            for j in range ( 1 , m ):
                ak[i] = ak[i] + ((2/p) * x[j] * np.cos( i*t[j]*wo ))*dt
                bk[i] = bk[i] + ((2/p) * x[j] * np.sin( i*t[j]*wo ))*dt
            
            nk[i] = ((ak[i]**2)*(bk[i]**2))**0.5
            fk[i] = -np.arctan(bk[i]/ak[i])
          

        sf = t1*0+A0

        for i in range (1 , n):
            sf = sf + ak[i] * np.cos(i*wo*t1) + bk[i] * np.sin(i*wo*t1)
            mxf = max(sf)/max(x)
            xf = sf/mxf

        ax[0,0].plot( t , x , color = 'k' )
        ax[0,0].set_ylabel('Amplitud')
        ax[0,0].set_xlabel('Tiempo')
        ax[0,0].set_title(f'{sig}')
        ax[0,1].plot( t1 , xf , color = 'b')
        ax[0,1].set_xlabel('Tiempo')
        ax[0,1].set_title('Representación por series de Fourier')
        ax[1,0].stem( nwo , nk , linefmt = 'darkslateblue' , markerfmt = ('darkslateblue','s'), basefmt = 'k')
        ax[1,0].set_xlabel('Espectro de amplitud')
        ax[1,1].set_xlabel('Espectro de fase')
        ax[1,1].stem( nwo , fk , linefmt = 'rebeccapurple', markerfmt = ('rebeccapurple','s') , basefmt = 'k')
        st.pyplot(fig)
with tab3:

    st.subheader('Transformada de Fourier')
    
    if op == 'Transformada':
        
        ch = 1+(max(f1,f2,f3))
        if (Fs < (2*ch)) :
            Fs = 4*ch
        if (N <ch) :
            st.warning('Recuerde considerar que pocas muestras reducen la precisión')
        
        T = 1/Fs
        v = 1/ch
        x = np.linspace(0.0, N*T , N)
        y = A1*np.sin(f1 *2*np.pi*x) + A2*np.sin(f2 *2*np.pi*x) + A3*np.sin(f3 *2*np.pi*x)
        func = f'x(t) = {A1}sen({f1}t) + {A2}sen({f2}t) + {A3}sen({f3}t)'
        st.latex (f' {func} ')
        v = 1/ch
        ax0.plot(x,y, color = 'darkcyan')
        plt.xlim(0,50*v)
        plt.ylabel('Amplitud')
        plt.xlabel('Tiempo')
        y_f = np.fft.fft(y)
        x_f = np.linspace(0.0, 1.0/(2.0*T), N//2)
        ax1[0].plot(x_f, 2.0/N * np.abs(y_f[:N//2]) , color = 'red')
        ax1[0].set_xlabel('Espectro en Amplitud')
        ax1[1].plot(x_f, 2.0/N * np.angle(y_f[:N//2]) , color = 'darkred')
        ax1[1].set_xlabel('Espectro en Fase')
        st.pyplot(fig0)
        st.pyplot(fig1)
        