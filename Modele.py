#application streamlit
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import seaborn as sns
import streamlit as st
from mpl_toolkits.mplot3d import Axes3D


st.title("How does human behavior influence the evolution of virulence?")
st.write("This model was produced during an M2 internship by Clément MONAURY and supervised by Frédéric Hamelin.")
st.header("System overview:")

st.latex(r'''
    \left  \{
    \begin{array}{r c l}
        \frac{d s}{d \tau}  & =  1 - (1 -x)b(\upsilon) s  i -  s \,,\\
        \frac{ d i}{d \tau}   & =  (1 -x)b(\upsilon) s i   - (\upsilon +1 ) i \,,\\
        \frac{d \upsilon}{d \tau} &=  a \upsilon [(1 -x) b^{'}(\upsilon)(s + \sigma i) - 1] \,,\\
        \frac{d x}{d  \tau}   & =  c x (1 -x) (i - \kappa) \,.
   \end{array}
   \right.
             ''' 
)

user = st.selectbox(
    "Parameters modification ",
    ("Simple" , "Acess to all parameters"))

st.header("Parameters")

st.write("Time")
if user=="Simple":
    tmax=100
else:
    tmax = st.slider("Simulation time",1,1000)



if user=="Acess to all parameters":
    col1, col2 = st.columns(2)
    with col1: 

        st.write("Shape parameters")
        sig =st.slider("Sampling rate : c",min_value = 0.0, max_value = 10.0,step = 0.1)
        supinfec = st.slider("Superinfection : sigma",min_value = 0.0, max_value = 1.0,step = 0.01)
        
        A=st.slider("Variance in virulence a",min_value = 0.1, max_value = 10.0,step = 0.1)
        


    with col2: 
        st.write("Equilibrum parameters")
        pay = st.slider("Ratio of cooperators payments to non-cooperators payments : kappa",min_value = 0.0, max_value = 1.0,step = 0.01)
        
        c = st.slider("Trade-off : constant infection capacity ",min_value = 0.0, max_value = 10.0,step = 0.1)
        k= st.slider("Trade-off shape parameter",min_value = 0.10, max_value = 1.0,step = 0.0001)






if user=="Simple":
    col1, col2 = st.columns(2)
    with col1: 
        st.write("Shape parameters")
        vitesse=st.selectbox("Ratio of virulence to learning rate of human",("Virulence evolution 10 times faster than human behaviour","Same speed","Human behaviour 10 times faster than pathogen"))
        if vitesse =="Virulence evolution 10 times faster than human behaviour":
            A = 10
            sig = 1
        if vitesse =="Same speed":
            A = 5
            sig = 5  
        if vitesse =="Human behaviour 10 times faster than pathogen":
            A = 1
            sig = 10   

        supinfec1 = st.number_input("Superinfection probability ",min_value = 0.0, max_value = 100.0,step = 1.0,value=10.0)
        supinfec = supinfec1/100
        

        


    with col2: 
        st.write("Equilibrum parameters")
        
        payement=st.selectbox("Ratio of cooperators payments to non-cooperators payments",("Low payout ratio (cooperators equilibrium)","Strong payout ratio"))
        if payement =="Low payout ratio (cooperators equilibrium)":
            pay = 0.05
        if payement =="Strong payout ratio":
            pay = 0.8        
        
        
        
        
        tradouf= st.selectbox("Choose a trade-off",("Strong (endemic equilibrium)","Medium (Disease free equilibrium with virulence)","Weak (Disease free equilibrium without virulence)"))
        if tradouf =="Strong (endemic equilibrium)":
            c=3.5
            k=0.5
        if tradouf =="Medium (Disease free equilibrium with virulence)":
            c=1.5
            k=0.7
        if tradouf =="Weak (Disease free equilibrium without virulence)":
            c=.5
            k=0.9
        
 



pas = 0.01
nbr_pas = int(tmax/pas)

###########################################Choix du trade offs

plot1 = st.selectbox("Display trade-off graphic ? :chart_with_upwards_trend:",("Yes","No"))


if user=="Simple":
    def beta(x,c, k ):
        return(c*x**k)
    def beta2(x,c, k ):
        return(c*k*x**(k-1))
   
    droite1 = np.zeros(2000)
    droite2 = np.zeros(2000)
    droite3 = np.zeros(2000)
    for y in range(2000):
        droite1[y] = y/100
        droite2[y] = 1+y/100
        droite3[y]=beta(y/100,c,k)
    figtrade, ax1 = plt.subplots()
    ax1.plot(droite1,droite1,"red")
    ax1.plot(droite1,droite2,"black")
    ax1.plot(droite1,droite3,"purple")
    ax1.set_xlabel('Virulence')
    ax1.set_ylabel('Transmission')
            

else:
    trade_choix = st.selectbox("Choose a trade-off relationship",["cx^k","(x*c)/(k+x)"])
    if trade_choix == "cx^k":
        def beta(x,c, k ):
            return(c*x**k)
        def beta2(x,c, k ):
            return(c*k*x**(k-1))

        droite1 = np.zeros(2000)
        droite2 = np.zeros(2000)
        droite3 = np.zeros(2000)
        for y in range(2000):
            droite1[y] = y/100
            droite2[y] = 1+y/100
            droite3[y]=beta(y/100,c,k)
        figtrade, ax1 = plt.subplots()
        ax1.plot(droite1,droite1,"red")
        ax1.plot(droite1,droite2,"black")
        ax1.plot(droite1,droite3,"purple")
        ax1.set_xlabel('Virulence')
        ax1.set_ylabel('Transmission')
            

    if trade_choix == "(x*c)/(k+x)":
        def beta(x,c, k ):
            return((x*c)/(k+x))
        def beta2(x,c, k ):
            return((c*k)/(k+x)**2)

        droite1 = np.zeros(2000)
        droite2 = np.zeros(2000)
        droite3 = np.zeros(2000)
        for y in range(2000):
            droite1[y] = y/100
            droite2[y] = y/100 +1
            droite3[y]=beta(y/100,c,k)
        figtrade, ax1 = plt.subplots()
        ax1.plot(droite1,droite1,"red")
        ax1.plot(droite1,droite2,"black")
        ax1.plot(droite1,droite3,"purple")
        ax1.set_xlabel('Virulence')
        ax1.set_ylabel('Transmission')        
        
        
        

if plot1=="Yes":
    st.subheader("Trade-off:")
    st.pyplot(figtrade)


################################################# Coeur du modèle


def model(Y0, t , c, k, A, supinfec,sig,pay) :
    #S, I , alpha, x = Y0
    S =Y0[0]
    I =Y0[1]
    alpha =Y0[2]
    x =Y0[3]
    N = S + I

    dS = 1 - (1 - x) * beta(alpha ,c , k) * I * S - S
    dI = (1 - x) * beta(alpha ,c , k) * I * S - alpha * I - I
    dalpha = A  *alpha*((1-x)*beta2(alpha,c, k) * (S + I* supinfec) - 1 )  
    dx =  sig*x * (1-x)*(I - pay)
    return(dS,dI,dalpha,dx)



def model_sanscoop(Y0, t , c, k,  A, supinfec,sig,pay) :
    #S, I , alpha, x = Y0
    S =Y0[0]
    I =Y0[1]
    alpha =Y0[2]

    N = S + I

    dS = 1 - beta(alpha ,c , k) * I * S - 1 * S
    dI =  beta(alpha ,c , k) * I * S - alpha * I - I 
    dalpha = A  *alpha*(beta2(alpha,c, k) * (S + I* supinfec) - 1 )  
    return(dS,dI,dalpha)








#valeurs de départ
if user=="Simple":
    s0,i0,c0,x0=0.9,0.1,0.5,0.5
    
else:    
    st.write("Initial values")
    col221,col21,col22,col23 = st.columns(4)
    with col221:
        s0 = st.slider("Suceptible",min_value = 0.01,max_value = 1.0, step = 0.01)
    with col21:
        i0 = st.slider("Infected ",min_value = 0.01,max_value = 1.0, step = 0.01)
    with col22:
        c0 = st.slider("Virulence",0.01,10.0)
    with col23:
        x0 = st.slider("Coopertators",min_value = 0.01,max_value = 1.00, step = 0.01)


###########################PLOT2D
st.subheader("Dynamiques des 3 compartiments en fonction du temps")

plot2 = st.selectbox("Display times plots graphics ? :chart_with_upwards_trend:",("Yes","No"))

temps = np.linspace(0,tmax,nbr_pas)
sol = odeint(model, y0 = [s0,i0 , c0,x0], t=temps,args = ( c, k,  A, supinfec,sig,pay))


sol_sanscoop = odeint(model_sanscoop, y0 = [s0,i0 , c0], t=temps,args = (c, k, A, supinfec,sig,pay))



temps = np.linspace(0,tmax,nbr_pas)

fig1, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(temps,sol[:,0],"green")
ax1.plot(temps,sol[:,1],"red")
ax1.plot(temps,sol[:,2],"purple")
ax2.plot(temps,sol[:,3],"black")


ax1.set_xlabel('Times')
ax1.set_ylabel('Density and virulence')
ax2.set_ylabel('Cooperators', color='black')






fignocoop, ax1 = plt.subplots()
ax1.plot(temps,sol_sanscoop[:,0],"green")
ax1.plot(temps,sol_sanscoop[:,1],"red")
ax1.plot(temps,sol_sanscoop[:,2],"purple")




ax1.set_xlabel('Times')
ax1.set_ylabel('Density and virulence')
ax2.set_ylabel('Cooperators', color='black')



if plot2 == "Yes":

    st.subheader("model coupling behavior and evolution")
    st.pyplot(fig1)


    st.subheader("model with only evolution ")
    st.pyplot(fignocoop)








############################## PLOT 



st.subheader("Virulence evolution")
plot3 = st.selectbox("Display virulence evolution ? :chart_with_upwards_trend:",("Yes","No"))


temps = np.linspace(0,tmax,nbr_pas)
pay1 = st.slider("Ratio of cooperators payments to non-cooperators payments 1",min_value = 0.0, max_value = 1.0,step = 0.01)
pay2 = st.slider("Ratio of cooperators payments to non-cooperators payments 2",min_value = 0.0, max_value = 1.0,step = 0.01)


sol_coop = odeint(model, y0 = [s0,i0 , c0,x0], t=temps,args = ( c, k,  A, supinfec,sig,pay1))
sol_coop2 = odeint(model, y0 = [s0,i0 , c0,x0], t=temps,args = ( c, k, A, supinfec,sig,pay2))

sol_sanscoop = odeint(model_sanscoop, y0 = [s0,i0 , c0], t=temps,args = ( c, k, A, supinfec,sig,pay))



temps = np.linspace(0,tmax,nbr_pas)

fig1, ax1 = plt.subplots()
ax1.plot(temps,sol_coop[:,2],color="blue")
ax1.plot(temps,sol_coop2[:,2],color = "blue",linestyle='dotted')
ax1.plot(temps,sol_sanscoop[:,2],"black")



ax1.set_xlabel('Times')
ax1.set_ylabel('Virulence', color='red')
ax2.set_ylabel('Coopérateurs', color='black')

if plot3 == "Yes":
    st.pyplot(fig1)
