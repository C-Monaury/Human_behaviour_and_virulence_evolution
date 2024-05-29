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
        \frac{d \upsilon}{\d \tau} &=  a \upsilon [(1 -x) b^{'}(\upsilon)(s + \sigma i) - 1] \,,\\
        \frac{d x}{d  \tau}   & =  a x (1 -x) (i - \kappa) \,.
   \end{array}
   \right.
             ''' 
)



st.header("Parameters")
col1, col2 = st.columns(2)

 st.write("Time")
    tmax = st.slider("Simulation time",1,1000)

with col1: 

    st.write("Shape parameters")
    sig =st.slider("Sampling rate",min_value = 0.0, max_value = 1000.0,step = 0.1)
    supinfec = st.slider("Superinfection",min_value = 0.0, max_value = 1.0,step = 0.01)
    
    A=st.slider("Variance in virulence",min_value = 0.1, max_value = 10.0,step = 0.1)
    


with col2: 
    st.write("Equilibrum parameters")
    pay = st.slider("Ratio of cooperator payments to non-cooperators payments",min_value = 0.0, max_value = 1.0,step = 0.01)
    
    c = st.slider("Constant infection capacity",min_value = 0.0, max_value = 10.0,step = 0.1)
    k= st.slider("Trade-off shape parameter",min_value = 0.10, max_value = 1.0,step = 0.0001)

pas = 0.01
nbr_pas = int(tmax/pas)

###########################################Choix du trade offs

plot1 = st.selectbox("Display trade-off graphic ? :chart_with_upwards_trend:",("Yes","No"))



trade_choix = st.selectbox("Choix du trade-off",["cx^k","(x*c)/(k+x)"])
if trade_choix == "cx^k":
    def beta(x,c, k ):
        return(c*x**k)
    def beta2(x,c, k ):
        return(c*k*x**(k-1))
    st.subheader("Forme du trade off")
    droite1 = np.zeros(200)
    droite2 = np.zeros(200)
    droite3 = np.zeros(200)
    for y in range(200):
        droite1[y] = y/10
        droite2[y] = 1+y/10
        droite3[y]=beta(y,c,k)
    figtrade, ax1 = plt.subplots()
    ax1.plot(range(200),droite1,"red")
    ax1.plot(range(200),droite2,"black")
    ax1.plot(range(200),droite3,"purple")
    ax1.set_xlabel('Clairance')
    ax1.set_ylabel('Transmission')
        

if trade_choix == "(x*c)/(k+x)":
    def beta(x,c, k ):
        return((x*c)/(k+x))
    def beta2(x,c, k ):
        return((c*k)/(k+x)**2)
    st.subheader("Forme du trade off")
    droite1 = np.zeros(20)
    droite2 = np.zeros(20)
    droite3 = np.zeros(20)
    for y in range(20):
        droite1[y] = y
        droite2[y] = y+1
        droite3[y]=beta(y,c,k)
    figtrade, ax1 = plt.subplots()
    ax1.plot(range(20),droite1,"red")
    ax1.plot(range(20),droite2,"black")
    ax1.plot(range(20),droite3,"purple")
    ax1.set_xlabel('Clairance')
    ax1.set_ylabel('Transmission')

if plot1=="Yes":
    st.pyplot(figtrade)


################################################# Coeur du modèle

#Ancienne version non fonctionelle



        



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

st.write("Valeurs initiales")
col221,col21,col22,col23 = st.columns(4)
with col221:
    s0 = st.slider("Sains initiale",min_value = 0.0,max_value = 1.0, step = 0.01)
with col21:
    i0 = st.slider("Infectés ",min_value = 0.0,max_value = 1.0, step = 0.01)
with col22:
    c0 = st.slider("Virulence initiale",0.0,10.0)
with col23:
    x0 = st.slider("Coopérateurs",min_value = 0.01,max_value = 1.00, step = 0.01)


###########################PLOT2D
st.subheader("Dynamiques des 3 compartiments en fonction du temps")

plot2 = st.selectbox("Display graphic ? :chart_with_upwards_trend:",("Yes","No"))

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


ax1.set_xlabel('Temps')
ax1.set_ylabel('Prévalence', color='red')
ax2.set_ylabel('Coopérateurs', color='black')






fignocoop, ax1 = plt.subplots()
ax1.plot(temps,sol_sanscoop[:,0],"green")
ax1.plot(temps,sol_sanscoop[:,1],"red")
ax1.plot(temps,sol_sanscoop[:,2],"purple")



ax1.set_xlabel('Temps')
ax1.set_ylabel('Prévalence', color='red')
ax1.set_ylabel('Virulence', color='purple')


if plot2 == "Yes":
    row1, row2  =st.columns(2)
   
    with row1:
        st.header("model coupling behavior and evolution")
        st.pyplot(fig1)

    with row2:
        st.header("model with only evolution ")
        st.pyplot(fignocoop)








############################## PLOT 

plot3 = st.selectbox("Display virulence evolution ? :chart_with_upwards_trend:",("Yes","No"))


st.subheader("Virulence evolution")


temps = np.linspace(0,tmax,nbr_pas)
pay1 = st.slider("Rapport du payement des coopérateurs sur celui des défecteurs1",min_value = 0.0, max_value = 1.0,step = 0.01)
pay2 = st.slider("Rapport du payement des coopérateurs sur celui des défecteurs2",min_value = 0.0, max_value = 1.0,step = 0.01)


sol_coop = odeint(model, y0 = [s0,i0 , c0,x0], t=temps,args = ( c, k,  A, supinfec,sig,pay1))
sol_coop2 = odeint(model, y0 = [s0,i0 , c0,x0], t=temps,args = ( c, k, A, supinfec,sig,pay2))

sol_sanscoop = odeint(model_sanscoop, y0 = [s0,i0 , c0], t=temps,args = ( c, k, A, supinfec,sig,pay))



temps = np.linspace(0,tmax,nbr_pas)

fig1, ax1 = plt.subplots()
ax1.plot(temps,sol_coop[:,2],color="blue")
ax1.plot(temps,sol_coop2[:,2],color = "blue",linestyle='dotted')
ax1.plot(temps,sol_sanscoop[:,2],"black")



ax1.set_xlabel('Temps')
ax1.set_ylabel('Virulence', color='red')
ax2.set_ylabel('Coopérateurs', color='black')

if plot3 == "Yes":
    st.pyplot(fig1)
