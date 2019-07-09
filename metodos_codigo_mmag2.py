import os
import matplotlib.pyplot as pyplot
from sympy import sympify

incHorizontal = 0.001

arq = open('saida.txt', 'w')


def calcPontospeloMetodo(YdosMets,y,x0,ordem,inc,expr, texto,metodo):

    t0 = float(x0)
    

    #Euler
    if(metodo == 1 or metodo == 6 or metodo == 11):
        if(metodo == 1):
            texto.append('Metodo Adan-Bashforth por Euler\n')
        
        elif(metodo == 6):
            texto.append('Metodo Adan-Multon por Euler\n')
        
        elif(metodo == 11):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler\n')

        texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )
        texto.append("h = " + str(inc) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):
            
            y = float(y + inc*expr.subs([("t", t0) , ("y", float(y))]))
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+inc
            YdosMets.append(y)
            cont = cont+1


    elif(metodo == 2 or metodo == 7 or metodo == 12):
        if(metodo == 2):
            texto.append('Metodo Adan-Bashforth por Euler Inverso\n')
            
        elif(metodo == 7):
            texto.append('Metodo Adan-Multon por Euler Inverso\n')
            
        elif(metodo == 12):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler Inverso\n')
        

        texto.append('y('+ str(x0) + ") = " + str(y)+'\n' )
        texto.append("h = " + str(inc) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):
            yporEuler = float(y + inc*expr.subs([("t", t0) , ("y", float(y))]))
            y = float(y + inc*expr.subs([("t", t0 + inc) , ("y", float(yporEuler))]))
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+inc
            YdosMets.append(y)
            cont = cont+1

    elif(metodo == 3  or metodo == 8 or metodo ==13):
        if(metodo == 3):
            texto.append('Metodo Adan-Bashforth por Euler Aprimorado\n')
            
        elif(metodo == 8):
            texto.append('Metodo Adan-Multon por Euler Aprimorado\n')
            
        elif(metodo == 13):
            texto.append('Metodo Formula Inversa de Diferenciacao por Euler Aprimorado\n')
            

        texto.append('y('+ str(x0) + ") = " + str(y)+"\n" )
        texto.append("h = " + str(inc) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        
        while (cont < ordem-1 ):
            yporEuler = float(y + inc*expr.subs([("t", t0) , ("y", float(y))]))            
            y = y + ( expr.subs([("t", t0 + inc) , ("y", float(yporEuler))]) +  expr.subs([("t", t0) , ("y", float(y))])) * inc *0.5
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+inc
            YdosMets.append(y)
            cont = cont+1

    elif(metodo == 4 or metodo == 9 or metodo ==14):
        if(metodo == 4):
            texto.append('Metodo Adan-Bashforth por Runge-Kuuta\n')
            
        elif(metodo == 9):
            texto.append('Metodo Adan-Multon por Runge-Kuuta\n')
            
        elif(metodo == 14):
            texto.append('Metodo Formula Inversa de Diferenciacao por Runge-Kutta\n')
            

        texto.append('y('+ str(x0) + ") = " + str(y) +"\n")
        texto.append("h = " + str(inc) + '\n')
        cont=1
        texto.append(str(0) + " "+ str(y) +'\n')
        while (cont < ordem-1 ):

            k1 = float(expr.subs([("t", t0) , ("y", float(y))]))
            k2 = float(expr.subs([("t", t0 + inc/2) , ("y", float(y + inc*k1/2))]))
            k3 = float(expr.subs([("t", t0 + inc/2) , ("y", float(y + inc*k2/2))]))
            k4 = float(expr.subs([("t", t0 + inc) , ("y", float(y + inc*k3))]))
            
            y = y + ((k1+2*k2+2*k3+k4)*inc)/6
            texto.append(str(cont) + " "+ str(y) +'\n')
            t0 = t0+inc
            YdosMets.append(y)
            cont = cont+1

def metAdamMulton(y,x0,inc, qtdPassos, funcaoDerivada, ordem, metodo,id):
    pyplot.clf()
    expr = sympify(funcaoDerivada)
    pyplot.title('Adam Multon')
    texto = []

    t0 = float(x0)
    
    YdosMets = [0.0]
    YdosMets.pop(0)
    YdosMets.append(y)
    if(metodo==5):
        texto.append('Metodo Adan-Multon\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(inc) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YdosMets = y
    
    
    calcPontospeloMetodo(YdosMets,y,x0,ordem,inc,expr, texto,metodo)
    #coloca ponto n+1
    YdosMets.append(0.0)
        
    t0 = float(x0)
    t0 = t0+(ordem-1)*inc
    
    contador = ordem-1

    
    while(contador<= qtdPassos):
        #Euler
        if(metodo == 5 or metodo == 6):
            YdosMets[-1] = YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc
        
        #Euler inverso
        if(metodo ==7 ):
            #usando metodo de Euler
            yn1 = YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc

            #usando metodo de Euler inverso  
            YdosMets[-1] = YdosMets[-2] + expr.subs([("t", t0) , ("y", float(yn1))])*inc
        
        #Euler aprimorado
        if(metodo ==8 or metodo == 9) :
            #usando metodo de Euler
            yn1 =YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc

            #usando metodo de Euler aprimorado  
            YdosMets[-1] = YdosMets[-2] + ( expr.subs([("t", t0) , ("y", float(yn1))]) +  expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))]))*inc*0.5    
            
        if(ordem == 2):
            parte1 = ( (1.0/2.0) * float( expr.subs([("t", t0) , ("y", YdosMets[1])])    ) )
            parte2 =  ( (1.0/2.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 + parte2 )   
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 3):
            parte1 = ( (5.0/12.0) * float( expr.subs([("t", t0) , ("y", YdosMets[2])])    ) )
            parte2 =  ( (2.0/3.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[1])])) )
            parte3 = ( (-1.0/12.0) * float(expr.subs([("t", t0-2*inc) , ("y", YdosMets[0])])) ) 
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2 + parte3)   
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 4):
            parte1 = ( (3.0/8.0) * float( expr.subs([("t", t0) , ("y", YdosMets[3])])    ) )
            parte2 =  ( (19.0/24.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[2])])) )
            parte3 = ( (-5.0/24.0) * float(expr.subs([("t", t0-2*inc) , ("y", YdosMets[1])])) ) 
            parte4 = ( (1.0/24.0) * float(expr.subs([("t", t0-3*inc) , ("y", YdosMets[0])])) ) 
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 +parte4)  
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)

        elif(ordem == 5):
            parte1 = ( (251.0/720.0) * float( expr.subs([("t", t0) , ("y", YdosMets[4])])    ) )
            parte2 =  ( (323.0/360.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[3])])) )
            parte3 = ( (-11.0/30.0) * float(expr.subs([("t", t0-2*inc) , ("y", YdosMets[2])])) ) 
            parte4 = ( (53.0/360.0) * float(expr.subs([("t", t0-3*inc) , ("y", YdosMets[1])])) ) 
            parte5 = ( (-19.0/720.0) * float(expr.subs([("t", t0 - 4*inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4 + parte5 ) 
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 6):
            parte1 = ( (95.0/288.0) * float( expr.subs([("t", t0) , ("y", YdosMets[5])])    ) )
            parte2 =  ( (1427.0/1440.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[4])])) )
            parte3 = ( (-133.0/240.0) * float(expr.subs([("t", t0 - 2*inc) , ("y", YdosMets[3])])) ) 
            parte4 = ( (241.0/720.0) * float(expr.subs([("t", t0 - 3*inc) , ("y", YdosMets[2])])) ) 
            parte5 = ( (-173.0/1440.0) * float(expr.subs([("t", t0 - 4*inc) , ("y", YdosMets[1])])) )
            parte6 = ( (3.0/160.0) * float(expr.subs([("t", t0 - 5*inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6)
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
            
        elif(ordem == 7):
            parte1 = ( (19087.0/60480.0) * float( expr.subs([("t", t0) , ("y", YdosMets[6])])    ) )
            parte2 =  ( (2713.0/2520.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[5])])) )
            parte3 = ( (-15487.0/20160.0) * float(expr.subs([("t", t0 - 2*inc) , ("y", YdosMets[4])])) ) 
            parte4 = ( (586.0/945.0) * float(expr.subs([("t", t0 - 3*inc) , ("y", YdosMets[3])])) ) 
            parte5 = ( (-6737.0/20160.0) * float(expr.subs([("t", t0 - 4*inc) , ("y", YdosMets[2])])) )
            parte6 = ( (263.0/2520.0) * float(expr.subs([("t", t0 - 5*inc) , ("y", YdosMets[1])])) )
            parte7 = ( (-863.0/60480.0) * float(expr.subs([("t", t0 - 6*inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6 + parte7)
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
            
        texto.append(str(contador) + " "+ str(YdosMets[-2]) +'\n')
        contador = contador +1
        t0 = t0+ inc

        pyplot.savefig("graficos/" + str(id) + ".png")
  
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metInverso(y,x0,inc, qtdPassos, funcaoDerivada, ordem, metodo,id):
    
    expr = sympify(funcaoDerivada)
    texto = []

    t0 = float(x0)
    
    YdosMets = [0.0]
    YdosMets.pop(0)
    YdosMets.append(y)
    if(metodo==10):
        texto.append('Metodo Formula Inversa de Diferenciacao\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(inc) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YdosMets = y
    
    
    calcPontospeloMetodo(YdosMets,y,x0,ordem,inc,expr, texto,metodo)
    #coloca ponto n+1
    YdosMets.append(0.0)
        
    t0 = float(x0)
    t0 = t0+(ordem-1)*inc
    
    contador = ordem-1

    
    while(contador<= qtdPassos):
        #Euler
        if(metodo == 10 or metodo == 11):
            YdosMets[-1] = YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc
        
        #Euler inverso
        if(metodo ==12 ):
            #usando metodo de Euler
            yn1 =YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc

            #usando metodo de Euler inverso  
            YdosMets[-1] = YdosMets[-2] + expr.subs([("t", t0) , ("y", float(yn1))])*inc
            
        #Euler aprimorado
        if(metodo ==13 or metodo == 14) :
            #usando metodo de Euler
            yn1 =YdosMets[-2] + expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))])*inc

            #usando metodo de Euler Aprimorado
            YdosMets[-1] = YdosMets[-2] + ( expr.subs([("t", t0) , ("y", float(yn1))]) +  expr.subs([("t", t0-inc) , ("y", float(YdosMets[-2]))]))*inc*0.5    
                    
        if(ordem == 2):
            parte1 = ( (1.0/1.0) * float( expr.subs([("t", t0) , ("y", YdosMets[1])])    ) )*inc
            parte2 = ( (1.0/1.0) * YdosMets[0])
            novoValor =  ( parte1 + parte2 )   
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 3):
            parte1 = ( (2.0/3.0) * float( expr.subs([("t", t0) , ("y", YdosMets[2])])    ) )*inc
            parte2 =  ( (4.0/3.0) * YdosMets[1])
            parte3 = ( (-1.0/3.0) * YdosMets[0])
            novoValor = ( parte1 +parte2 + parte3)   
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)

        elif(ordem == 4):
            parte1 = ( (6.0/11.0) * float( expr.subs([("t", t0) , ("y", YdosMets[3])])    ) )*inc
            parte2 =  ( (18.0/11.0) * YdosMets[2])
            parte3 = ( (-9.0/11.0) * YdosMets[1])
            parte4 = ( (2.0/11.0) * YdosMets[0])
            novoValor = ( parte1 +parte2  + parte3 +parte4)  
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)

        elif(ordem == 5):
            parte1 = ( (12.0/25.0) * float( expr.subs([("t", t0) , ("y", YdosMets[4])])    ) )*inc
            parte2 = ( (48.0/25.0) * YdosMets[3])
            parte3 = ( (-36.0/25.0) * YdosMets[2])
            parte4 = ( (16.0/25.0) * YdosMets[1])
            parte5 = ( (-3.0/25.0) * YdosMets[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 ) 
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 6):
            parte1 = ( (60.0/137.0) * float( expr.subs([("t", t0) , ("y", YdosMets[5])])    ) )*inc
            parte2 =  ( (300.0/137.0) * YdosMets[4])
            parte3 = ( (-300.0/137.0) * YdosMets[3])
            parte4 = ( (200.0/137.0) * YdosMets[2])
            parte5 = ( (-75.0/137.0) * YdosMets[1])
            parte6 = ( (12.0/137.0) * YdosMets[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6)
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
        
        elif(ordem == 7):
            parte1 = ( (60.0/147.0) * float( expr.subs([("t", t0) , ("y", YdosMets[6])])    ) )*inc
            parte2 = ( (360.0/147.0) * YdosMets[5])
            parte3 = ( (-450.0/147.0) * YdosMets[4])
            parte4 = ( (400.0/147.0) * YdosMets[3])
            parte5 = ( (-225.0/147.0) * YdosMets[2])
            parte6 = ( (72.0/147.0) * YdosMets[1])
            parte7 = ( (-10.0/147.0) * YdosMets[0])
            novoValor =  ( parte1 +parte2  + parte3 + parte4 + parte5 + parte6 + parte7)
            YdosMets.pop(0)
            YdosMets.pop(-1)
            YdosMets.append(novoValor)
            YdosMets.append(0.0)
            
        texto.append(str(contador) + " "+ str(YdosMets[-2]) +'\n')
        contador = contador +1
        t0 = t0+ inc


    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################




def metAdamBashforth(y,x0,inc, qtdPassos, funcaoDerivada, ordem, metodo,id):
    expr = sympify(funcaoDerivada)
    texto = []

    t0 = float(x0)
    
    YdosMets = [0.0]
    YdosMets.pop(0)
    YdosMets.append(y)
    if(metodo==0):
        texto.append('Metodo Adan-Bashforth\n')
        texto.append('y('+ str(x0) + ") = " + str(y[0])+"\n" )
        texto.append("h = " + str(inc) + '\n')

        cont=0
        while(cont<len(y)):
            a = y[cont]            
            texto.append(str(cont) + " "+ str(a) +'\n')
            cont = cont+1
        YdosMets = y
    
    
    calcPontospeloMetodo(YdosMets,y,x0,ordem,inc,expr, texto,metodo)
        
    t0 = float(x0)
    t0 = t0+(ordem-2)*inc
    contador = ordem-1
    while(contador<= qtdPassos):
        

        if(ordem == 2):
            parte1 = ( (1.0/1.0) * float( expr.subs([("t", t0) , ("y", YdosMets[0])])    ) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 )   
            YdosMets.pop(0)
            YdosMets.append(novoValor)
        
        elif(ordem == 3):
            parte1 = ( (3.0/2.0) * float( expr.subs([("t", t0) , ("y", YdosMets[1])])    ) )
            parte2 =  ( (-1.0/2.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2 )   
            YdosMets.pop(0)
            YdosMets.append(novoValor)
        
        elif(ordem == 4):
            parte1 = ( (23.0/12.0) * float( expr.subs([("t", t0) , ("y", YdosMets[2])])    ) )
            parte2 =  ( (-4.0/3.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[1])])) )
            parte3 = ( (5.0/12.0) * float(expr.subs([("t", t0-2*inc) , ("y", YdosMets[0])])) ) 
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3)  
            YdosMets.pop(0)
            YdosMets.append(novoValor)
        
        elif(ordem == 5):
            parte1 = ( (55.0/24.0) * float( expr.subs([("t", t0) , ("y", YdosMets[3])])    ) )
            parte2 =  ( (-59.0/24.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[2])])) )
            parte3 = ( (37.0/24.0) * float(expr.subs([("t", t0-2*inc) , ("y", YdosMets[1])])) ) 
            parte4 = ( (-3.0/8.0) * float(expr.subs([("t", t0-3*inc) , ("y", YdosMets[0])])) ) 
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4) 
            YdosMets.pop(0)
            YdosMets.append(novoValor)
        
        elif(ordem == 6):
            parte1 = ( (1901.0/720.0) * float( expr.subs([("t", t0) , ("y", YdosMets[4])])    ) )
            parte2 =  ( (-1387.0/360.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[3])])) )
            parte3 = ( (109.0/30.0) * float(expr.subs([("t", t0 - 2*inc) , ("y", YdosMets[2])])) ) 
            parte4 = ( (-637.0/360.0) * float(expr.subs([("t", t0 - 3*inc) , ("y", YdosMets[1])])) ) 
            parte5 = ( (251.0/720.0) * float(expr.subs([("t", t0 - 4*inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4 + parte5)
            YdosMets.pop(0)
            YdosMets.append(novoValor)
            
        elif(ordem == 7):
            parte1 = ( (4277.0/1440.0) * float( expr.subs([("t", t0) , ("y", YdosMets[5])])    ) )
            parte2 =  ( (-2641.0/480.0) * float(expr.subs([("t", t0-inc) , ("y", YdosMets[4])])) )
            parte3 = ( (4991.0/720.0) * float(expr.subs([("t", t0 - 2*inc) , ("y", YdosMets[3])])) ) 
            parte4 = ( (-3649.0/720.0) * float(expr.subs([("t", t0 - 3*inc) , ("y", YdosMets[2])])) ) 
            parte5 = ( (959.0/480.0) * float(expr.subs([("t", t0 - 4*inc) , ("y", YdosMets[1])])) )
            parte6 = ( (-95.0/288.0) * float(expr.subs([("t", t0 - 5*inc) , ("y", YdosMets[0])])) )
            novoValor = YdosMets[ordem-2] + inc * ( parte1 +parte2  + parte3 + parte4 + parte5)
            YdosMets.pop(0)
            YdosMets.append(novoValor)
            
        texto.append(str(contador) + " "+ str(YdosMets[ordem-2]) +'\n')
        contador = contador +1
        t0 = t0+ inc 

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################



def metRungeKutta(y0,x0,inc, qtdPassos, funcaoDerivada,id):
    expr = sympify(funcaoDerivada)
    texto = []
    texto.append('Metodo de Runge-Kuuta\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(inc) + '\n')

    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        

        k1 = float(expr.subs([("t", t0) , ("y", float(y0))]))
        k2 = float(expr.subs([("t", t0 + inc/2) , ("y", float(y0 + inc*k1/2))]))
        k3 = float(expr.subs([("t", t0 + inc/2) , ("y", float(y0 + inc*k2/2))]))
        k4 = float(expr.subs([("t", t0 + inc) , ("y", float(y0 + inc*k3))]))
        
        texto.append(str(contador) + " "+ str(y0) +'\n')

        y0 = y0 + ((k1+2*k2+2*k3+k4)*inc)/6

        contador = contador +1
        t0 = t0+ inc

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metEuler(y0,x0,inc, qtdPassos, funcaoDerivada,id):
    
    expr = sympify(funcaoDerivada)

    texto = []
    texto.append('Metodo de Euler\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(inc) + '\n')


    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        texto.append(str(contador) + " "+ str(y0) +'\n')
        y0 = y0 + expr.subs([("t", t0) , ("y", float(y0))])*inc
        contador = contador +1
        t0 = t0+ inc

    #escreve no arq
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metEulerInverso(y0,x0,inc, qtdPassos, funcaoDerivada,id):
    expr = sympify(funcaoDerivada)
    texto = []
    texto.append('Metodo de Euler Inverso\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(inc) + '\n')

    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):

        texto.append(str(contador) + " "+ str(y0) +'\n')

        #apilando metodo de Euler
        yn1 =y0 + expr.subs([("t", t0) , ("y", float(y0))])*inc

        #usando metodo de Euler inverso  
        y0 = y0 + expr.subs([("t", t0+inc) , ("y", float(yn1))])*inc
        
        contador = contador +1
        t0 = t0 + inc

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def metEulerAprimorado(y0,x0,inc, qtdPassos, funcaoDerivada,id):
    expr = sympify(funcaoDerivada)
    texto = []
    texto.append('Metodo de Euler Aprimorado\n')
    texto.append('y('+ str(x0) + ') = ' + str(y0)+ '\n')
    texto.append("h = " + str(inc) + '\n')


    t0 = float(x0)

    contador =0
    while(contador<=qtdPassos):
        
        texto.append(str(contador) + " "+ str(y0) +'\n')
        #usando metodo de Euler
        yn1 = y0 + expr.subs([("t", t0) , ("y", float(y0))])*inc

        #usando metodo de Euler inverso  
        y0 = y0 + ( expr.subs([("t", t0 + inc) , ("y", float(yn1))]) +  expr.subs([("t", t0) , ("y", float(y0))]))*inc*0.5

        contador = contador + 1
        t0 = t0 + inc

    #escreve no arq
    
    arq.writelines(texto)
    arq.write('\n')
    
    #######################


def main():

    numeroLinha = 0
    
    try:
        with open('entrada.txt') as arq:
            for linha in arq:
                numeroLinha = numeroLinha +1
                entrada = linha.split()
                if(entrada[0] == 'euler'):
                    metEuler(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'euler_inverso'):
                    metEulerInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'euler_aprimorado'):
                    metEulerAprimorado(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'runge_kutta'):
                    metRungeKutta(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],numeroLinha)
                elif(entrada[0] == 'adam_bashforth'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metAdamBashforth(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),0,numeroLinha)

                elif(entrada[0] == 'adam_bashforth_by_euler'):
                    metAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),1,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_euler_inverso'):
                    metAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),2,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_euler_aprimorado'):
                    metAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),3,numeroLinha)
                elif(entrada[0] == 'adam_bashforth_by_runge_kutta'):
                    metAdamBashforth(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),4,numeroLinha)
                elif(entrada[0] == 'adam_multon'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metAdamMulton(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),5,numeroLinha)

                elif(entrada[0] == 'adam_multon_by_euler'):
                    metAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),6,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_euler_inverso'):
                    metAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),7,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_euler_aprimorado'):
                    metAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),8,numeroLinha)
                elif(entrada[0] == 'adam_multon_by_runge_kutta'):
                    metAdamMulton(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),9,numeroLinha)
                
                elif(entrada[0] == 'formula_inversa'): 
                    listaY = [0.0]
                    listaY.pop(0)
                    cont=1
                    while cont < int(entrada[-1]):
                        listaY.append(float(entrada[cont]))
                        cont = cont+1
                    metInverso(listaY, float(entrada[-5]), float(entrada[-4]), int(entrada[-3]), entrada[-2],int(entrada[-1]),10,numeroLinha)

                elif(entrada[0] == 'formula_inversa_by_euler'):
                    metInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),11,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_euler_inverso'):
                    metInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),12,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_euler_aprimorado'):
                    metInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),13,numeroLinha)
                elif(entrada[0] == 'formula_inversa_by_runge_kutta'):
                    metInverso(float(entrada[1]), float(entrada[2]), float(entrada[3]), int(entrada[4]), entrada[5],int(entrada[6]),14,numeroLinha)
                else:
                    print ('O metodo ' + str(entrada[0] + 'nao existe nesse programa') ) 
                print ('Linha ' + str(linha) + ' terminou')
        arq.close()
    except IOError:
       print ('Arquivo entrada.txt nao existe')
    
main()
