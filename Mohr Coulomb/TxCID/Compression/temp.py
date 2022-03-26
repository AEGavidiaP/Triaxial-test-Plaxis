import numpy as np

i = 0;
Eps1_o = [1,2,3,4,5];
consolidation = [50,100];

tal = float(round((Eps1_o)*1.1,0))

if i==0:
   Dedo = np.zeros([tal, len(consolidation)])
   print(Dedo)

if len(Eps1_o)<len(Dedo[:,0]):
    valor = len(Dedo[:,0]) - len(Eps1_o)
    for j in range(valor):
        Eps1_o.append(Eps1_o[-1])
        j += 1

Dedo[:,i] =  Eps1_o
print(" ")
print(Dedo)



