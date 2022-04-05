import numpy as np
import matplotlib.pyplot as plt

def Triaxial_0(q, e1, pE, eV, e, q_data, e1_data, pE_data, eV_data, e_data, size_titulo, size_ejes, size_leyenda):
    
    #e1-q
    plt.subplot(221)
    l1 = plt.plot(e1, q, "b")
    l2 = plt.plot(e1_data,q_data, "k")
    plt.title('Stress-Strain', fontsize = size_titulo)
    plt.xlabel('$\epsilon_1$', fontsize = size_ejes)
    plt.ylabel('q', fontsize = size_ejes)
    plt.xlim(-0.2,0)
    plt.grid(True)
    
    #p-q
    plt.subplot(222)
    plt.plot(pE, q, "b")
    plt.plot(pE_data,q_data, "k")
    plt.title('p\'-q', fontsize = size_titulo)
    plt.xlabel('p\'', fontsize = size_ejes)
    plt.ylabel('q', fontsize = size_ejes)
    plt.grid(True)

    #e1-ev
    plt.subplot(223)
    plt.plot(e1, eV, "b")
    plt.plot(e1_data,eV_data, "k")
    plt.title('Axial-volumetric strain', fontsize = size_titulo)
    plt.xlabel('$\epsilon_1$', fontsize = size_ejes)
    plt.ylabel('$\epsilon_V$', fontsize = size_ejes)
    plt.xlim(-0.2,0)
    plt.grid(True)
    
    #p-e
    plt.subplot(224)
    plt.plot(pE, e, "b")
    plt.plot(pE_data,e_data, "k")
    plt.title('p\'-e', fontsize = size_titulo)
    plt.xlabel('p\'', fontsize = size_ejes)
    plt.ylabel('e', fontsize = size_ejes)
    plt.grid(True)
    
    return l1, l2

    
def Triaxial_1(q, e1, pE, eV, e, q_data, e1_data, pE_data, eV_data, e_data, size_titulo, size_ejes, size_leyenda):
    
    #e1-q
    plt.subplot(221)
    l3 = plt.plot(e1, q, "r")
    l4 = plt.plot(e1_data,q_data, "brown")
    
    #p-q
    plt.subplot(222)
    plt.plot(pE, q, "r")
    plt.plot(pE_data,q_data, "brown")

    #e1-ev
    plt.subplot(223)
    plt.plot(e1, eV, "r")
    plt.plot(e1_data,eV_data, "brown")
    
    #p-e
    plt.subplot(224)
    plt.plot(pE, e, "r")
    plt.plot(pE_data,e_data, "brown")
    
    return l3, l4


def Triaxial_2(q, e1, pE, eV, e, q_data, e1_data, pE_data, eV_data, e_data, size_titulo, size_ejes, size_leyenda):
    
    #e1-q
    plt.subplot(221)
    l5 = plt.plot(e1, q, "m")
    l6 = plt.plot(e1_data,q_data, "grey")
    
    #p-q
    plt.subplot(222)
    plt.plot(pE, q, "m")
    plt.plot(pE_data,q_data, "grey")

    #e1-ev
    plt.subplot(223)
    plt.plot(e1, eV, "m")
    plt.plot(e1_data,eV_data, "grey")
    
    #p-e
    plt.subplot(224)
    plt.plot(pE, e, "m")
    plt.plot(pE_data,e_data, "grey")
    
    return l5, l6

def graphics_triaxial(q, Eps1, pE, Evol, e, consolidation):
    
    size_titulo = 9.5
    size_ejes   = 7
    size_leyenda = 8
    
    q_data_0, e1_data_0, pE_data_0, eV_data_0, e_data_0 = np.loadtxt('Necessary File\q-e1-pE-eV-e_50.txt', unpack=True)
    q_data_1, e1_data_1, pE_data_1, eV_data_1, e_data_1 = np.loadtxt('Necessary File\q-e1-pE-eV-e_100.txt', unpack=True)
    q_data_2, e1_data_2, pE_data_2, eV_data_2, e_data_2 = np.loadtxt('Necessary File\q-e1-pE-eV-e_200.txt', unpack=True)
    
    plt.subplots()
    
    conso = len(consolidation)
    
    for i in range(0,conso):
        if i == 0:
            q_0, e1_0, pE_0, eV_0, e_0 = q[:,0], Eps1[:,0], pE[:,0], Evol[:,0], e[:,0]
            l1, l2 = Triaxial_0(q_0, e1_0, pE_0, eV_0, e_0, q_data_0, e1_data_0, pE_data_0, eV_data_0, e_data_0, size_titulo, size_ejes, size_leyenda)
            leyenda_1 = ["Plaxis-50kPa", "Test-50kPa"]
            plt.legend([l1, l2], labels=leyenda_1, bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0, fontsize = size_leyenda)
        if i == 1:
            q_1, e1_1, pE_1, eV_1, e_1 = q[:,1], Eps1[:,1], pE[:,1], Evol[:,1], e[:,1]
            l3, l4 = Triaxial_1(q_1, e1_1, pE_1, eV_1, e_1, q_data_1, e1_data_1, pE_data_1, eV_data_1, e_data_1, size_titulo, size_ejes, size_leyenda)
            leyenda_2 = ["Plaxis-100kPa", "Test-100kPa"]
            plt.legend([l3, l4], labels=leyenda_2, bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0, fontsize = size_leyenda)
        if i == 2:
            q_2, e1_2, pE_2, eV_2, e_2 = q[:,2], Eps1[:,2], pE[:,2], Evol[:,2], e[:,2]
            l5, l6 = Triaxial_2(q_2, e1_2, pE_2, eV_2, e_2, q_data_2, e1_data_2, pE_data_2, eV_data_2, e_data_2, size_titulo, size_ejes, size_leyenda)
            leyenda_3 = ["Plaxis-200kPa", "Test-200kPa"]
            plt.legend([l5, l6], labels=leyenda_3, bbox_to_anchor=(1.05, 0), loc=3, borderaxespad=0, fontsize = size_leyenda)
    
    #Save graph
    plt.subplots_adjust(top=0.94, bottom=0.1, left=0.13, right=0.77, hspace=0.43, wspace=0.33)

    plt.savefig("Results\Triaxial Test-Plaxis" +".png")