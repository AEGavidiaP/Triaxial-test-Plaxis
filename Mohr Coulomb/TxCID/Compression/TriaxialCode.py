from math import sin, pi

import matplotlib.pyplot as plt
import numpy as np

#Define function model properties

def properties_model(g_i, data):
    
    soil_configuration = [("Title", data["titulo"]),
                          ("Company", "Universidad de Chile"),
                          ("Comments", data["comentarios"]),
                          ("UnitForce", data["unidad_fuerza"]),
                          ("UnitLength", data["unidad_longitud"]),
                          ("UnitTime", data["tiempo"]),
                          ("UnitTemperature", "K"),
                          ("UnitEnergy", "kJ"),
                          ("UnitPower", "kW"),
                          ("WaterWeight", data["pesounitario_agua"]),
                          ("ReferenceTemperature", 293.15),
                          ("LiquidSpecificHeatCapacity", 4181.3),
                          ("LiquidThermalConductivity", 0.0006),
                          ("LiquidLatentHeat", 334000),
                          ("LiquidThermalExpansion", 0.00021),
                          ("LiquidTemperature", 293.15),
                          ("IceSpecificHeatCapacity", 2108),
                          ("IceThermalConductivity", 0.00222),
                          ("IceThermalExpansion", 5E-5),("VapourSpecificHeatCapacity", 1930),
                          ("VapourThermalConductivity", 2.5E-5),
                          ("VapourSpecificGasConstant", 461.5),
                          ("UseTemperatureDepWaterPropsTable", False),
                          ("ModelType", data["modelo"]),
                          ("ElementType", data["nodos"])]
    
    return g_i.setproperties(*soil_configuration)
    
#Define function to create polygon - Triaxial test (specimen)

def make_geometry(g_i, width, relacion):
    ancho = width
    high  = ancho*relacion*2
    medio = high/2
    
    top_points = []
    top_points.append([0, 0])
    top_points.append([ancho, 0])
    top_points.append([ancho, high])
    top_points.append([0, high])
    
    probeta = g_i.polygon(*top_points)
    
    line_v = g_i.line((0, 0), (0, high))[-1]
    g_i.linedispl(line_v, "Displacement_x", "Fixed", "Displacement_y", "Free")
    
    line_h = g_i.line((0, 0), (ancho, 0))[-1]
    g_i.linedispl(line_h, "Displacement_x", "Free", "Displacement_y", "Fixed")
    
    line_load_v = g_i.line((0, high), (ancho, high))[-1] 
    g_i.lineload(line_load_v, "qy_start", -1, "qx_start", 0)
    
    line_load_h = g_i.line((ancho, 0), (ancho, high))[-1] 
    g_i.lineload(line_load_h, "qx_start", -1, "qy_start", 0)
    
    return probeta, medio
    
 #Define function to create constitutive model - Mohr-Coulomb (drain)

def soilmaterial_mc_drain(g_i, suelo):

    soilmodel = 2    # Modelo Mohr-Coulomb se utiliza el valor "2"
    drainagetype = 0 #Tipo de drenaje, en el material drenado se utiliza el valor de "0"
    gamma = 0        #Para simular ensayos triaxiales el material se debe creo con un peso unitario igual a c
    
    soil_params = [("MaterialName", suelo["name"]),
                  ("SoilModel", soilmodel),
                  ("DrainageType", drainagetype),
                  ("gammaUnsat", gamma),
                   ("gammaSat", gamma),
                  ("Eref", suelo["E"]),
                  ("nu", suelo["nu"]),
                  ("cref", suelo["c"]),
                  ("phi", suelo["fi"]),
                  ("psi", suelo["dilatancia"])]

    soil_material = g_i.soilmat(*soil_params)
    
    return soil_material
    
def mesh_point(g_i, g_o, medio):
    
    g_i.gotomesh()
    g_i.mesh()
    
    #Select point
    output_port = g_i.selectmeshpoints()
    
    x_i= 0
    y_i= medio
    
    punto = g_o.addcurvepoint("stress point", (x_i, y_i))
    punto.Identification = "Central"
    
    curvepoints_central = g_o.CurvePoints.StressPoints
    print(g_o.tabulate(curvepoints_central, "x y"))
    
    g_o.update()
    
    return punto
    
    
def phase_construction(g_i, g_o, punto, suelo, consolidation):

    top_load   = g_i.LineLoad_1_1.qy_start
    right_load = g_i.LineLoad_2_1.qx_start
       
    # set up calculation
    
    g_i.gotostages()
    
    # initial phase
    
    model_conditions = (g_i.GroundwaterFlowBCs, g_i.Deformations, g_i.Dynamics, g_i.FieldStress, g_i.GroundwaterFlow, g_i.Water)
    
    inicial = g_i.InitialPhase
    g_i.InitialPhase.DeformCalcType = "K0 procedure"
    g_i.activate((g_i.Soils, g_i.lines), g_i.InitialPhase)
    g_i.deactivate(g_i.lineloads, g_i.InitialPhase)
    g_i.deactivate(model_conditions, g_i.InitialPhase)
    
    i = 0
    contador = []
    
    count_test = len(consolidation)
    
    for k in range(0,count_test):
        k = 2 + k*2
        
        conso = consolidation[i]
           
        # consolidation phase
            
        consolidation_phase = g_i.phase(g_i.InitialPhase)
        consolidation_phase.DeformCalcType = "Plastic"
        consolidation_phase.Identification = "Consolidacion [" + str(conso) + "kPa]"
        consolidation_phase.Deform.IgnoreUndrainedBehaviour = True
        g_i.activate((g_i.lineloads), consolidation_phase)
        top_load.set(consolidation_phase, -conso)
        right_load.set(consolidation_phase, -conso)
            
        # Calculo valor para el corte 
            
        angulo = suelo["fi"]*pi/180
        m_c    = 6*sin(angulo)/(3-sin(angulo))
        load_axial_failure  = (conso*((2/3)*m_c+1))/(1-(m_c/3))
        corte = 1.1*load_axial_failure
            
        # shear phase
            
        shear_phase = g_i.phase(consolidation_phase)
        shear_phase.DeformCalcType = "Plastic"
        shear_phase.Identification = "Corte [Consolidacion " + str(conso) + "kPa]"
        shear_phase.Deform.IgnoreUndrainedBehaviour = True
        shear_phase.Deform.ResetDisplacementsToZero = True
        top_load.set(shear_phase, -corte)
    
        g_i.calculate()
        g_i.view(shear_phase)
        q_o = [0]
        Eps1_o  = [0]
        pE_o = [-conso]
        Evol_o = [0]
        
        #Fases a considerar
        phaseorder = [g_o.Phases[k]]

        for phase in phaseorder:
            for step in phase.Steps:
                q_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.DeviatoricStress))
                Eps1_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.Eps1))
                pE_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.MeanEffStress))
                Evol_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.TotalVolumetricStrain))
        
        filas = round(len(q_o)*1.1,0)
        filas = int(filas)
        
        eo = suelo["e_o"]
        
        if i==0:
            q_total    = np.zeros([filas, count_test])
            Eps1_total = np.zeros([filas, count_test])
            pE_total   = np.zeros([filas, count_test])
            Evol_total = np.zeros([filas, count_test])
            e_total    = np.zeros([filas, count_test])
        
        valor = len(q_total[:,0]) - len(q_o)
        for j in range(valor):
            q_o.append(q_o[-1])
            Eps1_o.append(Eps1_o[-1])
            pE_o.append(pE_o[-1])
            Evol_o.append(Evol_o[-1])
            j += 1
    
        q_total[:,i]    = q_o
        Eps1_total[:,i] = Eps1_o
        pE_total[:,i]   = pE_o
        Evol_total[:,i] = Evol_o
        e_total[:,i]    = Evol_total[:,i]*(1+eo) + eo

        i += 1 
            
        with open("q - e1 - p' - eV [" + str(conso) + "kPa]-PLAXIS" + ".txt", "w") as file:
            file.writelines(["{}\t{}\t{}\t{}\n".format("q[kPa]" ,"e1[-]", "p'[kPa]", "eV[-]")])
            file.writelines(["{:.2f}\t{:.5f}\t{:.2f}\t{:.5f}\n".format(DeviatoricStress, Eps1, MeanEffStress, TotalVolumetricStrain)
                             for DeviatoricStress, Eps1, MeanEffStress, TotalVolumetricStrain in zip(q_o, Eps1_o, pE_o, Evol_o)])
        
        g_o.update()
    
    return q_total, Eps1_total, pE_total, Evol_total, e_total
    
    
def graphics_triaxial(q, Eps1, pE, Evol, e):
        
        #e1-q
        plt.subplot(221)
        plt.plot(Eps1[:,0], q[:,0], "b")
        plt.plot(Eps1[:,1], q[:,1], "r")
        plt.plot(Eps1[:,2], q[:,2], "m")
        x1_0, y1_0 = np.loadtxt('e1-q_50.txt', unpack=True)
        x1_1, y1_1 = np.loadtxt('e1-q_100.txt', unpack=True)
        x1_2, y1_2 = np.loadtxt('e1-q_200.txt', unpack=True)
        plt.plot(x1_0,y1_0, "k")
        plt.plot(x1_1,y1_1, "brown")
        plt.plot(x1_2,y1_2, "grey")
        plt.title('Carga-Deformacion')
        plt.xlabel('$\epsilon_1$')
        plt.ylabel('q')
        plt.xlim(-0.2,0)
        plt.grid(True)
        
        #p-q
        plt.subplot(222)
        plt.plot(pE[:,0], q[:,0], "b")
        plt.plot(pE[:,1], q[:,1], "r")
        plt.plot(pE[:,2], q[:,2], "m")
        x2_0, y2_0 = np.loadtxt('pE-q_50.txt', unpack=True)
        x2_1, y2_1 = np.loadtxt('pE-q_100.txt', unpack=True)
        x2_2, y2_2 = np.loadtxt('pE-q_200.txt', unpack=True)
        plt.plot(x2_0,y2_0, "k")
        plt.plot(x2_1,y2_1, "brown")
        plt.plot(x2_2,y2_2, "grey")
        plt.title('p\'-q')
        plt.xlabel('p\'')
        plt.ylabel('q')
        plt.grid(True)
        
        #e1-ev
        plt.subplot(223)
        plt.plot(Eps1[:,0], Evol[:,0], "b")
        plt.plot(Eps1[:,1], Evol[:,1], "r")
        plt.plot(Eps1[:,2], Evol[:,2], "m")
        x3_0, y3_0 = np.loadtxt('e1-eV_50.txt', unpack=True)
        x3_1, y3_1 = np.loadtxt('e1-eV_100.txt', unpack=True)
        x3_2, y3_2 = np.loadtxt('e1-eV_200.txt', unpack=True)
        plt.plot(x3_0,y3_0, "k")
        plt.plot(x3_1,y3_1, "brown")
        plt.plot(x3_2,y3_2, "grey")
        plt.title('Deformacion axial-volumetrica')
        plt.xlabel('$\epsilon_1$')
        plt.ylabel('$\epsilon_V$')
        plt.xlim(-0.2,0)
        plt.grid(True)
        
        #p-e
        plt.subplot(224)
        plt.plot(pE[:,0], e[:,0], "b")
        plt.plot(pE[:,1], e[:,1], "r")
        plt.plot(pE[:,2], e[:,2], "m")
        x4_0, y4_0 = np.loadtxt('pE-e_50.txt', unpack=True)
        x4_1, y4_1 = np.loadtxt('pE-e_100.txt', unpack=True)
        x4_2, y4_2 = np.loadtxt('pE-e_200.txt', unpack=True)
        plt.plot(x4_0,y4_0, "k")
        plt.plot(x4_1,y4_1, "brown")
        plt.plot(x4_2,y4_2, "grey")
        plt.title('p\'-e')
        plt.xlabel('p\'')
        plt.ylabel('e')
        plt.grid(True)
        
        #Guardar grafico
        plt.subplots_adjust(top=0.92, bottom=0.1, left=0.13, right=0.95, hspace=0.4, wspace=0.45)
        plt.savefig("Comparativa ensayos-modelo" +".png")
    