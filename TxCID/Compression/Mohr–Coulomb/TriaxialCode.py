from math import sin, pi

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

def make_geometry(g_i, width, relacion, str_ax):
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
    
    line_strain_h = g_i.line((0, high), (ancho, high))[-1]
    g_i.linedispl(line_strain_h, "uy_start", -high*(str_ax+.01))
    
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
    
    soil_params = [("MateriaLName", suelo["name"]),
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
    
    
def phase_construction(g_i):
    
    # set up calculation
    
    g_i.gotostages()
    
    top_load   = g_i.LineLoad_1_1.qy_start
    right_load = g_i.LineLoad_2_1.qx_start
    top_strain = g_i.LineDisplacement_3_1
    
    # initial phase
    
    model_conditions = (g_i.GroundwaterFlowBCs, g_i.Deformations, g_i.Dynamics, g_i.FieldStress, g_i.GroundwaterFlow, g_i.Water)
    
    inicial = g_i.InitialPhase
    g_i.InitialPhase.DeformCalcType = "K0 procedure"
    g_i.activate((g_i.Soils, g_i.lines), g_i.InitialPhase)
    g_i.deactivate(g_i.lineloads, g_i.InitialPhase)
    g_i.deactivate(top_strain, g_i.InitialPhase)
    g_i.deactivate(model_conditions, g_i.InitialPhase)
    
    return top_strain, top_load, right_load

def conso_shear(g_i, g_o, punto, suelo, consolidation, k, i, top_strain, top_load, right_load):
    
    eo = suelo["e_o"]   #Void Ratio Initial

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
        
    # shear phase
        
    shear_phase = g_i.phase(consolidation_phase)
    shear_phase.DeformCalcType = "Plastic"
    shear_phase.Identification = "Corte [Consolidacion " + str(conso) + "kPa]"
    shear_phase.Deform.IgnoreUndrainedBehaviour = True
    shear_phase.Deform.ResetDisplacementsToZero = True
    shear_phase.Deform.UseDefaultIterationParams = False
    shear_phase.Deform.ToleratedError = 0.0001
    shear_phase.Deform.MaxLoadFractionPerStep = 0.01
    g_i.activate(top_strain, shear_phase)

    g_i.calculate()
    
    g_i.view(consolidation_phase)
    Evol_conso   = []
    Evol_conso.append(g_o.getcurveresults(punto, g_o.Phases[k-1], g_o.ResultTypes.Soil.TotalVolumetricStrain))
    Evol_consolidation = Evol_conso[0]
    Evol_consolidation = float(Evol_consolidation)
    e_phase =  Evol_consolidation*(1+eo) +  eo
    
    q_o = [0]
    Eps1_o  = [0]
    pE_o = [-conso]
    Evol_o = [0]
    
    #Fases a considerar
    phaseorder  = [g_o.Phases[k]]

    for phase in phaseorder:
        for step in phase.Steps:
            q_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.DeviatoricStress))
            Eps1_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.Eps1))
            pE_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.MeanEffStress))
            Evol_o.append(g_o.getcurveresults(punto, step, g_o.ResultTypes.Soil.TotalVolumetricStrain))
    
    filas = round(len(q_o),0)
    
    q_total    = np.zeros([filas, 1])
    Eps1_total = np.zeros([filas, 1])
    pE_total   = np.zeros([filas, 1])
    Evol_total = np.zeros([filas, 1])
    e_total    = np.zeros([filas, 1])

    q_total[:,0]    = q_o
    Eps1_total[:,0] = Eps1_o
    pE_total[:,0]   = pE_o
    Evol_total[:,0] = Evol_o
    e_total[:,0]    = Evol_total[:,0]*(1+eo) + e_phase
    
    with open("Results\q-e1-pE-eV-e [" + str(conso) + "kPa]-PLAXIS" + ".txt", "w") as file:
        for j in range(filas):
            datos = "{:.2f}\t{:.5f}\t{:.2f}\t{:.5f}\t{:.5f}\n".format(q_total[j,0], Eps1_total[j,0], pE_total[j,0], Evol_total[j,0], e_total[j,0])
            file.writelines(datos)

    g_o.update()