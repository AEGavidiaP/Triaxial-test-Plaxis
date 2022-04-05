import TriaxialCode
import Graphics

from plxscripting.easy import *

l_input = 10000
l_ouput = 10001
pw = 'ixcMy6nxXr5/3iA%' #Cambiar segun usuario de PLAXIS

s_i, g_i = new_server('localhost', l_input, password=pw)
s_o, g_o = new_server('localhost', l_ouput, password=s_i.connection._password)

s_i.new()

# Modify model properties

model = {"titulo": "Calibracion del modelo",
        "comentarios": "Realizar calibracion de TxCID en compresion para modelo Mohr-Coulomb",
        "unidad_fuerza": "kN", #Trabajar con kiloNewton -> "kN"; trabajar con MegaNewton -> "MN"
        "unidad_longitud": "m", #Trabajar en milimetro -> "mm"; trabajar en centimetro -> "cm"; trabajar en metro -> "m"; trabajar en kilometro -> "m"
        "tiempo": "s", #Trabajar en segundos -> "s"; trabajar en minutos -> "min"; trabajar en horas -> "h"; trabajar en dias -> "day"
        "pesounitario_agua": 9.8, #Tener encuenta las unidades de fuerza y longitud: (unidad_fuerza)/((unidad_longitud)^3)
        "modelo": "Axisymmetry", #Trabajar en deformaciones planas ->"Planestrain"; trabajar con modelo axisimitrico "Axisymmetry"
        "nodos": "15noded" #Numero de nodos a trabajas: 6 nodos -> "6noded"; 15 nodos -> "15noded"
        }


#Define size (test specimen)

width = 1 #Specimen width
ratio = 2 #H/D ratio of the specimen
consolidation = [50, 100, 200] #Confinements

#Soil properties & assign material

suelo = {"name": "Suelo 1", #Nombre del suelo a evaluar
         "E": 20000,        #Modulo de elasticidad
         "nu": 0.3,         #Modulo de Poisson
         "c": 0,            #Cohesion
         "fi": 32,          #angulo de friccion
         "dilatancia": 0,   #angulo de dilatancia
         "e_o":0.76         #Indice de vacios inicial
        }

#Document properties        
TriaxialCode.properties_model(g_i, model)

#Specimen geometry
probeta_ensayo, medio = TriaxialCode.make_geometry(g_i, width, ratio)

#Constitutive Model
ensayo = TriaxialCode.soilmaterial_mc_drain(g_i, suelo)

#Assign constitutive model to the specimen
probeta_ensayo = g_i.setmaterial(ensayo)

#Generate & select mesh points
punto = TriaxialCode.mesh_point(g_i, g_o, medio)

#Construction phases for the triaxial test
q, Eps1, pE, Evol, e = TriaxialCode.phase_construction(g_i, g_o, punto, suelo, consolidation)

#Plot triaxial in 4 planes (Save in png)
Graphics.graphics_triaxial(q, Eps1, pE, Evol, e, consolidation)