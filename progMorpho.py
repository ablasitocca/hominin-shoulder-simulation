### Combinaison des changements morpho
# 4 paramètres

import opensim as osim
import numpy as np
import math
import os
from BiblioMorpho import *

#osim.LoadOpenSimLibrary('C:\OpenSim 4.1\plugins\ScapulothoracicJointPlugin40_WinX64')

fileMorpho = "paramMorpho.txt"

################################REMARQUES##################################################
#Les changements morpho du thorax s'applique dans un second temps les modèles pour lesquels
#il faut appliquer ces changements
###########################################################################################
## import des valeurs de morphologie #MAJ 2023
CL = choix_parm(fileMorpho)[0] #longueur clavicule (coeff mult)
HU = choix_parm(fileMorpho)[1] #torsion humérale (°)
GL = choix_parm(fileMorpho)[2] #orientation glen (°)



#####Applications des changements morpho de clavicule, humérus et orientation glen #MAJ 2023 #MAJ MoBL
for cl in CL:
    for hu in HU:
        for gl in GL:
            ## Définition de new model sur la base du model ref
            modelName = 'MoBL_ARMS_43.osim'
            model = osim.Model(modelName)
            newmodel = osim.Model(model)

            ## Application des changements morpho
            clav_long(modelName,newmodel,'setup_scaleClav.xml',cl)
            torsion_hum(newmodel,hu)
            orientation_glen(model,newmodel,gl)
            
            ## Création du nouveau model
            newmodel.setName('model_{}{}{}'.format(CL.index(cl),HU.index(hu),GL.index(gl)))
            newmodel.printToXML('model_{}{}{}.osim'.format(CL.index(cl),HU.index(hu),GL.index(gl)))
