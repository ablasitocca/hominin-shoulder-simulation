#### Bibliothèque des fonctions de changement morpho  ####
import opensim as osim
import numpy as np
import math
import xml.etree.ElementTree as ET

##############################################################################################################
####################### Changer la torsion humerale ##########################################################
##############################################################################################################
###### Changer insertion musculaire 
###### Changer orientation du coude
def calcul_torsion(x,y,z,torsion):
    """calcul les nouvelles coordonées des points d'insertion. pas de progressivité : toute la torsion est appliquée en prox
input : x : coordonnée cartésienne dans le repère de l'humérus (m)
        y : coordonnée cartésienne dans le repère de l'humérus (m)
        z : coordonnée cartésienne dans le repère de l'humérus (m)
        torsion : angle à ajouter à la torsion du modèle en °
output : x2,y2 nouvelles coordonnées cartésiennes de l'insertion avec progressivité linéaire de torsion
        theta2 différence entre les valeurs de torsion en radian avec progressivité linéaire de torsion"""

    #calcul des grandeurs initiales 
    r0 = math.sqrt(x**2 + z**2) #rayon de la section de l'humérus
    theta = math.acos(x/r0) #angle initiale #radians

    #calcul de l'angle à ajouter à l'angle d'insertion
    theta2 = ((torsion*np.pi)/180)#*coeff #en radian
    

    #calcul des nouvelles coordonnées cartésiennes
    x2 = r0*math.cos(theta-theta2) #angle en radian 
    z2 = r0*math.sin(theta-theta2)  

    return x2,z2,theta2


def torsion_hum(newmodel,torsion):
    """applique la nouvelle torsion humérale sur le modèle sélectionné : change insertions muscu et orientation avant-bras
input : newmodel : nom du model qui acceuillera les nouvelles props
        torsion: angle à ajouter à la torsion du modèle en °
output : retourne nouvelles prop du model"""
    ## liste contenant les muscles dont il faut changer l'insertion
    # coracobrachial, grand rond, pec maj, deltoide, grand dorsal
    #muscles = ['Coracobrachialis','TeresMajor','PectoralisMajorClavicle_S','PectoralisMajorThorax_I','PectoralisMajorThorax_M',
               #'DeltoideusClavicle_A','DeltoideusScapula_P','DeltoideusScapula_M','LatissimusDorsi_S','LatissimusDorsi_M',
               #'LatissimusDorsi_I','Brachial']
    muscles = ['PectoralisMajorClavicle_S','PectoralisMajorThorax_I','PectoralisMajorThorax_M','TeresMajor'] #'TMAJ','PECM1','PECM2','PECM3',

    ## variable state du model
    state = newmodel.initSystem()

    ## Changer insertion
    # accès au muscle
    for m in muscles:
        muscle = newmodel.getMuscles().get(m) #muscle class
        pps = muscle.getGeometryPath().getPathPointSet() #get path point set//tous les pts origin/insert

        #selction du la bonne insertion
        for i in range (0,pps.getSize()) : #parcourt les points d'insertion pour trouver celui sur l'humérus
            app = pps.get(i) #abstract path point = point d'insertion
            if app.getParentFrame().getName() == 'humerus': #on ne considère que les insertions sur l'humérus
                loc = app.getLocation(state)
                x,y,z = loc[0],loc[1],loc[2]

                #new location
                x2, z2, theta2 = calcul_torsion(x,y,z,torsion)
                newloc = osim.Vec3(x2,y,z2)

                #donne la nouvelle location au model
                ppt = osim.PathPoint.safeDownCast(app) #permet d'utiliser la fonction setLocation (pas utilisable avec les abastract path point)
                ppt.setLocation(newloc)
                


    ## Changer orientation avant-bras
    # attribution de l'angle maximal de torsion pour l'avant-bras
    theta2 = (torsion*np.pi)/180 #radian

    # accès au coude
    joint = newmodel.getJointSet().get('elbow') #joint class
    pof = joint.get_frames(0) #physical offset frame #humerus offset
    ortn = pof.get_orientation() #accès au vecteur orientation du modèle
    rx,ry,rz = ortn[0],ortn[1],ortn[2] #rotation en radian
    
    # changement de l'orientation
    ry2 = ry + theta2  #radian
    newOrtn = osim.Vec3(rx,ry2,rz)
    pof.set_orientation(newOrtn)
    
    return newmodel.initSystem()




##############################################################################################################
####################### Changer longueur clavicule ###########################################################
##############################################################################################################
##### Scale de la clavicule
def clav_long(model,newmodel,genericSetupScale,coeff_long):
    """change la longueur de la clavicule. Applique un scale factor (pour changer les insertions) puis change contrainte AC (pour déplacer les autres articulations)
input : Model : entrer nom du fichier .osim à mettre à l'échelle (modèle de référence)
        newmodel : model qui accueillera les nouvelles props
        genericSetupScale : fichier .xml de setup générique, qui contient par exemple le poids des marqueurs
        coeff_long : coeffficient de la longueur de clav (1 garde la même taille)
output : retourne les nouvelles prop du modèle"""
    ## Scale clavicule
    # model state
    state = newmodel.initSystem()
    newmodel.initSystem()

    #initialize scale tool
    scaleTool = osim.ScaleTool(genericSetupScale)

    #tell scale tool to use the model
    scaleTool.getGenericModelMaker().setModelFileName(model)

    #ScaleTOOL 
    scaleTool.getModelScaler().setApply(True)
    scaleTool.getModelScaler().processModel(newmodel)
    scaleTool.getModelScaler().setPreserveMassDist(True)

    # scale factor de la clavicule
    scale = scaleTool.getModelScaler().getScaleSet().get(0) #scale
    SF = scale.getScaleFactors() #scale factor clavicule
    x,y,z = SF[0],SF[1],SF[2]

    # nouveau scale factor
    z = coeff_long*z
    new_scale = osim.Vec3(x,y,z)
    scale.setScaleFactors(new_scale)

    #lancer le scale
    path2subject = scaleTool.getPathToSubject() ##ajouté le 09/07
    scaleTool.getModelScaler().processModel(newmodel,path2subject,26.2) ##ajouté le 09/07


    ## Changement de la contrainte AC
    # Accès à la contrainte AC
    #AC = newmodel.getConstraintSet().get('AC')

    # Accès à la location du point de contrainte
    #AC = osim.PointConstraint.safeDownCast(AC)
    #loc = AC.get_location_body_1()  
    #x,y,z = loc[0],loc[1],loc[2]

    # Changement de Z en fonction du coeff d'agrandissement de la clavicule
    #z = coeff_long*z
    #new_loc = osim.Vec3(x,y,z)

    # Attribution de la nouvelle location
    #AC.set_location_body_1(new_loc)
    
    return scaleTool.run(),newmodel.initSystem() #### scaleTool.run() ajouté le 09/07



##############################################################################################################
####################### Changer orientation glénoïde ######################################################### #ajout 2023
##############################################################################################################
def new_coords_euler_angles(alpha,A):
    """fonction qui calcule la nouvelle position B du centre articulaire A de l'articulation GH
après application d'une crânialisation d'angle alpha. la fonction exprime la position B dans le
repère fixe R0 utilisable par Opensim. Le calcul de la nouvelle position B se fait dans Rg, le
repère attaché à la glénoïde tel que le vecteur de translation entre les origines R0 et Rg O0-OG
vaut(-0.03,-0.04,-0.015) et les rotations pour passer de R0 à Rg : (0.15,0.6,-0.45)
input : alpha = flottant, valeur en ° ; A = vecteur dim 3, coordonnées du point A exprimées dans R0 (repère fixe
        de référence accroché à l'acromion). 
output : retourne les coordonnées de B exprimées dans R0 (directement utilisatible dans OpenSim)
"""
    ##conversion alpha en rad
    alpha = (alpha*np.pi)/180
    
    ## coordonnées de Og (origine de Rg) exprimées dans R0 #déterminé dans Opensim
    xg_r0,yg_r0,zg_r0 = -0.029, -0.022, -0.028 #-0.03,-0.04,-0.015

    ##coordonées du point A expirmées dans R0 #coordonnées données par scapula_offset de GH joint
    xa_r0, ya_r0,za_r0 = A[0],A[1],A[2]

    ## rotations pour passer de R0 -> Rg # angles d'Euler #déterminé dans Opensim
    rot_x,rot_y,rot_z = 0,0,-0.55 #0.15,0.6,-0.45

    ##vecteur OgA exprimé dans R0 #vecteur translation
    xOgA_r0,yOgA_r0,zOgA_r0 = xa_r0 - xg_r0, ya_r0 - yg_r0, za_r0 - zg_r0

    ##matrice rotation de R0 -> RG avec la séquence X(rot_x)Y(rot_y)Z(rot_z)
    #exprime RG dans le repère R0
    MR = np.zeros((3,3))
    MR[0,0] = np.cos(rot_y)*np.cos(rot_z)
    MR[1,0] = np.sin(rot_x)*np.sin(rot_y)*np.cos(rot_z)+np.cos(rot_x)*np.sin(rot_z) #2ème ligne, première colonne
    MR[2,0] = np.sin(rot_x)*np.sin(rot_z)-np.cos(rot_x)*np.sin(rot_y)*np.cos(rot_z)
    MR[0,1] = -np.sin(rot_z)*np.cos(rot_y)
    MR[1,1] = -np.sin(rot_x)*np.sin(rot_y)*np.sin(rot_z)+np.cos(rot_x)*np.cos(rot_z)
    MR[2,1] = np.sin(rot_x)*np.cos(rot_z)+np.cos(rot_x)*np.sin(rot_y)*np.sin(rot_z)
    MR[0,2] = np.sin(rot_y)
    MR[1,2] = -np.sin(rot_x)*np.cos(rot_y)
    MR[2,2] = np.cos(rot_x)*np.cos(rot_y)

    ##ETAPE 1 : calcul des coordonnées de A exprimées dans RG 
    #application de la translation
    xa_r0t,ya_r0t,za_r0t = xOgA_r0,yOgA_r0,zOgA_r0

    #application de la rotation : ATTENTION, c'est la matrice transposée de MR
    #xa_rg, ya_rg et za_rg = coordonnées de A exprimées dans Rg
    xa_rg = MR[0,0]*xa_r0t + MR[1,0]*ya_r0t + MR[2,0]*za_r0t ##matrice transposée
    ya_rg = MR[0,1]*xa_r0t + MR[1,1]*ya_r0t + MR[2,1]*za_r0t
    za_rg = MR[0,2]*xa_r0t + MR[1,2]*ya_r0t + MR[2,2]*za_r0t

    ## ETAPE 2 : calcul des coordonnées du nouveau centre articulaire de GH exprimée dans Rg: B_Rg
    d = np.sqrt(xa_rg**2+ya_rg**2+za_rg**2) #distance OG-A #distance entre la glénoïde et le centre articulaire de GH
    xb_rg = xa_rg
    yb_rg = d*np.cos(np.arctan(za_rg/ya_rg)+alpha)
    zb_rg = d*np.sin(np.arctan(za_rg/ya_rg)+alpha)

    ##ETAPE 3 : exprimer B dans R0
    #rotation #ATTENTION cette fois c'est bien la matrice MR non transposée
    xb_r0st = MR[0,0]*xb_rg + MR[0,1]*yb_rg + MR[0,2]*zb_rg 
    yb_r0st = MR[1,0]*xb_rg + MR[1,1]*yb_rg + MR[1,2]*zb_rg
    zb_r0st = MR[2,0]*xb_rg + MR[2,1]*yb_rg + MR[2,2]*zb_rg    

    #translation
    xb_r0 = xb_r0st + xg_r0
    yb_r0 = yb_r0st + yg_r0
    zb_r0 = zb_r0st + zg_r0
    
    return xb_r0,yb_r0,zb_r0


def orientation_glen(model,newmodel,alpha):
    """change orientation de la glénoïde
input : Model : entrer nom du fichier .osim = modèle de référence
        newmodel : model qui accueillera les nouvelles props
        alpha : d° à soustraire au model de ref pour atteindre la nouvelle orientation de la cavité glénoïde
output : retourne les nouvelles prop du modèle"""
    ## init model
    # model state
    state = newmodel.initSystem()
    newmodel.initSystem()

    ##accès aux coordonnées du centre de rotation de glenohumeral joint A exprimée dans R0 (scapula offset frame)
    #joint = newmodel.getJointSet().get('GlenoHumeral') #joint class
    joint = newmodel.getJointSet().get('shoulder0') #joint class
    pof = joint.get_frames(0) #scapphant_offset
    A_R0 = pof.get_translation() #coordonnées du point A exprimées dans R0 = translation de scapula_offset

    ## calcul des coordonnées du nouveau centre articulaire B exprimées dans R0
    B_R0 = new_coords_euler_angles(alpha,A_R0)
    xb_r0,yb_r0,zb_r0 = B_R0[0],B_R0[1],B_R0[2]
    B_R0 = osim.Vec3(xb_r0,yb_r0,zb_r0)

    ## application du nouveau centre articulaire
    pof.set_translation(B_R0)

    return newmodel.initSystem() 



##############################################################################################################
####################### Choix paramètres #####################################################################
##############################################################################################################
def choix_parm(file):
    """retourne les valeurs des paramètres morphologiques en fonction du fichier paramMorpho.txt
input : file paramMorpho.txt
output : les listes des valeurs des paramètres morpho"""
    CL = [] #longueur clavicule (coeff mult)
    HU = [] #torsion humérale (°)
    TH = [] #largeur du thorax (coeff mult)

    fpara = open(file,'r')
    fpara.readline()
    n = int(fpara.readline().split(' ')[2])

    lCL = fpara.readline().split('\t')
    lHU = fpara.readline().split('\t')
    lTH = fpara.readline().split('\t')

    for i in range (1,n+1):
        CL.append(float(lCL[i]))
        HU.append(float(lHU[i]))
        TH.append(float(lTH[i]))

    fpara.close()

    return [CL, HU, TH]

###############################################
################# Validation ##################
###############################################
##### Validation torsion
##print(calcul_torsion(1,0,150)) #donne les résultats attendus
##
##osim.LoadOpenSimLibrary('C:\OpenSim 4.1\plugins\ScapulothoracicJointPlugin40_WinX64')
##model = 'Seth41_modelRef.osim'
##newmodel = osim.Model(model)
##    
##### Validation torsion_hum
##torsion_hum(newmodel,20)
##newmodel.printToXML('Seth41_humTorsion_API.osim')
##
### Validation longCla_AC
##clav_long(model,newmodel,'setup_scaleClav.xml',0.5)
##newmodel.printToXML('Seth41_claLong_API.osim')
##
####### Validation orientation_glen
#alpha = -35
#A = [0,0,0] #[-0.00954994, -0.0340014, 0.00899682]
#print(new_coords_euler_angles(alpha,A))
##
##orientation_glen(model,newmodel,alpha)
##newmodel.setName('model_testglen.osim')
##newmodel.printToXML('model_testglen.osim')

##### Validation choix des paramètres
##print(choix_parm('paramMorpho.txt')[1])
