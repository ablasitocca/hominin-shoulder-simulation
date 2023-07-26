### bibliothèque de fonctions
import opensim as osim
import os
from scipy import signal
import numpy as np
from decimal import *
import matplotlib.pyplot as plt

########################################################################################
############################## SCALE TOOL ##############################################
########################################################################################
def adjustMarker(modelPath,trcDataFileName,genericSetupScale,outputFileName): #MAJ 2020 #ok 2023
    """ajuste les marqueurs du modèle. Pas de scale.
input : modelPath : entrer chemin du fichier .osim à mettre à l'échelle
        trcDataFileName : entrer nom du fichier qui contient essai statique .trc
        genericSetupScale : fichier .xml de setup générique, qui contient par exemple le poids des marqueurs
output : scaled_model.osim dans scaled_model"""
    #definition de l'essai utilisé pour l'ajustement des marqueurs
    E = trcDataFileName.split('_')[2]
    
    modelName = (os.path.split(modelPath)[1]).split(".")[0]
    trcPath = os.path.join('Data',trcDataFileName)

    #load model
    model = osim.Model(modelPath)
    state = model.initSystem()
    model.initSystem()

    #initialize scale tool
    scaleTool = osim.ScaleTool(genericSetupScale)

    #tell scale tool to use the model
    scaleTool.getGenericModelMaker().setModelFileName(modelPath)

    #import marker data and define initial and final time
    markerData = osim.MarkerData(trcPath)
    initial_time = markerData.getStartFrameTime()
    final_time = markerData.getLastFrameTime()

    #create an array double and set as the time range
    TimeArray = osim.ArrayDouble()
    TimeArray.set(0, initial_time)
    TimeArray.set(1,final_time)

    ##ScaleTOOL
    ## 
    scaleTool.getMarkerPlacer().setApply(True)
    scaleTool.getMarkerPlacer().setTimeRange(TimeArray)
    scaleTool.getMarkerPlacer().setStaticPoseFileName(trcPath)
    scaleTool.getMarkerPlacer().setOutputModelFileName(os.path.join('scaled_model','{}_E{}_markerAdjusted.osim'.format(outputFileName,E)))

    scaleTool.setPrintResultFiles(True)
    scaleTool.getMarkerPlacer().setMoveModelMarkers(True)
    scaleTool.setName('{}_E{}_markerAdjusted'.format(modelName,E))

    ##
    path2subject = scaleTool.getPathToSubject()

    ##
    scaleTool.getModelScaler().processModel(model,path2subject, 26.2)

    ## Run Scale Tool
    return scaleTool.run()


########################################################################################
############################## IK ######################################################
########################################################################################
def IKin (Model,trcDataFileName,genericSetupIK,resultsFolder,init_time,end_time,outputFileName) : #MAJ 2020 #ok 2023
    """Calcule la cinématique inverse
input : Model : entrer le nom du fichier modèle .osim
        trcDataFileName : entrer le nom du fichier .trc qui contient la cinématique des marqueurs
        genericSetupIK : entrer le nom du fichier .xml qui contient le setup pour IK général. La fonction entre seule les détails comme temps initial et final.
        resultsFolder : entrer le nom du dossier dans lequel on veut stocker les résultats (.mot et .sto)
output : motion file .mot et .sto dans ik_results"""
    ikTool = osim.InverseKinematicsTool(genericSetupIK)

    trcDataFile = os.path.join('Data',trcDataFileName)
    
    #load model
    model = osim.Model(os.path.join('scaled_model',Model))
    state = model.initSystem()
    model.initSystem()

    #set IK results in a results folder
    ikTool.setResultsDir(resultsFolder)

    #tell IK tool to use the model
    ikTool.setModel(model)

    #IK settings
    #tell the prog to consider trc file as the data file for time range
    markerData = osim.MarkerData(trcDataFile)

    #setup the IK tool
    ikTool.setMarkerDataFileName(trcDataFile)
    ikTool.setStartTime(init_time)
    ikTool.setEndTime(end_time)
    ikTool.setOutputMotionFileName(os.path.join('IK_results','{}_ik.mot'.format(outputFileName))) #stockage du fichier .mot

    #run IK
    return ikTool.run()


########################################################################################
############################## ID ######################################################
########################################################################################
def IDyn (Model,IKResultsFile,genericSetupID,resultsFolder,init_time,end_time,freq,outputFileName): #MAJ 2020 #ok 2023
    """Calcule la dynamique inverse
input : Model = nom du modele
        IKresultsfile = nom du fichier de résultats de l'IK en .mot
        genericSetupID = nom du fichier générique du setup .xml
        resultsFolders = nom du dossier qui va contenir les résultats
        init_time = temps de départ en seconde
        end_time = temps final du calcul en seconde. attention : il faut exclure les instants d'impact dans les calculs. 
        freq = fréquence de coupure
output : fichier des moments articulaires .sto. Les muscles sont exclus du calcul des moments. cf dans le .xml setup."""
    idTool = osim.InverseDynamicsTool(genericSetupID)

    motDataFile = os.path.join('IK_results',IKResultsFile)

    #load model
    model = osim.Model(os.path.join('scaled_model',Model))
    state = model.initSystem()
    model.initSystem()
  
    #set ID results in a results folder
    idTool.setResultsDir(resultsFolder)

    #tell ID tool to use the model
    idTool.setModel(model)

    #il faut exclure les muscles du calcul

    #setup the ID tool
    idTool.setCoordinatesFileName(motDataFile) ##
    idTool.setStartTime(init_time)
    idTool.setEndTime(end_time)
    idTool.setLowpassCutoffFrequency(freq)
    idTool.setOutputGenForceFileName('{}_id.sto'.format(outputFileName))#stockage du fichier .sto

    #run ID
    return idTool.run()


########################################################################################
############################## Filtrage ################################################
########################################################################################
def filtrage(dataInput,order,fc,fs): #MAJ 2020 #ok 2023
    """filtre une ligne de données avec un filtre zero phase (filtfilt)
input : dataInput : ligne de tableau à filtrer
        order : ordre du filtre
        fc : frequence de coupure
        fs : fréquence d'échantillonage
output : ligne de tableau filtré"""
    w = fc / (fs / 2) # Normalize the frequency
    b, a = signal.butter(order, w, 'low') #ordre,freq coup, type, analog ou digital
    out = signal.filtfilt(b,a,dataInput)
    return np.array(out)


def filtrageIK(IKresult,order,fc): #MAJ 2020 #ok 2023
    """filtre les résultats de la cinématique inverse
input : IKresult : nom du fichier de resultat IK .mot
    order : ordre du filtre de butterworth
    fc : fréquence de coupure
    fs : fréquence d'échantillonage (attention th de Shannon)
output : ligne des DDL (avec le temps), tableau des données filtrées"""
    path = os.path.join('IK_results',IKresult)
    
    f = open(path,'r')
    L1,L2,L3,L4 = f.readline(),f.readline(),f.readline(),f.readline()
    nl = int(L3.split('=')[1])
    nc = int(L4.split('=')[1])
    L5,L6,L7,L8,L9,L10 = f.readline(),f.readline(),f.readline(),f.readline(),f.readline(),f.readline()
    L11 = f.readline()
    DDL = L11.split('\t')

    #remplissage du tableau avec les données
    data = np.zeros((nl,nc))
    lines = f.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        for i in range (0,nc):
            data[l,i] = float(ligne[i])
        l += 1
    f.close()

    #calcul fréquence d'échantillonnage
    fs = 1/(data[2,0]-data[1,0])
    print("freq echantillonnage",fs)
    
    #filtrage des colonnes du tableau
    filteredData = np.zeros((nl,nc))
    filteredData[:,0] = data[:,0] #on remplit avec le temps NON FILTRE
    for i in range (1,nc): #on parcourt les colonnes sans la colonne temps
        filteredData[:,i] = filtrage(data[:,i],order,fc,fs)

    #création d'un fichier de résultat avec les données IK filtrées
    getcontext().prec = 8 #nombre de chiffres significatifs 
    fileName = IKresult.split('.')[0] + '_filtered.mot'
    path = os.path.join('IK_results', fileName)
    
    f = open(path,'w')
    f.write(L1),f.write(L2),f.write(L3),f.write(L4),f.write(L5),f.write(L6),f.write(L7),f.write(L8),f.write(L9),f.write(L10),f.write(L11)

    for i in range (0,nl): #on parcourt les lignes de données sans le header (qui est déjà écrit)
        #f.write('\t')
        for j in range (0,nc):
            f.write(str(round(filteredData[i,j],8))),f.write('\t')
            #f.write(str(Decimal(filteredData[i,j]))),f.write('\t')
        f.write('\n')
    f.close()
    
    return DDL,filteredData,data

########################################################################################
############################## Vitesse angulaire #######################################
########################################################################################
def vitesse_angulaire(IKresultsFilt): #Ajout 2023
    """calcul les vitesses angulaires à partir d'un fichier IK filtré.mot
input : fichier IK filtré "XXX_ik_filtered.mot"
output: crée un fichier ".mot" dans le dossier "IK_results" qui contient les vitesses angulaires en rad/s
        liste DDL, array à 1 colonne pour le temps, array vitesses angulaires : première colonne nulle correspondà la colonne temps
"""
    path = os.path.join('IK_results',IKresultsFilt)
    
    f = open(path,'r')

    L1,L2,L3,L4 = f.readline(),f.readline(),f.readline(),f.readline()
    nl = int(L3.split('=')[1])
    nc = int(L4.split('=')[1])
    L5,L6,L7,L8,L9,L10 = f.readline(),f.readline(),f.readline(),f.readline(),f.readline(),f.readline()
    L11 = f.readline()
    DDL = L11.split('\t')


    #remplissage du tableau avec les données
    data = np.zeros((nl,nc))
    lines = f.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        for i in range (0,nc):
            data[l,i] = float(ligne[i])
        l += 1
    f.close()

    ## paramètres de la dérivation
    time = data [:,0]
    dt = time[1] - time[0]

    ## calcul vitesse angulaire
    vitesse = np.zeros((nl,nc))
    vitesse[:,0] = time #on remplit la première colonne avec le temps
    for i in range (1,nc): #parcourt les colonnes pour les dériver
        vitesse[:,i] = np.gradient((data[:,i]*np.pi)/180,dt) #conversion en radian avant de dériver

        
    #création d'un fichier de résultat avec les données IK filtrées
    getcontext().prec = 8 #nombre de chiffres significatifs 
    fileName = IKresultsFilt.split('.')[0] + '_speed.mot'
    path = os.path.join('IK_results', fileName)
    
    f = open(path,'w')
    f.write(L1),f.write(L2),f.write(L3),f.write(L4),f.write(L5),f.write(L6),f.write(L7),f.write("Values are in rad per second\n"),f.write(L9),f.write(L10),f.write(L11)

    for i in range (0,nl): #on parcourt les lignes de données sans le header (qui est déjà écrit)
        #f.write('\t')
        for j in range (0,nc):
            f.write(str(round(vitesse[i,j],8))),f.write('\t')
            #f.write(str(Decimal(filteredData[i,j]))),f.write('\t')
        f.write('\n')
    f.close()
    
##    ## calcul accélération angulaire
##    acc = np.zeros((nl,nc))
##    for i in range (1,nc): #parcourt les colonnes pour les dériver
##      acc[:,i] = np.gradient(vitesse[:,i],dt)

    return DDL,vitesse,data



########################################################################################
############################## Puissance  ############################################## ajout 2023
########################################################################################
def puissance(fileID,fileSpeed): # ajout 2023
    """fonction qui calcule la puissance (scalaire) développée au niveau de l'articulation GH. Si les fichiers vitesse angulaire et ID n'ont pas la même longueur tq ID plus long de N lignes,
alors on ne prend pas en compte les N dernières lignes de l'iD. ça fait un petit décalage temporel. 
input : fileSpeed = fichier des vitesses angulaires 'XXX_ik_filtered_speed.mot". Doit être en rad/s
        fileID = fichiers des moments articulaires issus du calcul de l'ID 'XXX_id.sto"
output : créer un fichier .sto de résultats puissance dans le dossier "Power_results"
        dataPow = array des valeurs de puissance. Array écrit dans le fichier .sto de sortie"""

    ##import des données moment de l'articulation GH du fichier ID
    path = os.path.join('ID_results',fileID)

    fID = open(path,'r')

    L1,L2,L3,L4 = fID.readline(),fID.readline(),fID.readline(),fID.readline()
    nl = int(L3.split('=')[1])
    nc = int(L4.split('=')[1])
    L5,L6,L7 = fID.readline(),fID.readline(),fID.readline()
    DDL = L7.split('\t')

    #remplissage du tableau avec les données de temps et DDL de GH
    dataID = np.zeros((nl,4))
    lines = fID.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        dataID[l,0] = float(ligne[0]) #time
        dataID[l,1] = float(ligne[13]) #plane_elv
        dataID[l,2] = float(ligne[14])#shoulde_elv
        dataID[l,3] = float(ligne[15]) # axial_rot
        l += 1
        
    fID.close()
    
    ##import des données de vitesses angulaires de l'articultion GH du fichier speed
    path = os.path.join('IK_results',fileSpeed)

    fSp = open(path,'r')
    fSp.readline(),fSp.readline()
    nls,ncs = int(fSp.readline().split('=')[1]),int(fSp.readline().split('=')[1])
    fSp.readline(),fSp.readline(),fSp.readline(),fSp.readline(),fSp.readline(),fSp.readline(),fSp.readline()

    #remplissage d'un premier tableau avec toutes ls données de GH
    dataSp2 = np.zeros((nls,4))
    lines = fSp.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        dataSp2[l,0] = float(ligne[0]) #time
        dataSp2[l,1] = float(ligne[13]) #plane_elv
        dataSp2[l,2] = float(ligne[14])#shoulde_elv
        dataSp2[l,3] = float(ligne[15]) # axial_rot
        l += 1
        
    fSp.close()    

    ### Si les fichiers de vitesse angulaire et speed ont la même longueur
    if nls == nl : ##vérification que les fichiers ont la même longueur   
        ##calcul de la puissance
        dataPow = np.zeros((nl,2))
        for i in range (0,nl):
            dataPow[i,0] = dataSp2[i,0] #time
            dataPow[i,1] = dataID[i,1]*dataSp2[i,1] + dataID[i,2]*dataSp2[i,2] + dataID[i,3]*dataSp2[i,3] #puissance = produit scalaire du vecteur vitesse angulaire Gh et vecteur moment GH
            
        ##création d'un fichier de résultat avec les valeurs de la puissance articulaire
        getcontext().prec = 8 #nombre de chiffres significatifs 
        fileName = fileID.split('.')[0] + '_power.sto'
        path = os.path.join('Power_results', fileName)
        
        f = open(path,'w')
        f.write('Power glenohumeral joint - joint speed was in °/s for the computation of the power \n'),f.write(L2),f.write(L3),f.write(L4),f.write(L5),f.write(L6),f.write("time \t Power GH \n")

        for i in range (0,nl): #on parcourt les lignes de données sans le header (qui est déjà écrit)
            #f.write('\t')
            for j in range (0,2):
                f.write(str(round(dataPow[i,j],8))),f.write('\t')
            f.write('\n')
        f.close()

    ### Si les fichiers de vitesse angulaire et speed n'ont pas la même longueur
    else :
        print("les fichiers {} et {} n'ont pas le même nombre de lignes".format(fileID, fileSpeed))
        print(nls,nl)
        ##calcul de la puissance
        dataPow = np.zeros((nls,2))
        for i in range (0,nls):
            dataPow[i,0] = dataSp2[i,0] #time
            dataPow[i,1] = dataID[i,1]*dataSp2[i,1] + dataID[i,2]*dataSp2[i,2] + dataID[i,3]*dataSp2[i,3] #puissance = produit scalaire du vecteur vitesse angulaire Gh et vecteur moment GH
            
        ##création d'un fichier de résultat avec les valeurs de la puissance articulaire
        getcontext().prec = 8 #nombre de chiffres significatifs 
        fileName = fileID.split('.')[0] + '_power.sto'
        path = os.path.join('Power_results', fileName)
        
        f = open(path,'w')
        f.write('Power glenohumeral joint - joint speed was in °/s for the computation of the power \n'),f.write(L2),f.write(L3),f.write(L4),f.write(L5),f.write(L6),f.write("time \t Power GH \n")

        for i in range (0,nls): #on parcourt les lignes de données sans le header (qui est déjà écrit)
            #f.write('\t')
            for j in range (0,2):
                f.write(str(round(dataPow[i,j],8))),f.write('\t')
            f.write('\n')
        f.close()        
    
    return



########################################################################################
############################## Travail ################################################# #ajout 2023
########################################################################################
def travail(filePower): #ajout 2023
    """calcul le travail pour un fichier puissance. Equation 16 de Leboeuf et al., 2007 --> intégrale de la valeur abs
de la puissance pour prendre en compte les efforts moteurs et de freinage qui ont tous 2 un coût métabolique. 
input : filePower = fichier de puissance .sto
output : retourne la valeur du travail
"""
    ##import des données moment de l'articulation GH du fichier ID
    path = os.path.join('Power_results',filePower)

    f = open(path,'r')

    f.readline(),f.readline()
    nl = int(f.readline().split('=')[1])
    f.readline(),f.readline(),f.readline(),f.readline()

    #remplissage du tableau avec les données
    data = np.zeros((nl,2))
    lines = f.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        data[l,0] = float(ligne[0]) #time
        data[l,1] = float(ligne[1]) #puissance
        l += 1
        
    f.close()

    ##on enlève les valeurs de puissance nulle
    c=0 #compteur pour calculer l'index/ligne à partir duquel il y a une puissance nulle
    while c < nl and abs(data[c,1]) > 0.001 : 
        c+=1
    #print("c",c,"nl",nl)
    
    ##calcul du travail ##on va jusqu'à la ligne c car après c'est nul
    res = np.trapz(abs(data[0:c,1]),data[0:c,0]) #valeur absolue de la puissance
    
    return res

def travail_file(power_results):
    """fonction qui crée un fichier qui contient les données de trvail
pour tous les modèles et tous les cycle
input = liste. Liste contenant les noms des fichiers de puissance dont il faut calculer le travail.
        template du nom des éléments de la liste : 'model_000_E09C0_id_power.sto'
output = crée un fichier dans Power_results. Retourne un array des données de puissance"""
    
    
    #déterminer cb de lignes et de colonnes dans le tableau de stockage : les modèles en ligne et cycle en colonne
    list_model = []
    list_cycle = []
    for file in power_results:
        model= file.split('_')[1]
        cycle = file.split('_')[2]
        if model not in list_model:
            list_model.append(model)
        if cycle not in list_cycle:
            list_cycle.append(cycle)

    ##création du tableau qui contiendra les données de travail
    data = np.zeros((len(list_model),len(list_cycle)))

    #remplissage du tableau
    for m in range(0,len(list_model)):
        for c in range(0,len(list_cycle)):
            for r in power_results:
                if r.split('_')[1] == list_model[m] and r.split('_')[2] == list_cycle[c]:
                    #print("m",m,list_model[m])
                    #print("c",c,list_cycle[c])
                    #print("r",r)
                    data[m,c] = travail(r)
                    #print(travail(r))

    ##création d'un fichier qui contiendra les donnees de travail
    path = os.path.join('Power_results', 'GH joint work.txt')
    
    f = open(path,'w')
    
    #écriture du header
    f.write("model \t")
    for c in list_cycle:
        f.write(c),f.write("\t")
    f.write("\n")

    #écriture des données numériques
    for m in range (0,len(list_model)): #on parcourt les lignes de données sans le header (qui est déjà écrit)
        f.write(list_model[m]),f.write("\t")
        for c in range (0,len(list_cycle)):
            f.write(str(data[m,c])),f.write('\t')
        f.write('\n')
    f.close()
                    
    return data



##def travail_format_long(filetravail):
##    """importe un fichier de travail "GH_work.txt" en format large
##et retourne un fichier avec le tableau en format long"""
##    
##    ##import des données moment de l'articulation GH du fichier ID
##    path = os.path.join('Power_results',filetravail)
##
##    f = open(path,'r')
##
##    header = f.readline().split('\t')
##
##    data = np.zeros((27,12)) #nb de model,nb de cycle
##    models = []
##
##    lines = f.readlines()
##    l = 0
##    for line in lines:
##        ligne = line.split("\t")
##        model.append((ligne[0])) #model
##        for i in range(1,len(ligne)-1) : #on enlève la fin de la ligne, qui est juste un "\n"
##            data[l,i-1] = float(ligne[i]) #puissance
##        l += 1
##        
##    f.close()
##
##    ##met les données dans le bon sens : 1 colonne model, 1colonne cycle, 1 colonne valeur
##    dataReshape = np.zeros((27*12,3))
##    m=0
##    c=0
##    for model in models[0,-1]:
##        print(model)
##        for cycle in header[1,:]:
##            print(cycle)
##            dataReshape[,0]=
##            c+=1
##        m+=1        
##            
##    return

########################################################################################
############################## Valeurs temps ###########################################
########################################################################################
def valeurs_temps(file): #MAJ 2023
    """attribue les valeurs temporelles des cycles en fonction du fichier cycle_time.txt
input : fichier cycle_time.txt
output : nc9, nc11, nc13, nc19 = nb de cycles sélectionnés dans chaque essaie ;
        dataTime = array à 2 colonnes et nc9+nc11+nc13+nc19 lignes avec les valeurs de début et fin pour chaque cycle"""

    ftime = open(file,'r')
    ftime.readline()
    nc9 = int(ftime.readline().split(' ')[-1])
    nc11 = int(ftime.readline().split(' ')[-1])
    nc13 = int(ftime.readline().split(' ')[-1])
    nc19 = int(ftime.readline().split(' ')[-1])

    lines = ftime.readlines()
    n = len(lines)-4 #on enlève le nombre d'entêtes (### Essa..###) car on nevuet pas les mettre dans le tableau
    dataTime = np.zeros((n,2))

    i=0
    for line in lines :
        if list(line.split('\t')[0])[0] == '#':
            i=i
        else:
            dataTime[i,0]=float(line.split('\t')[1])
            dataTime[i,1]=float(line.split('\t')[2])
            i+=1
            
    return nc9, nc11, nc13, nc19, dataTime


def selection_cycle(file): #MAJ 2020
    """définit les cycles qui seront impliqués dans le calcul CMC. crée les listes
avec le nom des cycles
input = fichier cycle_time.txt
output = les 2 listes avec le numéro des cycles"""
    ftime = open(file,'r')

    ftime.readline()
    ftime.readline()
    ftime.readline()
    ftime.readline()
    
    lines = ftime.readlines()
    liste5 = []
    liste9 = []
    c = 5
    for line in lines:
        if list(line.split('\t')[0])[0] == '#' and c == 5:
            c = 9
        if list(line.split('\t')[0])[0:3] == ['c','y','c'] and c == 5:
            liste5.append(int(line.split('\t')[0].split('cycle')[1]))
        if list(line.split('\t')[0])[0:3] == ['c','y','c'] and c == 9:
            liste9.append(int(line.split('\t')[0].split('cycle')[1]))

    ftime.close()
        
    return liste5,liste9


########################################################################################
############################## Blocage DDL thorax ######################################
########################################################################################
##def lockDDLthorax(model): #MAJ2020
    """bloque les DDL en rotX et rotZ de l'articulation ground/thorax
input : fichier model.osim
output : nouveau fichier .osim avec les DDL bloqués"""
    Smodel = osim.Model(model)
    state = Smodel.initSystem()
    Smodel.initSystem()

    #accès aux DDL de l'articulation
    joint = Smodel.getJointSet().get('ground_thorax')
    rotX = joint.get_coordinates(0)
    rotZ = joint.get_coordinates(2)

    #blocage
    rotX.set_locked(True)
    rotZ.set_locked(True)

    Smodel.initSystem()
    
    return Smodel.printToXML(model)



###################################################################################################################
############################## Compute Muscle Moment Arm and Muscle moment potential############################### #Ajout 2023
###################################################################################################################
def computeMMA(Model,genericSetupAnlz,IKresult,resultsFolder,init_time,end_time) : 
    """Calcule les BDL musculaires au cours d'une cinématique réelle
input : Model : entrer le nom du fichier modèle .osim
        IkResult : entrer le nom du fichier .mot qui contient la cinématique
        genericSetupAnlz : entrer le nom du fichier .xml qui contient le setup gnl. La fonction entre seule les détails comme temps initial et final.
        resultsFolder : entrer le nom du dossier dans lequel on veut stocker les résultats (.mot et .sto)
output : moment arms files"""
    AnlzTool = osim.AnalyzeTool(genericSetupAnlz)

    IkFile = os.path.join('IK_results',IKresult)
    
    #load model
    model = osim.Model(os.path.join('scaled_model',Model))
    state = model.initSystem()
    model.initSystem()

    #set Analyze tool name
    AnlzTool.setName(IKresult.split('_')[0]+'_'+IKresult.split('_')[1]+'_'+IKresult.split('_')[2])
    
    #set AnalyzeTool results in a results folder
    AnlzTool.setResultsDir(resultsFolder)

    #tell AnalyzeTool tool to use the model
    AnlzTool.setModel(model)

    #setup the IK tool
    AnlzTool.setInitialTime(init_time)
    AnlzTool.setStartTime(init_time)
    AnlzTool.setFinalTime(end_time)
    #AnlzTool.setEndTime(end_time)
    AnlzTool.setCoordinatesFileName(IkFile)

    #run IK
    return AnlzTool.run()

def lectureAnalysis(file):
    path = os.path.join('MuscleAnalysis_results',file)
    
    f = open(path,'r')
    L1,L2,L3,L4 = f.readline(),f.readline(),f.readline(),f.readline()
    nl = int(L3.split('=')[1])
    nc = int(L4.split('=')[1])
    L5,L6,L7,L8,L9,L10,L11 = f.readline(),f.readline(),f.readline(),f.readline(),f.readline(),f.readline(),f.readline()
    L12 = f.readline()
    Muscles = L12.split('\t')

    #remplissage du tableau avec les données
    data = np.zeros((nl,nc))
    lines = f.readlines()
    l = 0
    for line in lines:
        ligne = line.split("\t")
        for i in range (0,nc):
            data[l,i] = float(ligne[i])
        l += 1
    f.close()
    
    return Muscles, data

def computeMMP(model,cycle,MuscleListe, MuscleListeUnite):
    """calcul le muscle moemnt potential. L'utilisateur choisit 1 model, cycle et 1 liste de muscle (la même que celle de analyze setting).
La fonction importe les 3 fichiers de DDL pour calcul le potentiel de moment résultant
input : model = string avec les 3 chiffres du model (exemple "000", "222").
        cycle = string avec l'essai et le cycle (exemple "E09C0")
        MuscleListe = liste des muscles sélectionnés pour étude
        MuscleListeUnite = liste des unités musculaires sélectionnées, suivant la liste de muscle donnée. 
output : écrit un fichier résultat pour 1 modele, 1cycle contenant le muscle moment potential au cours du temps
        retourne le tableau résultat de muscle moment potential"""

    #nombre de muscles étudiés (pas le nombre d'unités musculaires)
    N = len(MuscleListe)


    ###import modele ref pour extraire les Fisomax des unités muscualires étudiées
    ###############################################################################
    Flist = [] #tableau qui contient les noms des unités musculaires et la force iso max associée
    modelRef = osim.Model(os.path.join('model','Seth41_modelRef.osim'))
    musclesU = modelRef.getMuscles()
    NmusclesU = musclesU.getSize()
    for i in range(0,NmusclesU): #boucle qui parcourt les unités musculaires, i = entier
        MuscleU = musclesU.get(i)
        MuscleNameU = MuscleU.getName()
        for MuscleName in MuscleListe:
            if MuscleName in MuscleNameU:
                Flist.append(MuscleU.get_max_isometric_force()) #on stocke les valeurs des Fiso max dans l'ordre d'apparition du forceSet du model de référence. 
    #print(Flist)
                
    ###import des 3 fichiers de MMA (1 pour chaque DDL) #import des données brutes
    ###############################################################################
    #attribution du nom du fichier pour les 3 fichiers de MMA pour chaque DDL, nom des fichiers calculés avec la fonction AnalyzeTool
    Ddl_axial_rot = "model_"+model+"_"+cycle+"_MuscleAnalysis_MomentArm_axial_rot.sto"
    Ddl_plane_elv = "model_"+model+"_"+cycle+"_MuscleAnalysis_MomentArm_plane_elv.sto"
    Ddl_shoulder_elv = "model_"+model+"_"+cycle+"_MuscleAnalysis_MomentArm_shoulder_elv.sto"

    #appel des tableaux de données nommés au-dessus, calculés avec AnalyzeTool
    dataAxRot = lectureAnalysis(Ddl_axial_rot)[1] ##tableaux remplis des données de MMA
    dataPlaneElv = lectureAnalysis(Ddl_plane_elv)[1]
    dataShouldElv = lectureAnalysis(Ddl_shoulder_elv)[1]
    #print(dataAxRot)
    #print(dataPlaneElv)
    #print(dataShouldElv)
    
    #on fait une liste avec les tableaux de données brutes de MMA
    datalist = [dataAxRot,dataPlaneElv,dataShouldElv]
    #print(datalist)
    
    time = dataAxRot[:,0] #on définit la liste temps qui est commune à tous les tableaux de données car on travaille sur un seul cycle/descente
    #nombre de lignes des fichiers importés (ils ont tous la meme longueur car on travaille sur un seul cycle)
    nl = len(dataAxRot[:,0])

    
    ###calcul du MMP par DDL : Fisomax(unité 1)*BDL(unité 1) + Fisomax(unité 2)*BDL(unité 2)+ ...
    ##############################################################################################
    #création des tableaux vides qui acceuilleront les données de MMP pour chaque degré de liberté
    MMP_DDLAxRot = np.zeros((nl,N+1)) #tableau vide avec nl lignes (=longueur de la liste temps) et N+1 colonnes = nb de muscle + 1 colonne temps (N est bien associé au nb de muscles, pas d'unités musculaires)
    MMP_DDLPlaneElv = np.zeros((nl,N+1))
    MMP_DDLShouldElv = np.zeros((nl,N+1))
    
    #remplissage des premières colonnes avec les valeurs de temps
    MMP_DDLAxRot[:,0] = time
    MMP_DDLPlaneElv[:,0] = time
    MMP_DDLShouldElv[:,0] = time

    #création d'une liste qui contient les tableaux de données
    MMPDDL_list = [MMP_DDLAxRot,MMP_DDLPlaneElv,MMP_DDLShouldElv]

    #remplissage des tableaux (MMP_DDL..) qui acceuillent les moment potential DDL/DDL
    for jDDL in range (0,3): #jDDL = entier qui parcourt les DDL, si jDDL = 0, alors on travaille avec axial rot ; jDDL = 1, plane elv ; jDDL = 2, ShouldElv
        data = datalist[jDDL] #on sélectionne le bon tableau de données brutes de MMA associé au DDL jDDL
        MMPDDL = MMPDDL_list[jDDL] #on sélectionne le bon tableau vide qui acceuille les valeurs de MMP
        
        ##calculer moment potential d'1 muscle avec ses unités musculaires et Fiso max. 
        for muscle in MuscleListe: #parcourt des muscles
            iM = MuscleListe.index(muscle)
            for Umuscle in MuscleListeUnite: #parcourt des unités musculaires
                iMU = MuscleListeUnite.index(Umuscle)
                if muscle in Umuscle: #si l'unité musculaire fait partie du muscle, alors on l'intègre dans le calcul du moment potential du muscle 
                    Fisomax = Flist[iMU]
                    #print("jddl",jDDL, "muscle ",muscle)
                    #print("unite muscu", Umuscle, "index",iMU)
                    #print("fiso max",Fisomax)
                    #print("data",data[:,iMU+1])
                    for i in range(0,nl): #parcourt les lignes
                        MMPDDL[i,iM+1] = MMPDDL[i,iM+1] + Fisomax*data[i,iMU+1]
                    #print("MMPDDL",MMPDDL)
                    

    ### Calcul du moment potential final, avec les tableaux de chaque DDL : sqrt(MMP/DDL1**2 + MMP/DDL2**2 + MMP/DDL3**2)
    ########################################################################
    #création du tableau vide qui acceuillera les MMP
    MMP = np.zeros((nl,N+1))
    MMP[:,0] = time
    #remplissage du tableau
    for j in range (1,N+1): #comence à 1 car colonne temps déjà rempli, j représente muscle
        for i in range(0, nl): 
            MMP[i,j] = np.sqrt( MMP_DDLAxRot[i,j]**2 + MMP_DDLPlaneElv[i,j]**2 +MMP_DDLShouldElv[i,j]**2 )

    ##création d'un fichier qui contiendra les donnees de muscle moment potential
    path = os.path.join('MuscleAnalysis_results', 'model_{}_{}_Muscle Moment Potential.txt'.format(model, cycle))
    
    f = open(path,'w')
    
    #écriture du header
    f.write("Muscle Moment Potential"),f.write("\n")
    f.write("nRows: {}".format(nl)),f.write("\n")
    f.write("nColomns: {}".format(N+1)),f.write("\n")    
    f.write("Units: N.m"),f.write("\n")
    f.write("time \t")
    for c in MuscleListe:
        f.write(c),f.write("\t")
    f.write("\n")

    #écriture des données numériques
    for m in range (0,nl):
        for c in range (0,N+1):
            f.write(str(MMP[m,c])),f.write('\t')
        f.write('\n')
    f.close()
            
    return MMP

##########################################################################################
################################ Validation ##############################################
##########################################################################################
##osim.LoadOpenSimLibrary('C:\OpenSim 4.1\plugins\ScapulothoracicJointPlugin40_WinX64')
##model = 'model_000_E09_markerAdjusted.osim'
####lockDDLthorax(model)
####

####Filter
##data = 'model_000_E09_ik.mot'
##filt2 = filtrageIK(data,3,6,300)
##filt3 = filtrageIK(data,3,8,300)
##filt4 = filtrageIK(data,3,10,300)
##filt5 = filtrageIK(data,3,12,300)
##filt6 = filtrageIK(data,3,15,300)
##DDL = filt2[0]
##time = filt2[1][:,0]
##for ddl in range(1,len(DDL)):
##    plt.plot(time, filt2[2][:,ddl], 'k', label='raw_{}'.format(DDL[ddl]))
##    plt.plot(time, filt2[1][:,ddl], 'y--', label='filt Fc6_{}'.format(DDL[ddl]))
##    plt.plot(time, filt3[1][:,ddl], 'g--', label='filt Fc8_{}'.format(DDL[ddl]))
##    plt.plot(time, filt4[1][:,ddl], 'b--', label='filt Fc10_{}'.format(DDL[ddl]))
##    plt.plot(time, filt5[1][:,ddl], 'm--', label='filt Fc12_{}'.format(DDL[ddl]))
##    plt.plot(time, filt6[1][:,ddl], 'r--', label='filt Fc15_{}'.format(DDL[ddl]))
##    plt.legend()
##    plt.show()

#####vitesses angulaires
##data = 'model_000_E09_ik_filtered.mot'
##vit = vitesse_angulaire(data)
##DDL = vit[0]
##time = vit[1][:,0]
##for ddl in range(1,len(DDL)):
##    plt.plot(time, vit[2][:,ddl],'k',label="IK filtré")
##    plt.plot(time, vit[1][:,ddl], 'b--', label='vitesse ang_{}'.format(DDL[ddl]))
##    plt.legend()
##    plt.grid()
##    plt.show()
##    
####Puissance
##dataID = 'model_000_E09C0_id.sto'
##dataSp = 'model_000_E09C0_ik_filtered_speed.mot'
##puissance(dataID,dataSp)

##Travail
##filePower='model_000_E09C0_id_power.sto'
##print(travail(filePower))

###compute MMA
#genericSetupAnlz = "Analyze_setting.xml"
##IKresult = "model_000_E09C0_ik_filtered.mot"
##resultsFolder = "MuscleAnalysis_results"
##init_time = 14.30
##end_time = 14.44
##
#computeMMA(model,genericSetupAnlz,IKresult,resultsFolder,init_time,end_time) 
#MuscleListe = ["PectoralisMajor","TeresMajor","Infraspinatus","Subscapularis"]
#MuscleListeUnite = ["PectoralisMajorClavicle_S","PectoralisMajorThorax_I","PectoralisMajorThorax_M","TeresMajor","Infraspinatus_I","Infraspinatus_S","Subscapularis_S","Subscapularis_M","Subscapularis_I"]
#computeMMP("020","E13C8",MuscleListe,MuscleListeUnite)
