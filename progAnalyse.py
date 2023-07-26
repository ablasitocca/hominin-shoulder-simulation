#programme qui utilise la bibliothèque de fonctions bib.py pour lancer les calculs

from biblioAnalyse import *
import opensim as osim
import os
import numpy as np
from pathlib import Path
 

## Import plugin scapulothoracic joint
osim.LoadOpenSimLibrary('C:\OpenSim 4.1\plugins\ScapulothoracicJointPlugin40_WinX64')

## Import des valeurs temporelles
ftime ='cycle_time.txt' #'cycle_time_archive tous les cycles.txt'
nc9,nc11,nc13,nc19=valeurs_temps(ftime)[0],valeurs_temps(ftime)[1],valeurs_temps(ftime)[2],valeurs_temps(ftime)[3]
dataTime = valeurs_temps(ftime)[4]

############################# Adjusting markers #####################################MAJ 2023
##models = []
##files = os.listdir('model')
##for file in files:
##    if file.split('_')[0] == 'model':
##        models.append(file)
##
##genericSetupScale = 'scale_setting.xml'
##data = ['Homtech_20180913_09_staticR.trc','Homtech_20180913_11_static.trc','Homtech_20180913_13_static.trc','Homtech_20180914_19_static.trc']
##
##for modelName in models :
##    for static_trc_data in data :
##        E = static_trc_data.split('_')[2]
##        Model = os.path.join('model',modelName)
##        outputFileName = modelName.split('.')[0]
##
##        # lance la fonction d'ajustement des marqueurs
##        adjustMarker(Model,static_trc_data,genericSetupScale,outputFileName)
##
##        # lance une seconde fois l'ajustement des marqueurs sur le modèle ajusté 1 fois
##        modelName2 = modelName.split('.')[0] + '_E{}_markerAdjusted.osim'.format(E)
##        Model2 = os.path.join('scaled_model',modelName2)
##        adjustMarker(Model2,static_trc_data,genericSetupScale,outputFileName)
##print("Ajustement OK")


########################################## Inverse kinematics #############################MAJ 2023
##models = []
##files = os.listdir('scaled_model')
##for file in files:
##    if file.split('_')[0] == 'model':
##        models.append(file)
##
###ligne à enlever # pour faire que qqes models
###models = ["model_000_E09_markerAdjusted.osim","model_000_E11_markerAdjusted.osim","model_000_E13_markerAdjusted.osim","model_000_E19_markerAdjusted.osim","model_002_E09_markerAdjusted.osim","model_002_E11_markerAdjusted.osim","model_002_E13_markerAdjusted.osim","model_002_E19_markerAdjusted.osim","model_020_E09_markerAdjusted.osim","model_020_E11_markerAdjusted.osim","model_020_E13_markerAdjusted.osim","model_020_E19_markerAdjusted.osim","model_200_E09_markerAdjusted.osim","model_200_E11_markerAdjusted.osim","model_200_E13_markerAdjusted.osim","model_200_E19_markerAdjusted.osim","model_222_E09_markerAdjusted.osim","model_222_E11_markerAdjusted.osim","model_222_E13_markerAdjusted.osim","model_222_E19_markerAdjusted.osim"]
##
###generic setup file for IK
##genericSetupIK = 'IK_setting.xml'
##      
###results folder
##results_folder = 'IK_results'
##
###data
##data = ['Homtech_20180913_09_29sR.trc','Homtech_20180913_11.trc','Homtech_20180913_13.trc','Homtech_20180914_19_TrimmedRecomputed_3.trc']
##
##c = 0 #compter l'avancement    
##for model in models:
##    print("IK :",(c/len(models))*100,"%")
##    
##    if model.split('_')[2] == 'E09' :
##        #print("E09")
##        #kinemtaics data file 
##        trc_data = data[0]
##        
##        for i in range (0,nc9):
##            #initial and final time
##            init_time = dataTime[i,0]
##            end_time = dataTime[i,1]
##            print(init_time,end_time)
##            
##            #output prefix name
##            outputName = model.split('_marker')[0] + 'C{}'.format(i)
##            #print(outputName)
##            
##            # launch IK function
##            IKin(model,trc_data,genericSetupIK,results_folder,init_time,end_time,outputName)
##
##    elif model.split('_')[2] == 'E11' :
##        #print("E11")
##        #kinemtaics data file 
##        trc_data = data[1]
##
##        for i in range (nc9,nc9+nc11):
##            #initial and final time
##            init_time = dataTime[i,0]
##            end_time = dataTime[i,1]
##            print(init_time,end_time)
##        
##            #output prefix name
##            outputName = model.split('_marker')[0] + 'C{}'.format(i)
##            print(outputName)
##            
##            # launch IK function
##            IKin(model,trc_data,genericSetupIK,results_folder,init_time,end_time,outputName)
##
##    elif model.split('_')[2] == 'E13' :
##        #print("E13")
##        #kinemtaics data file 
##        trc_data = data[2]
##
##        for i in range (nc9+nc11,nc9+nc11+nc13):
##            #initial and final time
##            init_time = dataTime[i,0]
##            end_time = dataTime[i,1]
##            print(init_time,end_time)
##            
##            #output prefix name
##            outputName = model.split('_marker')[0] + 'C{}'.format(i)
##            print(outputName)
##            
##            # launch IK function
##            IKin(model,trc_data,genericSetupIK,results_folder,init_time,end_time,outputName)
##
##    elif model.split('_')[2] == 'E19' :
##        #print("E19")
##        #kinemtaics data file 
##        trc_data = data[3]
##
##        for i in range (nc9+nc11+nc13,nc9+nc11+nc13+nc19):
##            #initial and final time
##            init_time = dataTime[i,0]
##            end_time = dataTime[i,1]
##            #print(init_time,end_time)
##            
##            #output prefix name
##            outputName = model.split('_marker')[0] + 'C{}'.format(i)
##            #print(outputName)
##            
##            # launch IK function
##            IKin(model,trc_data,genericSetupIK,results_folder,init_time,end_time,outputName)
##    c+=1
##    
##print("IK OK")



############## Filtrage de IK ######################################MAJ 2023
########nettoyage du dossier avant de lancer l'étape de filtrage
##files = os.listdir('IK_results')
##for file in files:
##    if "filtered" in file.split('_') or "filtered.mot" in file.split('_'):
##        path = os.path.join('IK_results',file)
##        os.remove(path)
##        
##
####début de l'étape de filtrage
##results = []
##files = os.listdir('IK_results')
##for file in files:
##    if file.split('_')[0] == 'model':
##        results.append(file)
##
##results=["model_200_E09C2_ik.mot"] 
##
##order = 3 #ordre du filtre
##fc = 10 #fréquence de coupure
##
##for ik_result in results: 
##    filtrageIK(ik_result,order,fc)
##    
##print("filtrage IK OK")
##



############## Calcul vitesse angulaire ######################################Ajout 2023      
######etape de calcul des vitesses angulaires
##results = []
##files = os.listdir('IK_results')
##for file in files:
##    if file.split('_')[-1] == 'filtered.mot':
##        results.append(file)
##        
##
##for ik_result_filt in results: ## le fichier de vitesse angulaire retourné est en rad/s
##    vitesse_angulaire(ik_result_filt)
##    
##print("calcul vitesses angulaires OK")




################################ Inverse dynamics #######################################MAJ 2023
####ATTENTION, pour choisir quel type de filtrage on veut appliquer sur les données d'entrée (donc IK), choisir dans le code suivant la fréquence de coupure
######et choirsi si on prend ik.mot ou filtered.mot en données d'entrée. 
##models = []
##files = os.listdir('scaled_model')
##for file in files:
##    if file.split('_')[0] == 'model':
##        models.append(file)
##
##ik_results = []
##files = os.listdir('IK_results')
##for file in files:
##    if 'model' in file.split('_') and 'ik.mot' in file.split('_'): ##'ik.mot' au lieu de 'filtered.mot' pour lancer sur IK brute
##        ik_results.append(file)
##
###generic setup ID
##genericSetupID = 'ID_setting.xml'
###results folder
##results_folder = 'ID_results'
###frequence de coupure
##freq = 6 #application du filtre OS car on lance l'ID sur Ik brute // si pas de filtre : -1
##
##for ik_result in ik_results:
##    print("IK", ik_result)
##    for Model in models:
##        
##        if ik_result.split('_')[1] == Model.split('_')[1] and ik_result.split('_')[2].split('C')[0] == Model.split('_')[2] : #association de l'essai du model XXX avec le model XXX scalé pour l'essai
##            print("IK - Model",ik_result, Model)
##            #output prefix name
##            outputName = ik_result.split('_ik')[0]
##
##            #set initial time and final time for launching ID, depending on the cycle of ik_result
##            cycle = int(ik_result.split('_')[2].split('C')[1])
##            init_time = dataTime[cycle,0]
##            end_time = dataTime[cycle,1]
##            
##            IDyn(Model,ik_result,genericSetupID,results_folder,init_time,end_time,freq,outputName)
##
##print("ID OK")



############## Calcul Puissance ######################################Ajout 2023
#####création liste avec les fichiers de vitesse angulaire
##speed_results = []
##files = os.listdir('IK_results')
##for file in files:
##    if 'model' in file.split('_') and 'speed.mot' in file.split('_'):
##        speed_results.append(file)
##
##
###création liste avec les fichiers ID 
##ID_results = []
##files = os.listdir('ID_results')
##for file in files:
##    if 'model' in file.split('_') and 'id.sto' in file.split('_'):
##        ID_results.append(file)
##
##
##for fileID in ID_results:
##    for fileSpeed in speed_results:
##        if fileID.split('_')[1] == fileSpeed.split('_')[1] and fileID.split('_')[2] == fileSpeed.split('_')[2] : #association des fichiers ID et speed ensemble suivant leur nom 
##            puissance(fileID, fileSpeed)
##        
##print("calcul puissance OK")    

############## Calcul Travail ######################################Ajout 2023
power_results = []
files = os.listdir('Power_results')
for file in files:
    if 'model' in file.split('_') and 'power.sto' in file.split('_'):
        power_results.append(file)

travail_file(power_results)


###################### Calcul Muscle moment potential ######################################Ajout 2023
######Calcul d'abord des muscles moment arms
##############################################
##
##models = []
##files = os.listdir('scaled_model')
##for file in files:
##    if file.split('_')[0] == 'model':
##        models.append(file)
##        
##
##ik_results = []
##files = os.listdir('IK_results')
##for file in files:
##    if 'model' in file.split('_') and 'filtered.mot' in file.split('_'): ##'ik.mot' au lieu de 'filtered.mot' pour lancer sur IK brute
##        ik_results.append(file)
##
##
###generic setup Analyze tool
##genericSetupAnlz = 'Analyze_setting.xml'
###results folder
##results_folder = 'MuscleAnalysis_results'
##
##for ik_result in ik_results:
##    print("IK", ik_result)
##    for Model in models:
##        
##        if ik_result.split('_')[1] == Model.split('_')[1] and ik_result.split('_')[2].split('C')[0] == Model.split('_')[2] : #association de l'essai du model XXX avec le model XXX scalé pour l'essai
##            print("IK - Model",ik_result, Model)
##
##            #set initial time and final time for launching ID, depending on the cycle of ik_result
##            cycle = int(ik_result.split('_')[2].split('C')[1])
##            init_time = dataTime[cycle,0]
##            end_time = dataTime[cycle,1]
##            
##            computeMMA(Model,genericSetupAnlz,ik_result,results_folder,init_time,end_time)
##
####Nettoyage des résultats, on enlève les fichiers dont on n'a pas besoin. ATTENTION ! ne pas effacer le code d'affichage
##file_AnlzResults = os.listdir('MuscleAnalysis_results')
##for file in file_AnlzResults:
##    if 'MomentArm' not in file.split('_') and 'py' not in file.split('.'):
##        path = os.path.join('MuscleAnalysis_results',file)
##        os.remove(path)   
##
##print("Musle Moment Arm OK")
##


##########calcul muscle moment potential
#########################################
##models = ["000","001","002","010","020","100","200","111","222"]
##cycles = ["E09C0","E09C1","E09C2","E11C3","E11C4","E11C5","E13C6","E13C7","E13C8","E19C9","E19C10","E19C11"]
#####ATTENTION, les muscles doivent être dans l'ordre du model + il faut que la liste des unités musculaires soit la même que dans le fichier "analyse_setting.xml"
##MuscleListe = ["PectoralisMajor","TeresMajor","Infraspinatus","Subscapularis"]
##MuscleListeUnite = ["PectoralisMajorClavicle_S","PectoralisMajorThorax_I","PectoralisMajorThorax_M","TeresMajor","Infraspinatus_I","Infraspinatus_S","Subscapularis_S","Subscapularis_M","Subscapularis_I"]
##
##for m in models:
##    for c in cycles:
##        computeMMP(m,c,MuscleListe,MuscleListeUnite)
##
##
##print("Muscle Moment Potential OK")
