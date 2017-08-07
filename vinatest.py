import subprocess
import time
import commands
import string
import os

def configParameters(name, ligandname, proteinname):
    params = open("ligands-config.txt", "r")
    content = params.readlines()
    params.close()
    flag = 0
    data = []
    for line in content:
        lig = line.strip().split("_")
        lig[1] = lig[1].split(";")
        data.append(lig)
    for ligand in data:
        if name == ligand[0]: 
            center_x = ligand[1][0]
            center_y = ligand[1][1]
            center_z = ligand[1][2]
            flag = 1
            break
    if flag == 1:
        dataW = "receptor = protein/"+proteinname+"\n"
        dataW += "ligand = m-ligand/"+ligandname+"\n"
        dataW += "center_x = "+str(center_x)+"\n"
        dataW += "center_y = "+str(center_y)+"\n"
        dataW += "center_z = "+str(center_z)+"\n"
        dataW += "size_x = 20.0\nsize_y = 20.0\nsize_z = 20.0\n"
        dataW += "out = result.pdbqt"
        configFile = open("config.txt", "w")
        configFile.write(dataW)
        configFile.close()
    else:
        print "ERROR: no config parameters for ligand", name
        exit(2)

def calculateFreeEnergy():
    a = commands.getstatusoutput('./vina --config config.txt --log log-out.txt')
    energiaVina = string.split(string.split(a[1],"Affinity:")[1])[0]
    #print energiaVina
    return float(energiaVina)



name = "NAD"
ligand = "nad-modify.pdbqt"
protein = "1ENY.pdbqt"
configParameters(name, ligand, protein)
for i in range(10):
    print "Running Process " + str(i)
    p = subprocess.Popen('./vina --config config.txt --log log-out.txt --cpu 10', shell=True)
    p.communicate()
    time.sleep(10)
    os.rename("result.pdbqt", "run_"+str(i)+".pdbqt")
    os.rename("log-out.txt", "log-out_"+str(i)+".txt")
    


