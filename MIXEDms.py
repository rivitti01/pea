import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def MIXEDmsSolve(NStations, c, NcClasses, cD, cN, cZ, NoClasses, oD, ol):
    #print(NStations, c, NcClasses, cD, cN, cZ, NoClasses, oD, ol)

    ###### Compute the utilizations for the open part, and inflate the demands of the closed classes

    oUk  = np.zeros(NStations)

    for k in range(0, NStations):
        U = 0
        for ic in range(0, NoClasses):
            U = U + ol[ic] * oD[k, ic]
        if(c[k] > 0):
            oUk[k] = U / c[k]
        else:
            oUk[k] = U

    cDinf = np.zeros((NStations, NcClasses))
    for k in range(0, NStations):
        for ic in range(0, NcClasses):
            if(c[k] > 0):           # inflating is not meaninful for delay stations
                cDinf[k, ic] = cD[k, ic] / (1 - oUk[k])
            else:
                cDinf[k, ic] = cD[k, ic]
    
    #print(cDinf)
    
    ###### Solves the closed part

    maxRemJob = np.zeros(NcClasses)
    for i in range(0, NcClasses-1):
        maxRemJob[NcClasses-2-i] = maxRemJob[NcClasses-1-i] + cN[NcClasses-i-1]

    curConfName = "0"
    for i in range(1, NcClasses):
        curConfName = curConfName + "_0"
    cPki = np.zeros((NStations, 1))
    for k in range(0, NStations):
        cPki[k, 0] = 1
    olddict = {curConfName: {'cnt':0, 'cNk': np.zeros(NStations), 'cPki': cPki}}    
    #print(olddict)

    Ntot = int(sum(cN))
    for nt in range(1, Ntot+1):
    #    print("----> ", nt)
        
        Nit = np.zeros(NcClasses)
        srem = np.zeros(NcClasses)
        smaxJ = np.zeros(NcClasses)
        
        cit = 0
        nomore = 0
        rem = nt
        cnt = 0
        
        confdict = {}

        while (nomore == 0):
            # create the current conifguration
            while (cit < NcClasses-1):     # saturate the stations with the remaining jobs
                minJ = max(0, rem - maxRemJob[cit])
                maxJ = min(rem, cN[cit])
                rem = rem - minJ
                Nit[cit] = minJ
                
                smaxJ[cit] = maxJ
                srem [cit] = rem

                cit = cit + 1
            Nit[cit] = rem

            #solve configuration
            cRkc = np.zeros((NStations, NcClasses))
            cNkc = np.zeros((NStations, NcClasses))
            cXc = np.zeros(NcClasses)
            cRc = np.zeros(NcClasses)
            cNk = np.zeros(NStations)
            cPki = np.zeros((NStations, nt+1))
            
            for ic in range(0, NcClasses):
                #consider class c
                if (Nit[ic] > 0):    #class must have at least one job
                    # retrieves the solution of the prev. model
                    Npc = Nit.copy()
                    Npc[ic] = Npc[ic] -1 
                    pcName = str(int(Npc[0]))
                    for i in range(1, NcClasses):
                        pcName = pcName + "_" + str(int(Npc[i]))
                    PcNK  = olddict[pcName]['cNk']
                    PcPki = olddict[pcName]['cPki']
                    
                    for ik in range(0, NStations):
                        if (c[ik] == 0): # infinite server
                            cRkc[ik, ic] = cDinf[ik, ic]
                        elif (c[ik] >= 2.0):
                            cRkc[ik, ic] = 0
                            for j in range(1, nt+1):
                                cd = min(float(j), c[ik])
                                ### print(nt, j, ic, ik, Rkc[ik, ic], cd, float(j) / cd, PcPki[ik,j-1])
                                cRkc[ik, ic] = cRkc[ik, ic] + float(j) / cd * PcPki[ik, j-1]
                            cRkc[ik, ic] = cRkc[ik, ic] * cDinf[ik, ic]
                        else:
                            cRkc[ik, ic] = cDinf[ik, ic] * (1 + PcNK[ik])
                        cRc[ic] = cRc[ic] + cRkc[ik, ic]
                    cXc[ic] = Nit[ic] / (cZ[ic] + cRc[ic])
                    for ik in range(0, NStations):
                        if (c[ik] >= 2.0):
                            for j in range(nt, 0, -1):
                                cd = min(float(j), c[ik])
                                cPki[ik, j] = cPki[ik, j] + cXc[ic] * cDinf[ik, ic] / cd * PcPki[ik, j-1]

            for ik in range(0, NStations):
                for ic in range(0, NcClasses):
                    cNkc[ik, ic] = cXc[ic] * cRkc[ik, ic]
                    cNk[ik] = cNk[ik] + cNkc[ik, ic]
                spk = 0
                for j in range(0, nt+1):
                    spk = spk + cPki[ik, j]
                cPki[ik, 0] = 1 - spk

            #store configuration
            curConfName = str(int(Nit[0]))
            for i in range(1, NcClasses):
                curConfName = curConfName + "_" + str(int(Nit[i]))
            #print(curConfName, cnt, Nit)
            
            confdict[curConfName] = {'id': cnt, 'cNk': cNk, 'cPki': cPki}
            cnt = cnt + 1

            # advance to the next configuration in this level, if available
            cit = NcClasses-2      # chekcs if some class has not saturated yet
            while( (cit >= 0) and (Nit[cit] >= smaxJ[cit])):
                cit = cit - 1

            if(cit >= 0):
                Nit[cit] = Nit[cit] + 1
                rem = srem[cit] - 1
                srem [cit] = rem
                cit = cit + 1
            else:
                nomore = 1
        olddict = confdict
    
    ###### Solves the open part
    oRkc  = np.zeros((NStations, NoClasses))
    for k in range(0, NStations):
        U = oUk[k]
        
        if (c[k] <= 0.0):
            RF = 1
            AvN = 1
        elif (c[k] >= 2.0):
            Fc = 1
            Sm = 1
            Tr = 1
            for j in range(1, int(c[k])):
                Tr = Tr * U / j
                Fc = Fc * j / U
                Sm = Sm + Tr
            RF = (1 + 1 / (c[k] - U) / (1 + (c[k] - U) / U * Fc * Sm))
            AvN = 0
            for j in range(0, Ntot+1):
                cd = min(float(j + 1), c[k])
                AvN = AvN + float(j + 1) / cd * cPki[ik, j]
        else:
            RF = 1 / (1 - U)
            AvN = 1 + cNk[k]
        #print(k, c[k], cNk[k], RF, AvN)
        for ic in range(0, NoClasses):
            oRkc[k, ic] = oD[k, ic] * RF * AvN

    return {'cXc':cXc, 'cRc':cRc, 'cNkc':cNkc, 'cRkc':cRkc, 'cNk':cNk, 'oUk':oUk, 'oRkc':oRkc}

