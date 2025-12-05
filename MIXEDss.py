import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def MIXEDssSolve(NStations, c, NcClasses, cD, cN, cZ, NoClasses, oD, ol):
    #print(NStations, c, NcClasses, cD, cN, cZ, NoClasses, oD, ol)

    ###### Compute the utilizations for the open part, and inflate the demands of the closed classes

    oUk  = np.zeros(NStations)

    for k in range(0, NStations):
        U = 0
        for ic in range(0, NoClasses):
            U = U + ol[ic] * oD[k, ic]
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
    olddict = {curConfName: {'cnt':0, 'cNk': np.zeros(NStations)}}
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
            
            for ic in range(0, NcClasses):
                #consider class c
                if (Nit[ic] > 0):    #class must have at least one job
                    # retrieves the solution of the prev. model
                    Npc = Nit.copy()
                    Npc[ic] = Npc[ic] -1 
                    pcName = str(int(Npc[0]))
                    for i in range(1, NcClasses):
                        pcName = pcName + "_" + str(int(Npc[i]))
                    PcNK = olddict[pcName]['cNk']
                    
                    for ik in range(0, NStations):
                        if (c[ik] == 0): # infinite server
                            cRkc[ik, ic] = cDinf[ik, ic]
                        else:
                            cRkc[ik, ic] = cDinf[ik, ic] * (1 + PcNK[ik])
                        cRc[ic] = cRc[ic] + cRkc[ik, ic]
                    cXc[ic] = Nit[ic] / (cZ[ic] + cRc[ic])
            for ik in range(0, NStations):
                for ic in range(0, NcClasses):
                    cNkc[ik, ic] = cXc[ic] * cRkc[ik, ic]
                    cNk[ik] = cNk[ik] + cNkc[ik, ic]

            #store configuration
            curConfName = str(int(Nit[0]))
            for i in range(1, NcClasses):
                curConfName = curConfName + "_" + str(int(Nit[i]))
            #print(curConfName, cnt, Nit)
            
            confdict[curConfName] = {'id': cnt, 'cNk': cNk}
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
        else:
            RF = 1 / (1 - U)

        for ic in range(0, NoClasses):
           if(c[k] > 0):
                oRkc[k, ic] = oD[k, ic] * RF * (1 + cNk[k])
           else:
                oRkc[k, ic] = oD[k, ic]

    return {'cXc':cXc, 'cRc':cRc, 'cNkc':cNkc, 'cRkc':cRkc, 'cNk':cNk, 'oUk':oUk, 'oRkc':oRkc}

