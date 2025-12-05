import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg

def MVAmcssSolve(NStations, NClasses, D, c, N, Z):

    maxRemJob = np.zeros(NClasses)
    for i in range(0, NClasses-1):
        maxRemJob[NClasses-2-i] = maxRemJob[NClasses-1-i] + N[NClasses-i-1]
        
    print(N)
    print(maxRemJob)

    curConfName = "0"
    for i in range(1, NClasses):
        curConfName = curConfName + "_0"
    olddict = {curConfName: {'cnt':0, 'Nk': np.zeros(NStations)}}
    #print(olddict)

    Ntot = int(sum(N))
    for nt in range(1, Ntot+1):
    #    print("----> ", nt)
        
        Nit = np.zeros(NClasses)
        srem = np.zeros(NClasses)
        smaxJ = np.zeros(NClasses)
        
        cit = 0
        nomore = 0
        rem = nt
        cnt = 0
        
        confdict = {}

        while (nomore == 0):
            # create the current conifguration
            while (cit < NClasses-1):     # saturate the stations with the remaining jobs
                minJ = max(0, rem - maxRemJob[cit])
                maxJ = min(rem, N[cit])
                rem = rem - minJ
                Nit[cit] = minJ
                
                smaxJ[cit] = maxJ
                srem [cit] = rem

                cit = cit + 1
            Nit[cit] = rem

            #solve configuration
            Rkc = np.zeros((NStations, NClasses))
            Nkc = np.zeros((NStations, NClasses))
            Xc = np.zeros(NClasses)
            Rc = np.zeros(NClasses)
            Nk = np.zeros(NStations)
            
            for ic in range(0, NClasses):
                #consider class c
                if (Nit[ic] > 0):    #class must have at least one job
                    # retrieves the solution of the prev. model
                    Npc = Nit.copy()
                    Npc[ic] = Npc[ic] -1 
                    pcName = str(int(Npc[0]))
                    for i in range(1, NClasses):
                        pcName = pcName + "_" + str(int(Npc[i]))
                    PcNK = olddict[pcName]['Nk']
                    
                    for ik in range(0, NStations):
                        if (c[ik] == 0): # infinite server
                            Rkc[ik, ic] = D[ik, ic]
                        else:
                            Rkc[ik, ic] = D[ik, ic] * (1 + PcNK[ik])
                        Rc[ic] = Rc[ic] + Rkc[ik, ic]
                    Xc[ic] = Nit[ic] / (Z[ic] + Rc[ic])
            for ik in range(0, NStations):
                for ic in range(0, NClasses):
                    Nkc[ik, ic] = Xc[ic] * Rkc[ik, ic]
                    Nk[ik] = Nk[ik] + Nkc[ik, ic]

            #store configuration
            curConfName = str(int(Nit[0]))
            for i in range(1, NClasses):
                curConfName = curConfName + "_" + str(int(Nit[i]))
            #print(curConfName, cnt, Nit)
            
            confdict[curConfName] = {'id': cnt, 'Nk': Nk}
            cnt = cnt + 1

            # advance to the next configuration in this level, if available
            cit = NClasses-2      # chekcs if some class has not saturated yet
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
    
    return {'Xc':Xc, 'Rc':Rc, 'Nkc':Nkc, 'Rkc':Rkc, 'Nk':Nk}


# The idea of the enumartion algorithm:
# Always find the minimum and maximum number of jobs for a class,
# when there are [rem] jobs to distribute
# Start assigning each class its minimum.
# Then, try increase each class the number of jobs of one unit
# until the maximum is reached, and redistribute the remaining jobs - 1
# in the following. Two additional data structures, rembering the minimum
# number and the remaining number of jobs are used for back-trakcing

print(MVAmcssSolve(3,4,np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]]),np.array([1,1,1]),np.array([5,6,7,8]),np.array([0,0,0,0])))