import numpy as np

def read_WAVEDER(file='WAVEDER'):
        '''Read the WAVEDER file

        the matrix of the derivative of the cell periodic part
        of the wavefunctions with respect to i k is:      
        cder = CDER_BETWEEN_STATES(m,n,1:NK,1:ISP,j)= <u_m|  -i d/dk_j | u_n> = - <u_m| r_j |u_n>
        '''
        from scipy.io import FortranFile
        file = FortranFile(file, 'r')
        nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
        nodesn_i_dielectric_function = file.read_record(dtype= np.float)
        wplasmon = file.read_record(dtype= np.float).reshape(3,3)
        cder = file.read_record(dtype= np.complex64).reshape(3, ispin,nkpts,nbands_cder,nb_tot)


        return cder, nodesn_i_dielectric_function, wplasmon

def readnum(txt,i):
        a=''
        while txt[i].isdigit()==False and txt[i] not in ['-']:
                i=i+1
        while txt[i].isdigit()==True or txt[i] in ['-','+','e','E','.']:
                a=a+txt[i]
                i=i+1
        return float(a),i

def read_EIGENVAL(fileName='EIGENVAL'):
        with open(fileName,'r') as file:
                lines = [line.rstrip() for line in file]
        tmp=np.array([float(line) for line in lines[5].split()])
        nk, nb = int(tmp[1]), int(tmp[2])
        rowLen = np.array([float(line) for line in lines[8].split()]).shape[0]
        if rowLen==3:
                BS = np.zeros((nk,nb))
                KP = np.zeros((nk,3))
                for k in range(nk):
                        for b in range(nb):
                                BS[k, b] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb + (b + 1)].split()])[1]
                        KP[k, :] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + (k + 1) * nb].split()])[0:3]
        else:
                print('Strange rowLen')

        return KP,BS

def calcP(nb1_start,nb1_end,nb2_start,nb2_end,km,BS):
        P0x=0.
        P0y=0.
        P0z=0.
        for nb in range(nb1_start,nb1_end+1):
                for nb2 in range(nb2_start,nb2_end+1):
                        '''print(abs(km[0,0,0,nb,nb2]))
                        print(abs(km[1,0,0,nb,nb2]))
                        print(abs(km[2,0,0,nb,nb2]))
                        print(abs(np.array([km[0,0,0,nb,nb2],km[1,0,0,nb,nb2],km[2,0,0,nb,nb2]])))
                        print((abs(km[0,0,0,nb,nb2])+abs(km[1,0,0,nb,nb2])+abs(km[2,0,0,nb,nb2]))/3.)'''
                        P0x=P0x+(abs(km[0,0,0,nb,nb2])*(BS[0,nb]-BS[0,nb2]))**2
                        P0y=P0y+(abs(km[1,0,0,nb,nb2])*(BS[0,nb]-BS[0,nb2]))**2
                        P0z=P0z+(abs(km[2,0,0,nb,nb2])*(BS[0,nb]-BS[0,nb2]))**2
        return(P0x**0.5,P0y**0.5,P0z**0.5)

def read_POSCAR(fileName='POSCAR'):
        with open(fileName,'r') as file:
                lines = [line.rstrip() for line in file]
        acell=float(lines[1])
        a1=np.array([float(line) for line in lines[2].split()])
        a2=np.array([float(line) for line in lines[3].split()])
        a3=np.array([float(line) for line in lines[4].split()])
        natom=int(np.array([float(line) for line in lines[6].split()]).sum())
        coords=np.zeros((natom,3))
        for atom in range(natom):
                coords[atom]=np.array([float(line) for line in lines[8+atom].split()])
        return acell*a1, acell*a2, acell*a3, coords

def read_PROCAR(file='PROCAR',SO=False):
        f = open(file)
        txt = f.read()
        txt=txt+'!!!'
        f.close()
        it = 0
        for header in range(2):
                it = readnum(txt, it)[1]
        nk, it = readnum(txt, it)
        nb, it = readnum(txt, it)
        nion, it = readnum(txt, it)
        nk, nb, nion = int(nk), int(nb), int(nion)
        BS = np.zeros((nk, nb))
        KP = np.zeros((nk, 3))
        FB = np.zeros((nk, nb, nion, 10))
        for k in range(nk):
                it = readnum(txt, it)[1]
                kp, it = readnum(txt, it)
                KP[k, 0] = kp
                kp, it = readnum(txt, it)
                KP[k, 1] = kp
                kp, it = readnum(txt, it)
                KP[k, 2] = kp
                it = readnum(txt, it)[1]
                for band in range(nb):
                        it = readnum(txt, it)[1]
                        BS[k, band], it = readnum(txt, it)
                        it = readnum(txt, it)[1]

                        #read weigths
                        fb = []
                        it = readnum(txt, it)[1]
                        it = readnum(txt, it)[1]
                        for ions in range(nion):
                                it = readnum(txt, it)[1]
                                for i in range(10):
                                        fband, it = readnum(txt, it)
                                        FB[k, band, ions, i] = fband
                        for i in range(10):
                                it = readnum(txt, it)[1]
                        #if SO skip magnetization
                        if SO:
                                for ii in range(3):
                                        for ions in range(nion):
                                                it = readnum(txt, it)[1]
                                                for i in range(10):
                                                        it = readnum(txt, it)[1]
                                        for i in range(10):
                                                it = readnum(txt, it)[1]

                        #skip phases
                        it = readnum(txt, it)[1]
                        it = readnum(txt, it)[1]
                        for ions in range(nion):
                                for i in range(20):
                                        it = readnum(txt, it)[1]
                        for i in range(10):
                                it = readnum(txt, it)[1]
        return KP, BS, FB

def read_DOSCAR(fileName='DOSCAR'):
        with open(fileName,'r') as file:
                lines = [line.rstrip() for line in file]
        Natom, partialDos = np.array([float(line) for line in lines[0].split()])[1:3]
        emin, emax, Np, Efermi, tmp = np.array([float(line) for line in lines[5].split()])
        Natom, Np = int(Natom), int(Np)
        rowLen = np.array([float(line) for line in lines[6].split()]).shape[0]
        DOS = np.zeros((Np,rowLen))
        for i in range(Np):
                DOS[i,:] = np.array([float(line) for line in lines[6+i].split()])
        DOS[:,0]=DOS[:,0]-Efermi
        if partialDos:
                rowLen = np.array([float(line) for line in lines[7+Np].split()]).shape[0]
                PDOS=np.zeros((Np,rowLen,Natom))
                for atom in range(Natom):
                        for i in range(Np):
                                PDOS[i, :, atom] = np.array([float(line) for line in lines[6 + (Np+1)*(atom+1) + i].split()])
                PDOS[:, 0] = PDOS[:, 0] - Efermi
                return DOS,PDOS
        return DOS

def findString(txt,string,last=True,skipNum = 0):
        if last:
                i=txt.count(string)
                txt=txt.replace(string,'',i-1)
        i=txt.find(string)
        i = i + len(string)
        for it in range(skipNum+1):
                a,i=readnum(txt,i)
        return a
