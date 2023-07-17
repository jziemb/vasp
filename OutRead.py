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

def read_EIGENVAL(file='EIGENVAL'):
    f = open(file)
    txt = f.read()
    it = 0
    for header in range(11):
        it = readnum(txt, it)[1]
    nk, it = readnum(txt, it)
    nb, it = readnum(txt, it)
    nk, nb = int(nk), int(nb)
    BS = []
    KP = []
    for k in range(nk):
        kpoint = []
        kp, it = readnum(txt, it)
        kpoint.append(kp)
        kp, it = readnum(txt, it)
        kpoint.append(kp)
        kp, it = readnum(txt, it)
        kpoint.append(kp)
        it = readnum(txt, it)[1]
        bands = []
        for band in range(nb):
            it = readnum(txt, it)[1]
            band, it = readnum(txt, it)
            it = readnum(txt, it)[1]
            bands.append(band)
        BS.append(bands)
        KP.append(kpoint)
    f.close()
    return np.array(KP), np.array(BS)

def calcP(nb1_start,nb1_end,nb2_start,nb2_end):
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

