import numpy as np
from itertools import islice


def read_WAVEDER(file='WAVEDER'):
    '''Read the WAVEDER file

        the matrix of the derivative of the cell periodic part
        of the wavefunctions with respect to i k is:
        cder = CDER_BETWEEN_STATES(m,n,1:NK,1:ISP,j)= <u_m|  -i d/dk_j | u_n> = - <u_m| r_j |u_n>
        '''
    from scipy.io import FortranFile
    file = FortranFile(file, 'r')
    nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype=np.int32)
    nodesn_i_dielectric_function = file.read_record(dtype=float)
    wplasmon = file.read_record(dtype=float).reshape(3, 3)
    cder = file.read_record(dtype=np.complex64).reshape(3, ispin, nkpts, nbands_cder, nb_tot)

    return cder, nodesn_i_dielectric_function, wplasmon


def readnum(txt, i, max_len=32):
    ii = max_len
    a = ''
    while txt[i].isdigit() == False and txt[i] != '-':
        i = i + 1
    while (txt[i].isdigit() == True or txt[i] in ('-', '+', '.', 'E')) and ii > 0:
        a = a + txt[i]
        if txt[i] not in ['+', '-']:
            ii -= 1
        i = i + 1
    try:
        a = float(a)
    except:
        a, i = readnum(txt, i, max_len)
    else:
        a = float(a)
    return a, i


def read_EIGENVAL(fileName='EIGENVAL', shift=True):
    with open(fileName, 'r') as file:
        lines = [line.rstrip() for line in file]
    tmp = np.array([float(line) for line in lines[5].split()])
    nk, nb = int(tmp[1]), int(tmp[2])
    rowLen = np.array([float(line) for line in lines[8].split()]).shape[0]
    if rowLen == 3:
        BS = np.zeros((nk, nb))
        OC = np.zeros((nk, nb))
        KP = np.zeros((nk, 4))
        for k in range(nk):
            for b in range(nb):
                BS[k, b] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb + (b + 1)].split()])[1]
                OC[k, b] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb + (b + 1)].split()])[2]
            KP[k, :] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb].split()])[0:4]
        if shift:
            BS -= BS[0, int(round(OC[0, :].sum())) - 1]
    elif rowLen == 5:
        BS = np.zeros((nk, nb, 2))
        OC = np.zeros((nk, nb, 2))
        KP = np.zeros((nk, 4))
        for k in range(nk):
            for b in range(nb):
                BS[k, b, :] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb + (b + 1)].split()])[1:3]
                OC[k, b, :] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb + (b + 1)].split()])[3:5]
            KP[k, :] = np.array([float(line) for line in lines[5 + 2 * (k + 1) + k * nb].split()])[0:4]
        if shift:
            BS[:, :, 0] -= BS[0, int(round(OC[0, :, 0].sum())) - 1, 0]
            BS[:, :, 1] -= BS[0, int(round(OC[0, :, 1].sum())) - 1, 1]
    else:
        print('Strange rowLen')
    return KP, BS, OC


def calcP(k, nb1_start, nb1_end, nb2_start, nb2_end, km, BS):
    P0x = 0.
    P0y = 0.
    P0z = 0.
    for nb in range(nb1_start, nb1_end + 1):
        for nb2 in range(nb2_start, nb2_end + 1):
            '''print(abs(km[0,0,0,nb,nb2]))
                        print(abs(km[1,0,0,nb,nb2]))
                        print(abs(km[2,0,0,nb,nb2]))
                        print(abs(np.array([km[0,0,0,nb,nb2],km[1,0,0,nb,nb2],km[2,0,0,nb,nb2]])))
                        print((abs(km[0,0,0,nb,nb2])+abs(km[1,0,0,nb,nb2])+abs(km[2,0,0,nb,nb2]))/3.)'''
            P0x = P0x + (abs(km[0, 0, k, nb, nb2]) * (BS[k, nb] - BS[k, nb2])) ** 2
            P0y = P0y + (abs(km[1, 0, k, nb, nb2]) * (BS[k, nb] - BS[k, nb2])) ** 2
            P0z = P0z + (abs(km[2, 0, k, nb, nb2]) * (BS[k, nb] - BS[k, nb2])) ** 2
    return (P0x ** 0.5, P0y ** 0.5, P0z ** 0.5)


def read_POSCAR(fileName='POSCAR'):
    with open(fileName, 'r') as file:
        lines = [line.rstrip() for line in file]
    acell = float(lines[1])
    a1 = np.array([float(line) for line in lines[2].split()])
    a2 = np.array([float(line) for line in lines[3].split()])
    a3 = np.array([float(line) for line in lines[4].split()])
    atom = np.array([f for f in lines[5].split()])
    Natom = np.array([int(line) for line in lines[6].split()])
    natom = Natom.sum()
    coords = np.zeros((natom, 3))
    for i in range(natom):
        coords[i] = np.array([float(line) for line in lines[8 + i].split()])
    return np.concatenate(((acell * a1)[:, np.newaxis], (acell * a2)[:, np.newaxis], (acell * a3)[:, np.newaxis]),
                          axis=1).T, atom, Natom, coords


def write_POSCAR(fileName='POSCAR', a=[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], atom=['H'], Natom=[1],
                 coords=[[0., 0., 0.]], coordsSys='D'):
    with open(fileName, 'w') as file_write:
        txt = '''auto
1.0
%f %f %f
%f %f %f
%f %f %f
''' % (a[0, 0], a[0, 1], a[0, 2], a[1, 0], a[1, 1], a[1, 2], a[2, 0], a[2, 1], a[2, 2])
        for a in atom:
            txt = txt + a + ' '
        txt = txt + '\n'
        for n in Natom:
            txt = txt + str(n) + ' '
        txt = txt + '\n'
        txt = txt + coordsSys + '\n'
        for c in coords:
            txt = txt + '%f %f %f\n' % (c[0], c[1], c[2])
        file_write.write(txt)


def read_PROCAR(file='PROCAR', SO=False, G_only=False, Phases=False):
    f = open(file)
    txt = f.read()
    txt = txt + ' !!!'
    f.close()
    it = txt.index('k-points')
    # for header in range(2):
    #        it = readnum(txt, it)[1]
    nk, it = readnum(txt, it)
    if G_only == True:
        nk = 1
    nb, it = readnum(txt, it)
    nion, it = readnum(txt, it)
    nk, nb, nion = int(nk), int(nb), int(nion)
    BS = np.zeros((nk, nb))
    KP = np.zeros((nk, 3))
    FB = np.zeros((nk, nb, nion, 10))
    for k in range(nk):
        it = readnum(txt, it)[1]
        kp, it = readnum(txt, it, max_len=10)
        KP[k, 0] = kp
        kp, it = readnum(txt, it, max_len=10)
        KP[k, 1] = kp
        kp, it = readnum(txt, it, max_len=10)
        KP[k, 2] = kp
        it = readnum(txt, it)[1]
        for band in range(nb):
            it = readnum(txt, it)[1]
            BS[k, band], it = readnum(txt, it)
            it = readnum(txt, it)[1]

            # read weigths
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
            # if SO skip magnetization
            if SO:
                for ii in range(3):
                    for ions in range(nion):
                        it = readnum(txt, it)[1]
                        for i in range(10):
                            it = readnum(txt, it)[1]
                    for i in range(10):
                        it = readnum(txt, it)[1]

            # skip phases
            if Phases:
                it = readnum(txt, it)[1]
                it = readnum(txt, it)[1]
                for ions in range(nion):
                    for i in range(20):
                        it = readnum(txt, it)[1]
                for i in range(10):
                    it = readnum(txt, it)[1]
    return KP, BS, FB


def read_DOSCAR(fileName='DOSCAR', shift=True):
    with open(fileName, 'r') as file:
        lines = [line.rstrip() for line in file]
    Natom, partialDos = np.array([float(line) for line in lines[0].split()])[1:3]
    emin, emax, Np, Efermi, tmp = np.array([float(line) for line in lines[5].split()])
    Natom, Np = int(Natom), int(Np)
    rowLen = np.array([float(line) for line in lines[6].split()]).shape[0]
    DOS = np.zeros((Np, rowLen))
    for i in range(Np):
        DOS[i, :] = np.array([float(line) for line in lines[6 + i].split()])
    if partialDos:
        rowLen = np.array([float(line) for line in lines[7 + Np].split()]).shape[0]
        PDOS = np.zeros((Np, rowLen, Natom))
        for atom in range(Natom):
            for i in range(Np):
                PDOS[i, :, atom] = np.array([float(line) for line in lines[6 + (Np + 1) * (atom + 1) + i].split()])
        if shift:
            DOS[:, 0] = DOS[:, 0] - Efermi
            PDOS[:, 0] = PDOS[:, 0] - Efermi
            Efermi = 0.
        return DOS, PDOS, Efermi
    if shift:
        DOS[:, 0] = DOS[:, 0] - Efermi
        Efermi = 0.
    return DOS, Efermi


def findString(txt, string, last=True, skipNum=0):
    if last:
        i = txt.count(string)
        txt = txt.replace(string, '', i - 1)
    i = txt.find(string)
    i = i + len(string)
    for it in range(skipNum + 1):
        a, i = readnum(txt, i)
    return a


def read_LOCPOT(fileName='LOCPOT'):
    with open(fileName, 'r') as file:
        lines = [line.rstrip() for line in file]
    Natom = np.array([float(line) for line in lines[6].split()]).sum()
    Natom = int(Natom)
    Nx, Ny, Nz = np.array([float(line) for line in lines[Natom + 9].split()])
    Nx, Ny, Nz = int(Nx), int(Ny), int(Nz)
    POT = np.zeros(Nx * Ny * Nz)
    rowlen = len(np.array([float(line) for line in lines[Natom + 10].split()]))
    for i in range(Nx * Ny * Nz // rowlen):
        POT[i * rowlen:i * rowlen + rowlen] = np.array([float(line) for line in lines[i + Natom + 10].split()])
    j = Nx * Ny * Nz - Nx * Ny * Nz % rowlen
    POT[j:] = np.array([float(line) for line in lines[i + 1].split()])[:Nx * Ny * Nz % rowlen]
    POT = POT.reshape((Nz, Ny, Nx))
    return POT


def read_lines_between(filename, start_line, end_line):
    """
        Read lines from start_line to end_line (0-based, end_line exclusive).
        """
    with open(filename, 'r') as f:
        lines = list(islice(f, start_line, end_line))
    return lines


'''def read_BSEFATBAND(fileName='BSEFATBAND'):
        with open(fileName,'r') as file:
                lines = [line.rstrip() for line in file]
        Msize, EigNum = np.array([int(line) for line in lines[0].split()])
        BSEeigVal, IPeigVal = np.zeros(EigNum), np.zeros(EigNum)
        KP,BS,C = np.zeros((EigNum,Msize,3)), np.zeros((EigNum,Msize,2)), np.zeros((EigNum,Msize,1))
        for eig in range(EigNum):
                BSEeigVal[eig], IPeigVal[eig] = np.array([line for line in lines[1+eig*(Msize+1)].split()][2::2])
                for state in range(Msize):
                        KP[eig,state] = np.array([line for line in lines[1+eig*(Msize+1)+state+1].split()][:3],dtype=float)
                        BS[eig,state] = np.array([line for line in lines[1+eig*(Msize+1)+state+1].split()][3:5],dtype=float)
                        C[eig,state] = np.array([line for line in lines[1+eig*(Msize+1)+state+1].split()][5],dtype=float)
        return KP,BS,C,BSEeigVal,IPeigVal'''


def read_BSEFATBAND(fileName='BSEFATBAND', istate=0):
    with open(fileName, 'r') as f:
        lines = list(islice(f, 0, 1))
    Msize, EigNum = np.array([int(line) for line in lines[0].split()])
    with open(fileName, 'r') as f:
        lines = list(islice(f, 1 + (Msize + 1) * istate, 1 + (Msize + 1) * (istate + 1)))

    KP, BS, C = np.zeros((Msize, 3)), np.zeros((Msize, 2)), np.zeros((Msize))

    BSEeigVal, IPeigVal = np.array([line for line in lines[0].split()][2::2], dtype=float)

    for state in range(Msize):
        KP[state] = np.array([line for line in lines[1 + state].split()][:3], dtype=float)
        BS[state] = np.array([line for line in lines[1 + state].split()][3:5], dtype=float)
        C[state] = np.array([line for line in lines[1 + state].split()][5], dtype=float)
    return KP, BS, C, BSEeigVal, IPeigVal