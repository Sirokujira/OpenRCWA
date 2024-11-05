"""
input_data.py
"""

import math
import numpy as np

# 定数
def const():
    Parm = {}
    Parm['C']    = 2.99792458e8                      # 真空の光速
    Parm['MU0']  = 4 * math.pi * 1e-7                # 真空の透磁率
    Parm['EPS0'] = 1 / (Parm['C']**2 * Parm['MU0'])  # 真空の誘電率
    Parm['ETA0'] = Parm['C'] * Parm['MU0']           # 真空の波動インピーダンス
    Parm['EPS']  = 1e-6                              # 無次元微小量
    Parm['PEC']  = 1                                 # PEC(完全導体)

    return Parm

# データ入力
def read(fn, Parm):

    # error code
    #ierr = 0
    errmsg = "*** invalid %s data : line %d"

    # パラメーター追加
    Parm['title']      = ''
    Parm['source']     = 0
    Parm['planewave']  = [0] * 3
    Parm['rfeed']      = 0
    Parm['abc']        = [0] * 4
    Parm['pbc']        = [0] * 3
    Parm['prop']       = 0
    Parm['solver']     = [1000, 50, 1e-3]
    Parm['dt']         = 0
    Parm['tw']         = 0
    Parm['plot3dgeom'] = 0 

    # 配列宣言
    Xn = None
    Yn = None
    Zn = None
    Freq1 = None
    Freq2 = None

    # 予め全体を読み込む
    with open(fn, 'rt', encoding='utf-8') as fp:
        header = fp.readline()
        tlines = fp.readlines()

    # ヘッダー処理
    d = header.strip().split()
    if (len(d) < 3) or (d[0] != 'OpenRCWA'):
        print('*** Not OpenRCWA data : %s' % fn)
        return
    version = (100 * int(d[1])) + int(d[2])
    if version < 400:
        print('*** version of data is old : %s.%s < 4.0' % (d[1], d[2]))
        return

    # 配列の大きさを調べる
    nmaterial = 2
    ngeometry = 0
    nfeed = 0
    ninductor = 0
    npoint = 0
    for tline in tlines:
        if   tline.lstrip().startswith('material = '):
            nmaterial += 1
        elif tline.lstrip().startswith('geometry = '):
            ngeometry += 1
        elif tline.lstrip().startswith('feed = '):
            nfeed += 1
        elif tline.lstrip().startswith('load = '):
            d = tline.strip().split()
            if len(d) > 7:
                if (d[6] == 'R') or (d[6] == 'C'):
                    # RとCは等価な物性値を持つ物体と考える
                    nmaterial += 1
                    ngeometry += 1
                elif d[6] == 'L':
                    ninductor += 1
        elif tline.lstrip().startswith('point = '):
            npoint += 1

    # 観測点1+/1-を追加
    if npoint > 0:
        npoint += 2
    #print(nmaterial, ngeometry, nfeed, ninductor, npoint)

    # 配列作成
    iMaterial  = np.zeros( nmaterial,      'i4')  # m
    fMaterial  = np.zeros((nmaterial, 11), 'f8')  # epsr, esgm, amur, msg, ... 
    iGeometry  = np.zeros((ngeometry,  2), 'i4')  # m, shape
    fGeometry  = np.zeros((ngeometry,  8), 'f8')  # x1, x2, y1, y2, z1, z2, ...
    iFeed      = np.zeros((nfeed,      4), 'i4')  # idir, i, j, k
    fFeed      = np.zeros((nfeed,      9), 'f8')  # x, y, z, v, delay, z0, dx, dy, dz
    iInductor  = np.zeros((ninductor,  4), 'i4')  # idir, i, j, k
    fInductor  = np.zeros((ninductor, 10), 'f8')  # x, y, z, dx, dy, dz, L, ...
    iPoint     = np.zeros((npoint,     4), 'i4')  # idir, i, j, k
    fPoint     = np.zeros((npoint,     6), 'f8')  # x, y, z, dx, dy, dz

    # 物性値、0=真空、1=PEC
    iMaterial[0:2] = 1
    fMaterial[0:2, 0] = [1] * 2  # epsr
    fMaterial[0:2, 2] = [1] * 2  # amur
    #print(iMaterial, fMaterial)

    # データ本体
    nmaterial = 2
    ngeometry = 0
    nfeed = 0
    ninductor = 0
    npoint = 0
    iline = 1
    for tline in tlines:
        # 行番号
        iline += 1

        # token分解
        d = tline.strip().split()  #;print(d)
        
        # check
        if d[0].startswith('#'):  # コメント行は飛ばす
            continue
        elif (d[0] == 'end'):  # end行で終了
            break
        elif (len(d) < 3):  # データ数3個未満は飛ばす
            continue
        elif (d[1] != "="):  # 第2データが=以外は飛ばす
            continue

        # keyword
        key = d[0]  #;print(key)

        # データ配列=3番目以降
        p = d[2:]  #;print(p)

        # キーワードによる場合分け
        if key == "title":
            Parm['title'] = tline[len("title = "):].rstrip()  #;print(Parm['title'])
        elif key == 'xmesh':
            Nx, Xn = _makemesh(p)
        elif key == 'ymesh':
            Ny, Yn = _makemesh(p)
        elif key == 'zmesh':
            Nz, Zn = _makemesh(p)
        elif key == 'material':
            if (len(p) > 4):
                #print(p)
                iMaterial[nmaterial] = int(p[0])
                if   p[0] == '1':  # 通常媒質
                    fMaterial[nmaterial, 0] = float(p[1])
                    fMaterial[nmaterial, 1] = float(p[2])
                    fMaterial[nmaterial, 2] = float(p[3])
                    fMaterial[nmaterial, 3] = float(p[4])
                elif p[0] == '2':  # 分散性媒質
                    fMaterial[nmaterial, 0] = 1
                    fMaterial[nmaterial, 1] = 0
                    fMaterial[nmaterial, 2] = 1
                    fMaterial[nmaterial, 3] = 0
                    fMaterial[nmaterial, 4:8] = [float(pos) for pos in p[1:5]]
                nmaterial += 1
            else:
                print(errmsg % (key, iline))
        elif key == 'geometry':
            if (len(p) > 7):
                if len(p) == 8:
                    p += [0] * 2  # 座標データが6個のときは2個ダミーを加える
                iGeometry[ngeometry, :] = [int(p[0]), int(p[1])]
                fGeometry[ngeometry, :] = [float(pos) for pos in p[2:10]]
                ngeometry += 1
            else:
                print(errmsg % (key, iline))
        elif key== "feed":
            if len(p) > 6:
                Parm['source'] = 0
                iFeed[nfeed, 0] = ['X', 'Y', 'Z'].index(p[0])
                fFeed[nfeed, 0:3] = [float(p[1]), float(p[2]), float(p[3])]
                fFeed[nfeed, 3:5] = [float(p[4]), float(p[5])]
                fFeed[nfeed, 5] = float(p[6])
                nfeed += 1
            else:
                print(errmsg % (key, iline))
        elif key == "planewave":
            if len(p) > 2:
                Parm['source'] = 1
                Parm['planewave'][0:3] = [float(p[0]), float(p[1]), int(p[2])]
            else:
                print(errmsg % (key, iline))
        elif key == "point":
            if len(p) > 3:
                iPoint[npoint, 0] = ['X', 'Y', 'Z'].index(p[0])
                fPoint[npoint, 0:3] = [float(p[1]), float(p[2]), float(p[3])]
                if (npoint == 0) and (len(p) > 4):
                    # ポート1の伝搬方向=0...5)
                    Parm['prop'] = ['+X', '-X', '+Y', '-Y', '+Z', '-Z'].index(p[4])
                npoint += 1
            else:
                print(errmsg % (key, iline))
        elif key == "load":
            #print(p)
            if len(p) > 5:
                idir = ['X', 'Y', 'Z'].index(p[0])  # 0/1/2
                pos = [float(p[1]), float(p[2]), float(p[3])]  # (X,Y,Z)
                RCL = p[4]  # 'R'/'C'/'L'
                zload = float(p[5])
                if ((RCL == 'R') or (RCL == 'C')) and (zload > 0):
                    # R,Cのときは等価な電気定数を持つ物体に置き換える(形状番号=90...95)
                    # その他のデータは後で代入する
                    iGeometry[ngeometry, 0] = nmaterial
                    iGeometry[ngeometry, 1] = (90 if (RCL =='R') else 93) + idir
                    fGeometry[ngeometry, 0:3] = pos
                    fGeometry[ngeometry, 3] = zload  # R/C
                    ngeometry += 1
                    nmaterial += 1
                elif (RCL == 'L') and (zload > 0):
                    # inductor
                    iInductor[ninductor, 0] = idir
                    fInductor[ninductor, 0:3] = pos
                    fInductor[ninductor, 6] = zload  # L
                    ninductor += 1
            else:
                print(errmsg % (key, iline))
        elif key == "rfeed":
            Parm['rfeed'] = float(p[0])
        elif key == "abc":
            if   (p[0] == "0"):
                Parm['abc'] = [0] * 4
            elif (p[0] == "1") and (len(p) > 3):
                Parm['abc'] = [1, int(p[1]), float(p[2]), float(p[3])]
            else:
                print(errmsg % (key, iline))
        elif key == "pbc":
            if len(p) > 2:
                Parm['pbc'] = [int(p[0]), int(p[1]), int(p[2])]
            else:
                print(errmsg % (key, iline))
        elif key == "frequency1":
            if len(p) > 2:
                [f0, f1, fdiv] = [float(p[0]), float(p[1]), int(p[2])]
                if fdiv >= 0:
                    df = (f1 - f0) / fdiv if (fdiv > 0) else 0
                    Freq1 = np.empty(fdiv + 1, np.float64)
                    for i in range(fdiv + 1):
                        Freq1[i] = f0 + (i * df)
            else:
                print(errmsg % (key, iline))
        elif key == "frequency2":
            if len(p) > 2:
                [f0, f1, fdiv] = [float(p[0]), float(p[1]), int(p[2])]
                if fdiv >= 0:
                    df = (f1 - f0) / fdiv if (fdiv > 0) else 0
                    Freq2 = np.empty(fdiv + 1, np.float64)
                    for i in range(fdiv + 1):
                        Freq2[i] = f0 + (i * df)
            else:
                print(errmsg % (key, iline))
        elif key == "solver":
            if len(p) > 2:
                Parm['solver'] = [int(p[0]), int(p[1]), float(p[2])]
                # 最大反復回数は出力間隔の整数倍
                if (Parm['solver'][0] % Parm['solver'][1]) != 0:
                    Parm['solver'][0] = (Parm['solver'][0] // Parm['solver'][1]) * Parm['solver'][1]
            else:
                print(errmsg % (key, iline))
        elif key == "timestep":
            Parm['dt'] = float(p[0])
        elif key == "pulsewidth":
            Parm['tw'] = float(p[0])
        elif key == "plot3dgeom":
            Parm['plot3dgeom'] = int(p[0])

    # error check

    # Xメッシュ
    if (Xn is None) or (len(Xn) < 2):
        print("*** no xmesh data")
        #ierr = 1
    else:
        for i in range(len(Xn) - 1):
            if Xn[i] >= Xn[i + 1]:
                print("*** invalid xmesh data")
                #ierr = 1
                break

    # Yメッシュ
    if (Yn is None) or (len(Yn) < 2):
        print("*** no ymesh data")
        #ierr = 1
    else :
        for j in range(len(Yn) - 1):
            if Yn[j] >= Yn[j + 1]:
                print("*** invalid ymesh data")
                #ierr = 1
                break

    # Zメッシュ
    if (Zn is None) or (len(Zn) < 2):
        print("*** no zmesh data")
        #ierr = 1
    else:
        for k in range(len(Zn) - 1):
            if Zn[k] >= Zn[k + 1]:
                print("*** invalid zmesh data")
                #ierr = 1
                break

    # 物体形状の物性値番号のチェック
    for n in range(iGeometry.shape[0]):
        m = iGeometry[n, 0]
        if m >= len(iMaterial):
            print("*** too large material id = %d: geometry #%d" % (m , n + 1))
            #ierr = 1

    # 給電点と平面波入射はどちらか一方のみ
    if (nfeed == 0) and (Parm['source'] == 0):
        print("*** no source")
        #ierr = 1
    if (nfeed >  0) and (Parm['source'] == 1):
        print("*** feed and planewave")
        #ierr = 1

    # 収束判定条件
    if (Parm['solver'][0] <= 0) or (Parm['solver'][1] <= 0):
        print("*** invalid solver data")
        #ierr = 1

    # 周波数
    if Freq1 is None:
        print("*** no frequency1 data")
        #ierr = 1
    if Freq2 is None:
        print("*** no frequency2 data")
        #ierr = 1

    # warnings

    # PBCのときはMur-1stのみサポート
    if ((Parm['abc'][0] == 1) and \
        (Parm['pbc'][0] > 0 or Parm['pbc'][1] > 0 or Parm['pbc'][2] > 0)):
        print("*** warning : PBC -> Mur-1st")
        Parm['abc'][0] = 0

    #print(Parm)
    #print(iMaterial)
    #print(fMaterial)
    #print(iGeometry)
    #print(fGeometry)
    #print(iFeed)
    #print(fFeed)
    #print(iPoint)
    #print(fPoint)
    #print(iInductor)
    #print(fInductor)

    return \
    Nx, Ny, Nz, Xn, Yn, Zn, \
    iMaterial, fMaterial, iGeometry, fGeometry, iFeed, fFeed, \
    iPoint, fPoint, iInductor, fInductor, Freq1, Freq2

# (private) 節点データ作成
def _makemesh(p):
    EPS = 1e-12

    # error check
    if len(p) < 3:  # データ数は最低3個
        print('*** number of mesh data < 3')
        return None
    if (len(p) % 2 == 1) and (float(p[0]) > float(p[-1])):
        p = p[::-1]  # 座標が降順なら昇順に変換する

    # 区間配列作成
    pos = np.zeros(0, float)
    div = np.zeros(0, int)
    for i, d in enumerate(p):
        if i % 2 == 0:
            # 区間座標
            if (i > 0) and (abs(float(d) - pos[-1]) < EPS):
                # 区間幅が0ならその区間を削除する
                div = np.delete(div, i // 2 - 1)
                continue
            pos = np.append(pos, float(d))
        else:
            # 区間分割数
            if int(d) <= 0:
                # 区間分割数は正
                print('*** division of mesh <= 0')
                return None
            div = np.append(div, int(d))

    # error check
    if (len(pos) < 2) or (len(div) < 1):
        # 区間座標は最低2個、区間分割数は最低1個必要
        print('*** no valid mesh data')
        return None

    # 節点データ作成
    node = np.array(pos[0], float)  # 最初
    for i in range(len(div)):
        dl = (pos[i + 1] - pos[i]) / div[i]
        node = np.append(node, np.linspace(pos[i] + dl, pos[i + 1], div[i]))
    #print(node)
    
    return len(node) - 1, node
