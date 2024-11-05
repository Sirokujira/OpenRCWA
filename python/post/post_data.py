# -*- coding: utf-8 -*-
"""
post_data.py
"""

import math
import matplotlib.pyplot as plt

# ポスト処理制御データを読み込む
def read(fn):

    # ポスト処理制御データ(辞書型)
    Post = { \
        'mloss'      : 0, \
        'iter'       : 0, \
        'feed'       : 0, \
        'point'      : 0, \

        'smith'      : 0, \
        'zin'        : [0, 0, 0, 0], \
        'yin'        : [0, 0, 0, 0], \
        'ref'        : [0, 0, 0, 0], \
        'spara'      : [0, 0, 0, 0], \
        'coupling'   : [0, 0, 0, 0], \
        'f0d'        : [0, 0, 0, 0, 0, 0], \
        'freqdiv'    : 10, \

        'f1ddir'     : [], \
        'f1ddiv'     : [], \
        'f1dangle'   : [], \
        'f1dstyle'   : 0, \
        'f1dcompo'   : [1, 0, 0], \
        'f1ddb'      : 1, \
        'f1dnorm'    : 0, \
        'f1dscale'   : [0, 0, 0, 0], \

        'f2d'        : [0, 0, 0], \
        'f2dcompo'   : [1, 0, 0, 0, 0, 0, 0], \
        'f2ddb'      : 1, \
        'f2dscale'   : [0, 0, 30], \
        'f2dobj'     : 0, \

        'n1dcompo'   : [], \
        'n1ddir'     : [], \
        'n1dpos'     : [], 
        'n1ddb'      : 0, \
        'n1dscale'   : [0, 0, 0, 0], \
        'n1dnoinc'   : 0, \

        'n2ddim'     : [1, 0],  # not used \
        'n2dframe'   : 0,  # not used \
        'n2dcompo'   : [], \
        'n2ddir'     : [], \
        'n2dpos'     : [], \
        'n2ddb'      : 0, \
        'n2dscale'   : [0, 0, 0, 20], \
        'n2dcontour' : 0, \
        'n2dobj'     : 1, \
        'n2dnoinc'   : 0, \
        'n2dzoom'    : [0, 0, 0, 0, 0], \

        'fscale'     : 0, \
        'funit'      : '', \

        'title'      : '', \
        'dt'         : 0, \
        'z0'         : [], \

        'C'          :  2.99792458e8, \
        'MU0'        :  4 * math.pi * 1e-7, \
        'ETA0'       :  0}
    Post['ETA0'] =  Post['C'] * Post['MU0']
    #print(Post)

    # 予め全体を読み込む(UTF-8)
    with open(fn, encoding='utf-8') as fp:
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

    # データ本体
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

        # キーワード
        key = d[0]  #;print(key)

        # データ配列=3番目以降
        p = d[2:]  #;print(p)

        # キーワードによる場合分け
        # 時間波形(2D)
        if   key == 'matchingloss':
            Post['mloss'] = int(p[0])
        elif key == 'plotiter':
            Post['iter'] = int(p[0])
        elif key == 'plotfeed':
            Post['feed'] = int(p[0])
        elif key == 'plotpoint':
            Post['point'] = int(p[0])
        # 周波数特性(2D)
        elif key == 'plotsmith':
            Post['smith'] = int(p[0])
        elif key == 'plotzin':
            Post['zin'][0] = int(p[0])
            if (len(p) > 3):
                Post['zin'][1:4] = [float(p[1]), float(p[2]), int(p[3])]
        elif key == 'plotyin':
            Post['yin'][0] = int(p[0])
            if (len(p) > 3):
                Post['yin'][1:4] = [float(p[1]), float(p[2]), int(p[3])]
        elif key == 'plotref':
            Post['ref'][0] = int(p[0])
            if (len(p) > 3):
                Post['ref'][1:4] = [float(p[1]), float(p[2]), int(p[3])]
        elif key == 'plotspara':
            Post['spara'][0] = int(p[0])
            if (len(p) > 3):
                Post['spara'][1:4] = [float(p[1]), float(p[2]), int(p[3])]
        elif key == 'plotcoupling':
            Post['coupling'][0] = int(p[0])
            if (len(p) > 3):
                Post['coupling'][1:4] = [float(p[1]), float(p[2]), int(p[3])]
        elif key == 'plotfar0d':
            Post['f0d'][0:3] = [float(p[0]), float(p[1]), int(p[2])]
            if (len(p) > 5):
                Post['f0d'][3:6] = [float(p[3]), float(p[4]), int(p[5])]
        elif key == 'freqdiv':
            Post['freqdiv'] = int(p[0])
        # 遠方界面内(2D)
        elif key == 'plotfar1d':
            if len(p) > 1:
                Post['f1ddir'].append(p[0])  # string
                Post['f1ddiv'].append(int(p[1]))
                if len(p) > 2:
                    Post['f1dangle'].append(float(p[2]))
                else:
                    Post['f1dangle'].append(0)
        elif key == 'far1dstyle':
            Post['f1dstyle'] = int(p[0])
        elif key == 'far1dcomponent':
            if len(p) > 2:
                Post['f1dcompo'][0:3] = [int(p[0]), int(p[1]), int(p[2])]
        elif key == 'far1ddb':
            Post['f1ddb'] = int(p[0])
        elif key == 'far1dnorm':
            Post['f1dnorm'] = int(p[0])
        elif key == 'far1dscale':
            if len(p) > 2:
                Post['f1dscale'][0] = 1
                Post['f1dscale'][1] = min(float(p[0]), float(p[1]))
                Post['f1dscale'][2] = max(float(p[0]), float(p[1]))
                Post['f1dscale'][3] = int(p[2])
        # 遠方界全方向(3D)
        elif key == 'plotfar2d':
            Post['f2d'][0] = 1
            if len(p) > 1:
                Post['f2d'][1:3] = [int(p[0]), int(p[1])]
        elif key == 'far2dcomponent':
            if len(p) > 6:
                Post['f2dcompo'] = [int(p[0]), int(p[1]), int(p[2]), int(p[3]), int(p[4]), int(p[5]), int(p[6])]
        elif key == 'far2ddb':
            Post['f2ddb'] = int(p[0])
        elif key == 'far2dscale':
            if len(p) > 1:
                Post['f2dscale'][0] = 1
                Post['f2dscale'][1:3] = [float(p[0]), float(p[1])]  # dB only
        elif key == 'far2dobj':
            Post['f2dobj'] = float(p[0])  # not used
        # 近傍界線上(2D)
        elif key == 'plotnear1d':
            if len(p) > 3:
                Post['n1dcompo'].append(p[0])  # string : 'E'...
                Post['n1ddir'].append(['X', 'Y', 'Z'].index(p[1]))
                Post['n1dpos'].append([float(p[2]), float(p[3])])
        elif key == 'near1ddb':
            Post['n1ddb'] = int(p[0])
        elif key == 'near1dscale':
            if len(p) > 2:
                Post['n1dscale'][0] = 1
                Post['n1dscale'][1] = min(float(p[0]), float(p[1]))
                Post['n1dscale'][2] = max(float(p[0]), float(p[1]))
                Post['n1dscale'][3] = int(p[2])
        elif key == 'near1dnoinc':
            Post['n1dnoinc'] = int(p[0])
        # 近傍界面上(2D)
        elif key == 'plotnear2d':
            if len(p) > 2:
                Post['n2dcompo'].append(p[0])  # string : 'E'...
                Post['n2ddir'].append(['X', 'Y', 'Z'].index(p[1]))
                Post['n2dpos'].append(float(p[2]))
        elif key == 'near2ddim':
            if len(p) > 1:
                Post['n2ddim'] = [int(p[0]), int(p[1])]  # not used (2D only)
        elif key == 'near2dframe':
            Post['n2dframe'] = int(p[0])  # not used
        elif key == 'near2ddb':
            Post['n2ddb'] = int(p[0])
        elif key == 'near2dscale':
            if len(p) > 2:
                Post['n2dscale'][0] = 1
                Post['n2dscale'][1] = min(float(p[0]), float(p[1]))
                Post['n2dscale'][2] = max(float(p[0]), float(p[1]))
                Post['n2dscale'][3] = int(p[2])
        elif key == 'near2dcontour':
            Post['n2dcontour'] = int(p[0])
        elif key == 'near2dobj':
            Post['n2dobj'] = int(p[0])
        elif key == 'near2dnoinc':
            Post['n2dnoinc'] = int(p[0])
        elif key == 'near2dzoom':
            Post['n2dzoom'][0] = 1
            if len(p) > 3:
                Post['n2dzoom'][1] = min(float(p[0]), float(p[1]))
                Post['n2dzoom'][2] = max(float(p[0]), float(p[1]))
                Post['n2dzoom'][3] = min(float(p[2]), float(p[3]))
                Post['n2dzoom'][4] = max(float(p[2]), float(p[3]))

    return Post

# 追加の設定
def setup(Post, Freq1, Freq2):
    #print(Freq1, Freq2)

    # 周波数のバンド
    freq0 = (Freq1[0] + Freq1[-1] + Freq2[0] + Freq2[-1]) / 4
    if   freq0 > 1e12:
        fscale = 1e-12
        funit = '[THz]'
    elif freq0 > 1e9:
        fscale = 1e-9
        funit = '[GHz]'
    elif freq0 > 1e6:
        fscale = 1e-6
        funit = '[MHz]'
    elif freq0 > 1e3:
        fscale = 1e-3
        funit = '[kHz]'
    else:
        fscale = 1
        funit = '[Hz]'
    #print(fscale, funit)

    Post['fscale'] = fscale
    Post['funit'] = funit

    # 遠方界の成分
    Post['farcomp'] = ['E-abs', 'E-theta', 'E-phi', 'E-major', 'E-minor', 'E-RHCP', 'E-LHCP']

    # 遠方界の名前
    if Post['source'] == 0:
        Post['farname'] = 'gain'
    else:
        Post['farname'] = 'cross section'

    # 遠方界の単位
    Post['f0dunit'] = ''
    Post['f1dunit'] = ''
    Post['f2dunit'] = ''
    if Post['source'] == 0:
        # 給電点
        Post['f0dunit'] = '[dBi]'  # dB only
        if Post['f1ddb'] == 1:
            Post['f1dunit'] = '[dBi]' if Post['f1dnorm'] == 0 else '[dB]'
        if Post['f2ddb'] == 1:
            Post['f2dunit'] = '[dBi]'
    else:
        # 平面波入射
        Post['f0dunit'] = '[dBm^2]'  # dB only
        if Post['f1ddb'] == 1:
            Post['f1dunit'] = '[dBm^2]' if Post['f1dnorm'] == 0 else '[dB]'
        else:
            Post['f1dunit'] = '[m^2]' if Post['f1dnorm']== 0 else ''
        if Post['f2ddb'] == 1:
            Post['f2dunit'] = '[dBm^2]'
        else:
            Post['f2dunit'] = '[m^2]'

    return Post
