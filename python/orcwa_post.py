# -*- coding: utf-8 -*-
"""
OpenRCWA (Python)
ost_post.py : post
"""

import matplotlib.pyplot as plt
import post.load_bin, post.post_data
import post.plot_iter, post.plot_fchar
import post.plot_far1d, post.plot_far2d
import post.plot_near1d, post.plot_near2d

def main():

    # [0-1] 計算結果をloadする
    fn_out = 'orcwa.npz'
    (title, dt, source, fPlanewave, z0,
    Ni, Nj, Nk, N0, NN,
    Xn, Yn, Zn, Freq1, Freq2,
    cEx, cEy, cEz, cHx, cHy, cHz,
    Eiter, Hiter, Iiter, Niter,
    VFeed, IFeed, VPoint, Ntime,
    Zin, Ref, Pin, Spara, Coupling,
    fSurface, cSurface, gline) \
        = post.load_bin.load(fn_out)

    # [0-2] ポスト処理制御データを読み込む
    fn_in = 'python.orcwa'
    Post = post.post_data.read(fn_in)

    # [0-3] 計算データの一部をポスト処理データにコピーする
    Post['title']  = title
    Post['dt']     = dt
    Post['source'] = source
    Post['z0']     = z0
    Post['Ni']     = Ni
    Post['Nj']     = Nj
    Post['Nk']     = Nk
    Post['N0']     = N0
    Post['NN']     = NN

    # [0-4] 追加の設定
    Post = post.post_data.setup(Post, Freq1, Freq2)
    #print(Post)

    # [0-5] ウィンドウの大きさ他
    Post['w2d'] = [6.0, 4.0]          # width, height [inch]
    Post['w3d'] = [5.0, 5.0, 60, 30]  # width, height [inch], theta, phi [deg]
    # フォント
    plt.rcParams['font.family'] = ['sans-serif', 'serif', 'monospace', 'MS gothic', 'MS mincho'][0]  # 0/1/2/3/4
    plt.rcParams['font.size'] = 10  # fontsize [pixel]

    # [1] 時間変化 (2D)
    post.plot_iter.plot(Post, Freq1, Eiter, Hiter, Iiter, Niter, VFeed, IFeed, VPoint, Ntime)

    # [2] 周波数特性　(2D)
    post.plot_fchar.plot(Post, Freq1, Freq2, Zin, Ref, Pin, Spara, Coupling, fSurface, cSurface)

    # [3] 遠方界面内 (2D)
    post.plot_far1d.plot(Post, Freq2, Pin, fSurface, cSurface)

    # [4] 遠方界全方向 (3D)
    post.plot_far2d.plot(Post, Freq2, Pin, fSurface, cSurface, gline)
    
    # [5] 近傍界線上 (2D)
    post.plot_near1d.plot(Post, fPlanewave, Xn, Yn, Zn, Freq2, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN)

    # [6] 近傍界面上 (2D)
    post.plot_near2d.plot(Post, fPlanewave, Xn, Yn, Zn, Freq2, cEx, cEy, cEz, cHx, cHy, cHz, Ni, Nj, Nk, N0, NN, gline)

# enyry point
if __name__ == "__main__":
    main()
