#!/usr/bin/env python
# coding: utf-8

'''参考文献に記載されているDelphiのコードの移植版
Notes
-----
参考文献: 多数の点を球面上に一様分布させるソフトウェアGSS Generator
https://www.jstage.jst.go.jp/article/geoinformatics/12/1/12_1_3/_article/-char/ja/
Examples
--------
# 600個の点を生成
$ python gss_generator.py 160
'''

from __future__ import division, absolute_import, print_function

import csv
import sys

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def generate_gss(n):
    '''一般化螺旋集合を生成する
    Parameters
    ----------
    n : int
        生成する一般化螺旋集合の点数
    Returns
    -------
    ndarray
        一般化螺旋集合
        生成された天底の点から天頂の点までの球面座標の配列
    '''
    # [[colatitude, longitude], ...]
    gss = np.zeros((n, 2))

    # 天底の点の座標を代入
    gss[0][0] = np.pi
    gss[0][1] = 0.0

    # 天底-天頂間に含まれる点の座標を順に算出
    for index in range(1, n - 1):
        # 参考資料に記載のDelphiのコードでは計算式におけるkを
        # そのまま配列インデックスとして使用するため,
        # 配列サイズを+1として先頭要素は使用せず他関数でも無視されている.
        # 紛らわしいため, このコードでは配列インデックスとkを区別する.
        k = index + 1

        # 参考資料 第(1)式
        # 第(2)式と第(3)式が依存しているパラメータの算出
        hk = -1.0 + 2.0 * (k - 1.0) / (n - 1.0)

        # 参考資料 第(2)式
        # 余緯度の算出
        gss[index][0] = np.arccos(hk)

        # 参考資料 第(3)式
        # 経度の算出
        prev_lon = gss[index - 1][1]
        lon = prev_lon + 3.6 / np.sqrt(n) / np.sqrt(1.0 - hk * hk)
        # 参考資料に記載のDelphiのコードでは経度の値は
        # 増加するままとしているが点数が増えた場合の浮動小数点数の精度に
        # 不安があるため丸め込む
        gss[index][1] = lon % (np.pi * 2)

    # 天頂の点の座標を代入
    gss[n - 1][0] = 0.0
    gss[n - 1][1] = 0.0

    return gss


def spherical2cartesian(coord):
    '''球面座標を直交座標に変換する
    Parameters
    ----------
    spherical_coord : array_like
        天底の点から天頂の点までの球面座標の配列
    Returns
    -------
    cartesian_coord : ndarray
        直交座標(x, y, z)
    '''
    col = coord[0]
    lon = coord[1]
    x = np.sin(col) * np.cos(lon)
    y = np.sin(col) * np.sin(lon)
    z = np.cos(col)
    return np.array([x, y, z])


def cartesian2spherical(coord):
    '''直交行座標を直交座標に変換する
    Parameters
    ----------
    cartesian_coord : array_like
        直交座標(x, y, z)
    Returns
    -------
    spherical_coord : ndarray
        球面座標(余緯度, 経度)
    '''
    x = coord[0]
    y = coord[1]
    z = coord[2]
    colat = np.arccos(z / np.sqrt(x * x + y * y + z * z))
    lon = np.arctan2(y, x)
    return np.array([colat, lon])


def relocate(gss):
    '''天頂と天底の点の座標を一様性が増すように調整する
    Parameters
    ----------
    gss : array_like
        一般化螺旋集合
    Returns
    -------
    relocated_gss : ndarray
        天頂と天底の点の座標が調整された一般化螺旋集合
    '''
    # 資料に第(4)式と第(5)式に当たるが5での除算が
    # 資料に記載のDelphiのコードでは行われておらず単なる総和になっている.
    # ただし, 5での除算を加えて平均の結果を確認したところ
    # 5点の直交座標の平均と総和の結果は異なるが,
    # それらを球面座標に変換した結果は同じなので問題ない.
    n = len(gss)
    relocated_gss = np.array(gss)

    # 天頂と天底の周囲の5点の座標の配列インデックス
    index_of_five_points = np.array([1, 2, 4, 5, 6])

    # 天底の点の座標を調整
    spherical_coords_for_nadir = gss[index_of_five_points]
    cartesian_coords_for_nadir = np.array(
        [spherical2cartesian(c) for c in spherical_coords_for_nadir])
    cartesian_coords_mean_for_nadir = cartesian_coords_for_nadir.sum(axis=0)
    relocated_gss[0] = cartesian2spherical(cartesian_coords_mean_for_nadir)

    # 天頂の点の座標を調整
    spherical_coords_for_zenith = gss[n - (index_of_five_points + 1)]
    cartesian_coords_for_zenith = np.array(
        [spherical2cartesian(c) for c in spherical_coords_for_zenith])
    cartesian_coords_mean_for_zenith = cartesian_coords_for_zenith.sum(axis=0)
    relocated_gss[-1] = cartesian2spherical(cartesian_coords_mean_for_zenith)

    return relocated_gss


def save_scatter2d(gss, filename):
    # 下半球に含まれる点のみを抽出
    n = 0
    for coord in gss:
        if coord[0] <= np.pi * 0.5:
            n += 1
    gss_of_hemisphere = gss[0:n]

    # 直交座標に変換
    cartesian_coords = np.array(
        [spherical2cartesian(c) for c in gss_of_hemisphere])

    # グラフを作成
    plt.figure(figsize=(6, 6))

    # 値の範囲
    plt.xlim((-1.1, 1.1))
    plt.ylim((-1.1, 1.1))

    # 散布図を描画
    plt.scatter(
        x=cartesian_coords.T[1],
        y=cartesian_coords.T[0],
        c='k',
        s=3)

    # 散布図を画像ファイルとして保存
    plt.savefig(filename)
    plt.close()


def save_scatter3d(gss, filename):
    # 直交座標に変換
    cartesian_coords = np.array(
        [spherical2cartesian(c) for c in gss])

    # グラフを作成
    figure = plt.figure(figsize=(6, 6))
    ax = Axes3D(figure)

    # 値の範囲
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(-1.1, 1.1)

    # 3Dの散布図を描画
    ax.scatter3D(
        cartesian_coords.T[0],
        cartesian_coords.T[1],
        cartesian_coords.T[2])

    # 3Dの散布図を画像ファイルとして保存
    plt.savefig(filename)
    plt.close()


def save_cartesian(gss, filename):
    # 直交座標に変換
    cartesian_coords = np.array([spherical2cartesian(c,1) for c in gss])
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        for coord in cartesian_coords:
            writer.writerow(coord)


def save_spherical(gss, filename):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        for coord in gss:
            writer.writerow(coord)


def main():
    if len(sys.argv) < 2:
        sys.stderr.write('n argment is required\n')
        sys.exit(1)

    try:
        n = int(sys.argv[1])
    except ValueError:
        sys.stderr.write('%s is illegal argment\n' % sys.argv[1])
        sys.exit(1)

    if (n < 20):
        sys.stderr.write('%s is lower than 20\n' % sys.argv[1])
        sys.exit(1)

    # 一般化螺旋集合を生成
    gss = generate_gss(n)

    # 天頂と天底の座標を調整した一般化螺旋集合を生成
    relocated_gss = relocate(gss)

    # 2Dの散布図を画像ファイルとして保存
    # 下半球を天底方向から眺めた様子
    save_scatter2d(gss, 'gss.png')
    save_scatter2d(relocated_gss, 'relocated_gss.png')

    # 3Dの散布図を画像ファイルとして保存
    # 球全体を眺めた様子
    save_scatter3d(gss, 'gss_3d.png')
    save_scatter3d(relocated_gss, 'relocated_gss_3d.png')

    # すべての点の球面座標をCSVファイルとして保存
    save_spherical(gss, 'gss_spherical_coords.csv')
    save_spherical(relocated_gss, 'relocated_gss_spherical_coords.csv')

    # すべての点の直交座標をCSVファイルとして保存
    save_cartesian(gss, 'gss_cartesian_coords.csv')
    save_cartesian(relocated_gss, 'relocated_gss_cartesian_coords.csv')

if __name__ == '__main__':
    main()
