/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.util;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月1日 上午5:40:21
 */
public class PARAMETER {

    /**
     * I 方向最大值
     */
    static public final int NXA = 1000;
    /**
     * J 方向最大值
     */
    static public final int NYA = 1000;
    /**
     * 数组最大值
     */
    static public final int NXYA = NXA * NYA;
    /**
     * Wall 最大网格数
     */
    static public final int NWA = 1000;
    /**
     * Symmetry 最大网格数
     */
    static public final int NSA = 1000;
    /**
     * Outlet 最大网格数
     */
    static public final int NOA = 1000;
    /**
     * Inlet 最大网格数
     */
    static public final int NIA = 1000;
    /**
     * O-C 内部网格最大网格数
     */
    static public final int NOCA = 1000;
    /**
     * 最大网格层数
     */
    static public final int MNG = 6;
    /**
     * 最小的数
     */
    static public final double SMALL = 1.E-40;
    /**
     * 最大的数
     */
    static public final double LAGRE = 1.E4;
//
    /**
     * 边界条件SOUTH INDEX
     */
    public static final int SOUTH = 1;

    /**
     * 边界条件NORTH INDEX
     */
    public static final int NORTH = 2;

    /**
     * 边界条件WEST INDEX
     */
    public static final int WEST = 1;

    /**
     * 边界条件EAST INDEX
     */
    public static final int EAST = 2;
//
    /**
     * 重力常数
     */
    public static final double GRAVITY = 9.8;
//
}
