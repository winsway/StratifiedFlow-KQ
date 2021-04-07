/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.util;

/**
 *
 * @author winsway
 */
public interface Flux {

    /**
     * 计算界面上的phi值 对的加权平均
     *
     * @param refValue 参考值
     * @param dValue 变化值
     * @param u 上游位置
     * @param f 界面位置
     * @param d 下游位置
     * @return 界面插值
     */
    static public double faceValue(double refValue, double dValue, double u, double f, double d) {
        double temp;
        temp = refValue + dValue * (f - u) / (d - u);
        return temp;
    }

}
