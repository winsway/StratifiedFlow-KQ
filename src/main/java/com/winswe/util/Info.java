/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.util;

/**
 *
 * @author winsw
 */
public class Info {

    public final double Uso;
    public final double Usw;
    public final double expPressureDrop;
    public double codePressureDrop;
    public double hl;

    public Info(double Uso, double Usw, double expPressureDrop) {
        this.Uso = Uso;
        this.Usw = Usw;
        this.expPressureDrop = expPressureDrop;
    }

    @Override
    public String toString() {
        return "Uso\t" + Uso + "\t"
                + "Usw\t" + Usw + "\t"
                + "expPressureDrop\t" + expPressureDrop + "\t"
                + "codePressureDrop\t" + codePressureDrop + "\t"
                + "hl\t" + hl;
    }

    public String fileName() {
        return "Uso = " + Uso + "_" + "Usw = " + Usw;
    }

}
