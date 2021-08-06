/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.onedimension;

import static java.lang.Math.PI;
import static java.lang.Math.pow;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月24日 上午10:03:40
 */
public class TwophaseFlowEquation {

    /**
     * volumetric flowrate m³/s
     */
    private final double volumeFlowrate;

    private final double densityW, densityO;
    private final double dynamicViscosityW, dynamicViscosityO;
    /**
     * diameter m
     */
    private final double diameter;
    /**
     * gravity
     */
    public final static double gravity = 9.8;
    /**
     * pipe cross area，㎡
     */
    public double pipeCrossArea;
    /**
     * average velocity，m/s
     */
    private final double averageVelocity;

    private final double hl;

    private double Si;

    private double So;

    private double Sw;

    private double Ao;

    private double Aw;

    private double Ho;

    private double Hw;

    private double Do;

    private double Dw;

    private double H;

    private double Uso;

    private double Usw;

    private double Uo;
    private double Uw;

    private double Qo;
    private double Qw;

    public TwophaseFlowEquation(
            double Usw,
            double Uso,
            double diameter,
            double densityW,
            double densityO,
            double dynamicViscosityW,
            double dynamicViscosityO,
            double hl) {
        this.Usw = Usw;
        this.Uso = Uso;

        this.densityW = densityW;
        this.densityO = densityO;
        this.dynamicViscosityW = dynamicViscosityW;
        this.dynamicViscosityO = dynamicViscosityO;
        this.diameter = diameter;
        this.hl = hl;
        pipeCrossArea = PI / 4.0 * diameter * diameter;
        this.Qo = this.Uso * pipeCrossArea;
        this.Qw = this.Usw * pipeCrossArea;
        this.volumeFlowrate = (Usw + Uso) * pipeCrossArea;
        averageVelocity = volumeFlowrate / pipeCrossArea;
        this.H = 2 * hl - 1.0;
        this.Si = diameter * Math.sqrt(1 - H * H);
        this.So = diameter * Math.acos(H);
        this.Sw = diameter * PI - So;
        this.Ao = 0.25 * diameter * (So - Si * H);
        this.Aw = pipeCrossArea - Ao;
        this.Ho = this.Ao / this.pipeCrossArea;
        this.Hw = this.Aw / this.pipeCrossArea;
        this.Uo = Uso / Ho;
        this.Uw = Usw / Hw;

        if (Math.abs(Uw - Uo) < 0.05) {
            this.Do = 4 * Ao / So;
            this.Dw = 4 * Aw / Sw;
        } else if (Uw > Uo) {
            this.Do = 4 * Ao / (So + Si);
            this.Dw = 4 * Aw / (Sw + Si);
        } else if (Uw < Uo) {
            this.Do = 4 * Ao / So;
            this.Dw = 4 * Aw / Sw;
        }

    }

    public double getOilReynoldsNumber() {
        return TwophaseFlowEquation.calculateReynoldsNumberByAverageVelocity(
                this.Uso,
                this.Do,
                this.densityO,
                this.dynamicViscosityO
        );
    }

    public double getWaterReynoldsNumber() {
        return TwophaseFlowEquation.calculateReynoldsNumberByAverageVelocity(
                this.Usw,
                this.Dw,
                this.densityW,
                this.dynamicViscosityW
        );
    }

    /**
     * U=4Q/(&pi;D<sup>2</sup>)
     *
     * @param volumeFlowrate volume flowrate m³/s
     * @param diameter pipe diamter m
     * @return average velocity m/s
     */
    public static double calculateAverageVelocityByVolumeFlowrate(
            double volumeFlowrate,
            double diameter
    ) {
        return 4 * volumeFlowrate / (PI * diameter * diameter);
    }

    /**
     * Q =<code>Re&PI;D&mu;/(4.0&rho;)</code>
     *
     * @param Re Reynolds number
     * @param diameter pipe diameter m
     * @param density density kg/m³
     * @param dynamicVisosity dynamic viscosity Pa.s
     * @return volume flowrate m³/s
     */
    public static double calculateVolumeFlowrateByReynoldsNumber(
            double Re,
            double diameter,
            double density,
            double dynamicVisosity
    ) {
        return Re * PI * diameter * dynamicVisosity / (4.0 * density);
    }

    /**
     * @param averageVelocity average velocity m/s
     * @param diameter pipe diameter m
     * @return volume flowrate m³/s
     */
    public static double calculateVolumeFlowrateByAverageVelocity(
            double averageVelocity,
            double diameter
    ) {
        return averageVelocity * PI * diameter * diameter / 4.0;
    }

    /**
     *
     * @param density density kg/m³
     * @param averageVelocity average velocity m/s
     * @param diameter pipe diameter m
     * @param mu dynamic viscosity Pa s
     * @return ReynoldsNumber
     */
    public static double calculateReynoldsNumberByAverageVelocity(
            double averageVelocity,
            double diameter,
            double density,
            double mu
    ) {
        return density * averageVelocity * diameter / mu;
    }

    /**
     *
     *
     * @param volumeFlowrate volmue flowrate m³/s
     * @param density density kg/m³
     * @param diameter pipe diameter m
     * @param mu dynamic viscosity Pa s
     * @return ReynoldsNumber
     */
    public static double calculateReynoldsNumberByDynamicViscosity(
            double volumeFlowrate,
            double diameter,
            double density,
            double mu
    ) {
        return (4 * volumeFlowrate * density) / (PI * diameter * mu);
    }

    /**
     *
     *
     * @param volumeFlowrate volume flowrate m³/s
     * @param nu kinematic viscosity ㎡/s
     * @param diameter pipe diameter m
     * @return ReynoldsNumber
     */
    public static double calculateReynoldsNumberByKinematicViscosity(
            double volumeFlowrate,
            double diameter,
            double nu
    ) {
        return (4 * volumeFlowrate) / (PI * diameter * nu);
    }

    /**
     *
     * @return
     */
    public double getVolumeFlowrate() {
        return volumeFlowrate;
    }

    /**
     *
     * @return
     */
    public double getDiameter() {
        return diameter;
    }

    /**
     *
     * @return
     */
    public static double getGravity() {
        return gravity;
    }

    /**
     *
     * @return
     */
    public double getPipeCrossArea() {
        return pipeCrossArea;
    }

    /**
     *
     * @return
     */
    public double getAverageVelocity() {
        return averageVelocity;
    }

    public double getQo() {
        return Qo;
    }

    public double getQw() {
        return Qw;
    }

}
