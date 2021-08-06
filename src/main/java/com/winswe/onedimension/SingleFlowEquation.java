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
public class SingleFlowEquation {

    /**
     * volumetric flowrate m³/s
     */
    private final double volumeFlowrate;
    /**
     * Density kg/m³
     */
    private final double density;
    /**
     * &nu;-dynamic viscosity Pa.s
     */
    private final double dynamicViscosity;
    /**
     * diameter m
     */
    private final double diameter;
    /**
     * wall rough intensity m
     */
    private final double roughness;
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

    /**
     * chose different method to calculate pressure drop
     */
    public enum Type {

        /**
         * darcy formula
         */
        Darcy,
        /**
         * repizon formula
         */
        Repinzon
    }

    /**
     *
     * constructure function
     *
     * @param volumeFlowrate volume flowrate m³/s
     * @param density density kg/m³
     * @param dynamicViscosity dynamic viscosity Pa.s
     * @param diameter diameter m
     * @param roughness roughness m
     */
    public SingleFlowEquation(
            double volumeFlowrate,
            double diameter,
            double density,
            double dynamicViscosity,
            double roughness) {
        this.volumeFlowrate = volumeFlowrate;
        this.density = density;
        this.dynamicViscosity = dynamicViscosity;
        this.diameter = diameter;
        this.roughness = roughness;
        pipeCrossArea = PI / 4.0 * diameter * diameter;
        averageVelocity = volumeFlowrate / pipeCrossArea;
    }

    /**
     *
     * @return Reynolds Number
     */
    public double getReynoldsNumber() {
        return SingleFlowEquation.calculateReynoldsNumberByKinematicViscosity(
                volumeFlowrate,
                diameter,
                dynamicViscosity / density
        );
    }

    /**
     *
     * @param length length m
     * @param type Darcy or Repinzon
     * @return h m
     */
    public double getHt(double length, Type type) {
        double Re;
        double lambda;
        double ht = 0;
        Re = getReynoldsNumber();
        if (type == Type.Darcy) {
            lambda = this.calculateLambda(Re, diameter, roughness);
            ht = this.ht(lambda, length, diameter, averageVelocity);
        } else if (type == Type.Repinzon) {
            ht
                    = this.RepinzonMethod(Re, volumeFlowrate, diameter, roughness, dynamicViscosity / density)
                    * length;
        }
        return ht;
    }

    /**
     *
     * h=&lambda;<sup>L</sup>&frasl;<sub>D</sub>(U<sup>2</sup>)/<sub>(2g)</sub>
     *
     * @param lambda friction factor
     * @param L pipe length，m
     * @param D pipe diameter，m
     * @param U average velocity，m/s
     * @return h，m
     */
    public double ht(double lambda, double L, double D, double U) {
        return lambda * L / D * U * U / 2 / gravity;
    }

    /**
     *
     *
     * @param Re renolds number
     * @param e roughness m
     * @param D diameter m
     * @return &lambda;-friction factor
     */
    public double calculateLambda(double Re, double D, double e) {
        double temp = 0;
        double epsilon = 2 * e / D;
        double Re1 = 59.7 / Math.pow(epsilon, 8 / 7);
        double Re2 = (665 - 7651 * Math.log10(epsilon)) / epsilon;
        if (Re <= 2000) {
            temp = lambda1(Re);
            System.out.println("Laminar Zone");
        } else if (Re > 2000 && Re < Re1) {
            temp = lambda2(Re);
            System.out.println("Hydraulically Smooth Zone");
        } else if (Re > Re1 && Re < Re2) {
            temp = lambda3(Re, D, e);
            System.out.println("Mixed Friction Zone");
        } else if (Re > Re2) {
            temp = lambda4(D, e);
            System.out.println("Rough Zone");
        }
        return temp;
    }

    /**
     * Laminar Zone
     *
     * @param Re renolds number
     * @return &lambda;-laminar Zone Friction Factor
     */
    public static double lambda1(double Re) {
        return 64.0 / Re;
    }

    /**
     * Hydraulically Smooth Zone
     *
     * @param Re renolds number
     * @return &lambda;-Hydraulically Smooth Zone Friction Factor
     */
    public static double lambda2(double Re) {
        return 0.3164 / Math.pow(Re, 0.25);
    }

    /**
     * Mixed Friction Zone
     *
     * @param Re renolds number
     * @param e roughness m
     * @param D diameter，m
     * @return &lambda;-Mixed Friction Zone Friction Factor
     */
    public static double lambda3(double Re, double D, double e) {
        double a = -1.81 * Math.log10(6.8 / Re + Math.pow(e / 3.7 / D, 1.11));
        return 1 / a / a;
    }

    /**
     * Rough Zone
     *
     * @param e roubhness m
     * @param D diameter，m
     * @return &lambda;-Rough Zone Friction Factor
     */
    public static double lambda4(double D, double e) {
        double epsilon = 2 * e / D;
        return 1.0 / (1.74 - 2 * Math.log10(epsilon)) / (1.74 - 2 * Math.log10(epsilon));
    }

    /**
     * 层流计算公式
     *
     * @param Q volume flowrate，m³/s
     * @param nu kinematic viscosity，㎡/s
     * @param D diameter，m
     * @return &lambda;-Friction Factor
     */
    public static double laminar(double Q, double D, double nu) {
        return 4.15 * Q * nu / Math.pow(D, 4);
    }

    /**
     *
     * @param Re renolds number
     * @param e roughness number m
     * @param D diameter m
     * @param Q volume flowrate m³/s
     * @param nu kinematic viscosity，㎡/s
     * @return h m
     */
    public double RepinzonMethod(double Re, double Q, double D, double e, double nu) {
        double temp = 0;
        double epsilon = 2 * e / D;
        double Re1 = 59.7 / Math.pow(epsilon, 8 / 7);
        double Re2 = (665 - 7651 * Math.log10(epsilon)) / epsilon;
        double beta = 0;
        double m = 1;
        double A;

        if (Re <= 2000) {
            beta = 4.15;
            A = 1;
            m = 1;
            System.out.println("Laminer Zone");
        } else if (Re > 2000 && Re < Re1) {
            A = 0.3164;
            m = 0.25;
            beta = 0.0246;
            System.out.println("Hydraulically Smooth Zone");
        } else if (Re > Re1 && Re < Re2) {
            A = pow(10, 0.127 * Math.log10(e / D) - 0.627);
            m = 0.123;
            beta = 0.0802 * A;
            System.out.println("Mixed Friction Zone");
        } else if (Re > Re2) {
            m = 0;
            beta = 0.0826 * lambda4(D, e);
            System.out.println("Rough Zone");
        }

        temp = beta * pow(nu, m) * pow(Q, 2 - m) / pow(D, 5 - m);
        return temp;
    }

    /**
     *
     * @param length length m
     * @param Type Darcy or Repinzon
     * @return pressure drop kg/(S^2m) Pa=kg/(S^2m^2)
     */
    public double getPressureDrop(double length, Type Type) {
        double ht = this.getHt(length, Type);
        return this.getPressureDrop(ht, density);
    }

    /**
     *
     * @param ht m
     * @param density density kg/m³
     * @return pressure drop kg/(S^2m)
     */
    private double getPressureDrop(double ht, double density) {
        return ht * density * gravity;
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
    public double getDensity() {
        return density;
    }

    /**
     *
     * @return
     */
    public double getDynamicViscosity() {
        return dynamicViscosity;
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
    public double getRoughness() {
        return roughness;
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

}
