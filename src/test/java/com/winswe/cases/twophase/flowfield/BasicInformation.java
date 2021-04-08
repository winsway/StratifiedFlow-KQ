/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.twophase.flowfield;

import com.winswe.onedimension.SingleFlowEquation;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年4月2日 下午10:46:03
 */
public class BasicInformation {

    final String position = "./winsway/";
    final String case_ = "Re7777Water";

    /**
     * Reynolds nubmer
     */
    final double Re = 100.0;

    final double density = 1000;

    final double dynamicViscosity = 1e-3;

    final double diameter = 0.056;

    final double roughness = 8e-5;

    final double length = 1.0;

    SingleFlowEquation water;

    double Uoil = 0.3225;
    double Uwater = 0.105;

    public void ondDcac() {

        double flowRateOil
                = SingleFlowEquation.calculateVolumeFlowrateByAverageVelocity(
                        Uoil,
                        diameter);
        double flowRateWater
                = SingleFlowEquation.calculateVolumeFlowrateByAverageVelocity(
                        Uwater,
                        diameter);

        System.out.println("Oil Volume Flow Rate = " + flowRateOil + " m3/s");
        System.out.println("Water Volume Flow Rate = " + flowRateWater + " m3/s");

        double re = SingleFlowEquation.calculateReynoldsNumberByDynamicViscosity(
                flowRateWater,
                diameter,
                density,
                dynamicViscosity
        );

        System.out.println("Re = " + re);
    }

    public static void main(String[] args) {
        BasicInformation temp = new BasicInformation();
        temp.ondDcac();
    }
}
