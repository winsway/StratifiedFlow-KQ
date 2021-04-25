/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.onedimension;

import com.winswe.io.DataFileWriter;

import java.io.File;
import java.io.IOException;

/**
 * 用于测试层流求解算法和解析解的区别
 *
 * @author winsway
 */
public class Laminar {

    final String position = "./tutorials/";
    final String case_ = "Re100Water";

    /**
     * Reynolds nubmer
     */
    final double Re = 1000.0;

    final double density = 1000;

    final double dynamicViscosity = 1e-3;

    final double diameter = 0.05;

    final double roughness = 2e-4;

    final double length = 1.0;

    SingleFlowEquation water;

    public void ondDcac() {
        double volumeFlowrate = SingleFlowEquation.calculateVolumeFlowrateByReynoldsNumber(
                Re,
                diameter,
                density,
                dynamicViscosity
        );
        System.out.println("volume flowrate  = " + volumeFlowrate + " m3/s");
        water = new SingleFlowEquation(volumeFlowrate, diameter, density, dynamicViscosity, roughness);
    }

    /**
     *
     *
     * @throws java.io.IOException
     */
    public void describle() throws IOException {
        File file = new File(position + case_ + "/", "readme.txt");
        try ( DataFileWriter des = new DataFileWriter(file)) {
            des.fileWriter.format("Volume Flowrate Q       = %.4e m3/s\n", water.getVolumeFlowrate());
            des.fileWriter.format("Average Velocity U      = %.4e m/s\n", water.getAverageVelocity());
            des.fileWriter.format("Re                      = %.2e \n", water.getReynoldsNumber());
            des.fileWriter.format("Dacy dp                 = %9.4f  pa\n", water.getPressureDrop(length, SingleFlowEquation.Type.Darcy));
            des.fileWriter.format("Repinzon dp             = %9.4f  pa\n", water.getPressureDrop(length, SingleFlowEquation.Type.Repinzon));
            des.close();
        }
    }

    public static void main(String[] args) throws IOException {
        Laminar test = new Laminar();
        test.ondDcac();
        test.describle();
    }
}
