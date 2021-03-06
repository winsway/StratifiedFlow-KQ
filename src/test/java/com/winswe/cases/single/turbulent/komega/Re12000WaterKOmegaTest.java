/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.single.turbulent.komega;

import com.winswe.io.DataFileWriter;
import com.winswe.onedimension.SingleFlowEquation;
import com.winswe.solver.*;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月31日 下午8:49:13
 */
public class Re12000WaterKOmegaTest {

    final String position = "./tutorials/case/single/Re12000Water";
    final String caseName = "komega";

    /**
     * Reynolds nubmer
     */
    final double Re = 12000.0;

    final double density = 1000;

    final double dynamicViscosity = 1e-3;

    final double diameter = 0.0243;

    final double roughness = 1e-4;

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
        File file = new File(position + "/" + caseName + "/", "readme.txt");
        try ( DataFileWriter des = new DataFileWriter(file)) {
            des.fileWriter.format("Volume Flowrate Q       = %.4e m3/s\n", water.getVolumeFlowrate());
            des.fileWriter.format("Average Velocity U      = %.4e m/s\n", water.getAverageVelocity());
            des.fileWriter.format("Re                      = %.2e \n", water.getReynoldsNumber());
            des.fileWriter.format("Dacy dp                 = %9.4f  pa\n", water.getPressureDrop(length, SingleFlowEquation.Type.Darcy));
            des.fileWriter.format("Repinzon dp             = %9.4f  pa\n", water.getPressureDrop(length, SingleFlowEquation.Type.Repinzon));
            des.close();
        }
    }

    public void flowField() {
        System.out.println("readConfigure");

        TwoPhaseSolverForOneCase twoPhaseSolver = new TwoPhaseSolverForOneCase(position, caseName);
        twoPhaseSolver.readConfigure();

        //for mesh
        twoPhaseSolver.createMesh();
        twoPhaseSolver.outPutMesh();

        //for field
        twoPhaseSolver.createField();
        twoPhaseSolver.startIteration();
        twoPhaseSolver.outPutFields();
    }

    public static void main(String[] args) throws IOException {
        Re12000WaterKOmegaTest temp = new Re12000WaterKOmegaTest();
        temp.ondDcac();
        temp.describle();
        temp.flowField();

    }

}
