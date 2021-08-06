/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.gridindependence;

import com.winswe.io.DataFileWriter;
import com.winswe.onedimension.SingleFlowEquation;
import com.winswe.solver.TwoPhaseSolverForOneCase;
import java.io.File;
import java.io.IOException;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 *
 * oil superficial velocity 0.14 m/s<br>
 * water superficial velocity 0.318 m/s<br>
 * Qoil = 7.0939e-05 m3/s<br>
 * Qwater = 1.6113e-04 m3/s<br>
 *
 * @see Yusuf, N., et al., Effect of oil viscosity on the flow structure and
 * pressure gradient in horizontal oil–water flow. Chemical Engineering Research
 * and Design, 2012. 90(8): p. 1019-1030.
 *
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月31日 下午8:49:13
 */
public class GridIndependenceSSTKomegaTest {

    static final String position_ = "./tutorials/case/grid independence/sstkomega";

    @BeforeClass
    static public void basicDiscribtion() throws IOException {
        System.out.println("grid independence test");
        File file = new File(position_ + "/", "readme.txt");
        double waterVolumeFlowRate = SingleFlowEquation.calculateVolumeFlowrateByAverageVelocity(0.318, 0.0254);
        double oilVolumeFlowRate = SingleFlowEquation.calculateVolumeFlowrateByAverageVelocity(0.14, 0.0254);
        try ( DataFileWriter des = new DataFileWriter(file)) {
            des.fileWriter.format("Volume Flowrate of Water Q       = %.4e m3/s\n", waterVolumeFlowRate);
            des.fileWriter.format("Average Velocity of Water U      = %.4e m/s\n", 0.318);
            des.fileWriter.format("Density of Water                 = %.4f kg/m3\n", 998.);
            des.fileWriter.format("Viscosity of Water               = %.4e Pa.s\n", 1E-3);
            des.fileWriter.print("\n");
            des.fileWriter.format("Volume Flowrate of Oil Q       = %.4e m3/s\n", oilVolumeFlowRate);
            des.fileWriter.format("Average Velocity of Oil U      = %.4e m/s\n", 0.318);
            des.fileWriter.format("Density of Oil                 = %.4f kg/m3\n", 875.);
            des.fileWriter.format("Viscosity of Oil               = %.4e Pa.s\n", 12E-3);
            des.close();
        }
    }

    @Test
    public void sstkomega20X20() {
        final String position = position_ + "/";
        final String caseName = "20X20";

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

    @Test
    public void sstkomega50X50() {
        final String position = position_ + "/";
        final String caseName = "50X50";

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

    @Test
    public void sstkomega80X80() {
        final String position = position_ + "/";
        final String caseName = "80X80";

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

    @Test
    public void sstkomega100X100() {
        final String position = position_ + "/";
        final String caseName = "100X100";

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
}
