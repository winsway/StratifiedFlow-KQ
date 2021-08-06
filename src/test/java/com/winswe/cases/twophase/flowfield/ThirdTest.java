/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.twophase.flowfield;

import com.winswe.solver.TwoPhaseSolverForOneCase;
import org.junit.Test;

/**
 * lambda 0.25 and Um 0.43 m/s<br>
 * oil superficial velocity 0.3225 m/s<br>
 * water superficial velocity 0.105 m/s<br>
 * Qoil 7.9432E-4 m3/s<br>
 * Qwater 2.5861E-4 m3/s<br>
 *
 * @see Kumara, W.A.S., et al., Comparison of Particle Image Velocimetry and
 * Laser Doppler Anemometry measurement methods applied to the oil–water flow in
 * horizontal pipe. Flow Measurement and Instrumentation, 2010. 21(2): p.
 * 105-117.
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年4月2日 下午11:10:21
 */
public class ThirdTest {

    String position_ = "ThirdPoint";

    @Test
    public void kepsilon() {
        final String position = "./tutorials/case/flowfield/" + position_;
        final String caseName = "kepsilon";

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
    public void komega() {
        final String position = "./tutorials/case/flowfield/" + position_;
        final String caseName = "komega";

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
    public void sstkomega() {
        final String position = "./tutorials/case/flowfield/" + position_;
        final String caseName = "sstkomega";

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
    public void laminar() {
        final String position = "./tutorials/case/flowfield/" + position_;
        final String caseName = "laminar";

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
