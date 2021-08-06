/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.twophase.flowfield;

import com.winswe.solver.TwoPhaseSolverForOneCase;
import org.junit.Test;

/**
 * lambda 0.5 and Um 0.5 m/s<br>
 * oil superficial velocity 0.25 m/s<br>
 * water superficial velocity 0.25 m/s<br>
 * Qoil 6.1575E-4 m3/s<br>
 * Qwater 6.1575E-4 m3/s<br>
 *
 * @see Kumara, W.A.S., B.M. Halvorsen and M.C. Melaaen, Particle image
 * velocimetry for characterizing the flow structure of oil–water flow in
 * horizontal and slightly inclined pipes. Chemical Engineering Science, 2010.
 * 65(15): p. 4332-4349.
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年4月2日 下午11:10:21
 */
public class FirstTest {

    String position_ = "FirstPoint";

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
