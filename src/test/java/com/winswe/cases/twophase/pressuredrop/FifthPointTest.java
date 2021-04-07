/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.twophase.pressuredrop;

import com.winswe.solver.TwoPhaseSolver;
import org.junit.Test;

/**
 * oil superficial velocity 0.33 m/s<br>
 * water superficial velocity 0.11 m/s<br>
 * Qoil 1.5304402129650454E-4 m3/s<br>
 * Qwater 5.101467376550151E-5 m3/s<br>
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年4月2日 下午11:10:21
 */
public class FifthPointTest {

    String position_ = "FifthPoint";

    @Test
    public void kepsilon() {
        final String position = "./tutorials/case/twophase/" + position_;
        final String caseName = "kepsilon";

        TwoPhaseSolver twoPhaseSolver = new TwoPhaseSolver(position, caseName);
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
        final String position = "./tutorials/case/twophase/" + position_;
        final String caseName = "komega";

        TwoPhaseSolver twoPhaseSolver = new TwoPhaseSolver(position, caseName);
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
        final String position = "./tutorials/case/twophase/" + position_;
        final String caseName = "sstkomega";

        TwoPhaseSolver twoPhaseSolver = new TwoPhaseSolver(position, caseName);
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
        final String position = "./tutorials/case/twophase/" + position_;
        final String caseName = "laminar";

        TwoPhaseSolver twoPhaseSolver = new TwoPhaseSolver(position, caseName);
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
