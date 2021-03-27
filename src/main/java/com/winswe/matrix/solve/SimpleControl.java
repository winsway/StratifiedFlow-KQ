/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.matrix.solve;

import com.winswe.io.IOobject;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月24日 下午3:42:49
 */
public class SimpleControl {

    /**
     * configure file
     */
    private final IOobject iOobject;

    /**
     * max iteration number
     */
    private final int maxNubmer = 10000;

    private int count;

    private List<SolverPerformance> result = new ArrayList<>();

    public SimpleControl(IOobject iOobject) {
        this.iOobject = iOobject;

    }

    public boolean loop() {
        boolean check = true;

        for (SolverPerformance solverPerformance : result) {
            if (!solverPerformance.convergence()) {
                check = false;
                break;
            }
        }

        if (count > maxNubmer) {
            check = true;
        }

        this.next();

        this.output();

        return !check;
    }

    public void next() {
        count++;
    }

    public void output() {
        System.out.printf("Loop Number = %d\n", this.count);
        for (SolverPerformance solverPerformance : result) {
            System.out.printf(
                    "Equation Name = %s, "
                    + "First Residual = %.4e, "
                    + "Relative Residual = %.4e, "
                    + "Convergence Criterion = %.4e, "
                    + "Sub Loop Number = %d\n",
                    solverPerformance.getName(),
                    solverPerformance.getFirstResidual(),
                    solverPerformance.getRelativeResidual(),
                    solverPerformance.getConvergenceCriterion(),
                    solverPerformance.getSubLoopNumber()
            );
        }
        System.out.println("");
    }

    public IOobject getiOobject() {
        return iOobject;
    }

    public int getMaxNubmer() {
        return maxNubmer;
    }

    public int getCount() {
        return count;
    }

    public List<SolverPerformance> getResult() {
        return result;
    }

}
