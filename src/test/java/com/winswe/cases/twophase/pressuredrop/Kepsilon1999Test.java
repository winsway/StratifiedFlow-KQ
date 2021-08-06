/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.cases.twophase.pressuredrop;

import com.winswe.io.DataFileWriter;
import com.winswe.solver.TwoPhaseSolverForArray;
import com.winswe.util.Info;
import java.io.File;
import java.io.IOException;
import org.junit.Test;

/**
 *
 * @see Angeli, P. and G.F. Hewitt, Pressure gradient in horizontal
 * liquid–liquid flows. International Journal of Multiphase Flow, 1999. 24(7):
 * p. 1183-1203.
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年4月2日 下午11:10:21
 */
public class Kepsilon1999Test {

    static String position_ = "./tutorials/case/twophase/PressureDrop/1999/kepsilon";
    static Info[] caseList = {
        new Info(0.11, 0.11, 33.52),
        new Info(0.22, 0.11, 72.30),
        new Info(0.11, 0.22, 83.06),
        new Info(0.22, 0.22, 124.67),
        new Info(0.33, 0.11, 118.54),
        new Info(0.33, 0.22, 171.57),
        new Info(0.11, 0.33, 154.47),
        new Info(0.22, 0.33, 183.02),
        new Info(0.33, 0.33, 243.80),
        new Info(0.44, 0.11, 164.16),
        new Info(0.44, 0.22, 228.24),
        new Info(0.11, 0.44, 209.90),
        new Info(0.22, 0.44, 262.87),
        new Info(0.55, 0.11, 238.70),
        new Info(0.11, 0.55, 280.84)
    };

//    @BeforeClass
//    static public void mkdir() {
//        for (int i = 0; i < caseList.length; i++) {
//            File file = new File(position_ + "/", caseList[i].fileName());
//            file.mkdir();
//        }
//    }
    @Test
    public void kepsilon() throws IOException {
        final String position = position_;

        File file = new File(position_ + "/", "readme.txt");
        try ( DataFileWriter des = new DataFileWriter(file)) {
            for (int i = 0; i < caseList.length; i++) {
                Info temp = caseList[i];
                String caseName = temp.fileName();
                TwoPhaseSolverForArray twoPhaseSolver = new TwoPhaseSolverForArray(position, caseName, temp);
                twoPhaseSolver.readConfigure();

                //for mesh
                twoPhaseSolver.createMesh();
                twoPhaseSolver.outPutMesh();

                //for field
                twoPhaseSolver.createField();
                twoPhaseSolver.startIteration();
                twoPhaseSolver.outPutFields();

                des.fileWriter.print(temp.toString() + "\n");
            }
            des.fileWriter.close();
        }
    }
}
