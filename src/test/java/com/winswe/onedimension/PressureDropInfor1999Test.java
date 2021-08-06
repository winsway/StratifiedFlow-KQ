/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.onedimension;

import com.winswe.util.Info;
import org.junit.Test;

/**
 *
 * @author winsw
 */
public class PressureDropInfor1999Test {

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
    static double diameter = 0.0243;
    static double densityO = 801, densityW = 1000;
    static double dynamicViscosityW = 1e-3, dynamicViscosityO = 1.6e-3;

    double[] hl = {
        0.491886246,
        0.376697953,
        0.610751636,
        0.494240677,
        0.314474972,
        0.424709266,
        0.675777362,
        0.565386276,
        0.49543657,
        0.274021065,
        0.376880018,
        0.718382828,
        0.614783798,
        0.24505026,
        0.748849991
    };

    @Test
    public void ondDcac() {

        for (int i = 0; i < caseList.length; i++) {
            TwophaseFlowEquation tfe = new TwophaseFlowEquation(
                    caseList[i].Usw, caseList[i].Uso, diameter,
                    densityW, densityO,
                    dynamicViscosityW, dynamicViscosityO,
                    hl[i]);
            System.out.println("Oil Re \t " + tfe.getOilReynoldsNumber() + "\t" + "Water Re \t " + tfe.getWaterReynoldsNumber());

        }

    }

}
