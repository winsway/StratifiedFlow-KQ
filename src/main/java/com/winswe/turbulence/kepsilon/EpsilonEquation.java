/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.kepsilon;

import com.winswe.mesh.Structed2D;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月25日 下午9:01:38
 */
public class EpsilonEquation {

    /**
     * von Karman Constant
     */
    final public double CAPPA = 0.41;

    /**
     *
     */
    final public double ELOG = 8.342;

    /**
     *
     * @param mesh mesh
     * @param velocity velocity m/s
     * @param S strain rate tensor
     */
    public void calculateS(
            Structed2D mesh,
            double[] velocity,
            double[] S
    ) {

    }

}
