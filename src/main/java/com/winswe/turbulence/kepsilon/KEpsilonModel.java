/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.kepsilon;

import com.winswe.field.VolScalarField;
import com.winswe.turbulence.Turbulence;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月25日 下午9:01:09
 */
public class KEpsilonModel implements Turbulence {

    private KEquation kEquation;
    private EpsilonEquation epsilonEquation;

    /**
     * eddy viscosity
     */
    private VolScalarField mut;
    
    

    
}
