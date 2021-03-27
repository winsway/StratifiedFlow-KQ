/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.boundary;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月7日 下午6:15:58
 */
public class NoSlip extends RobinBC {

    public NoSlip() {
        super(0, 1.0, 0, 1.0);
    }

}
