/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.mesh;

/**
 * moving mesh interface
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月13日 上午9:38:31
 */
public interface Move {

    /**
     * change mesh
     *
     * @param hl liquid hold up
     */
    public void move(double hl);
}
