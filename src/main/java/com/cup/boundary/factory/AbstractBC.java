/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.boundary.factory;

import com.cup.boundary.BoundaryCondition;

/**
 *
 * @author winsw
 */
public abstract class AbstractBC implements BoundaryCondition {

    int x, y, z;

    RobinBC robinbc;

    int Type;

    public RobinBC getRobinBC() {
        return robinbc;
    }

    @Override
    public AbstractBC getAbstractBC() {
        return this;
    }

    public void setType(int value, RobinBC robinbc) {
        Type = value;
        robinbc = new RobinBC(robinbc);
    }
}
