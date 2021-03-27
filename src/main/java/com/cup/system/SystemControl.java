/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.system;

import com.cup.mesh.Mesh;

/**
 *
 * @author winsway
 */
public class SystemControl {

    public Mesh cell;
    public ControlDict conDict;
    public FvSolution fvSolution;

    public SystemControl(Mesh blockMesh, ControlDict controlDict, FvSolution fvSolution) {
        this.cell = blockMesh;
        this.conDict = controlDict;
        this.fvSolution = fvSolution;
    }

    public SystemControl(Mesh blockMesh, ControlDict controlDict) {
        this.cell = blockMesh;
        this.conDict = controlDict;
    }

    @Override
    public SystemControl clone() throws CloneNotSupportedException {
        SystemControl clone = new SystemControl(cell, conDict, fvSolution);
        clone.cell = this.cell.clone();
        clone.conDict = this.conDict.clone();
        clone.fvSolution = this.fvSolution.clone();
        return clone;
    }
}
