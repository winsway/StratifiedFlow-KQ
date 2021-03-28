package com.winswe.field;

import com.alibaba.fastjson.JSONObject;
import com.winswe.boundary.NoSlip;
import com.winswe.boundary.RobinBC;
import com.winswe.boundary.ZeroGradient;
import com.winswe.io.IOobject;
import com.winswe.mesh.Structed2D;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 * 体标量场
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午3:00:55
 */
public class VolScalarField {

    private double[] internalField;

    private List<double[]> boudaryField;

    private final double[] FI;

    private final String name;

    private final Structed2D mesh;

    private final int[] boundaryConditionType;

    private final RobinBC[] boundaryConditionParameter;

    private final int size;

    private final IOobject iOobject;

    private final JSONObject jsonObject;

    public VolScalarField(String name, Structed2D mesh, IOobject iOobject) {
        this.name = name;
        this.mesh = mesh;
        this.iOobject = iOobject;
        this.size = (mesh.getNX() + 2) * (mesh.getNY() + 2);
        this.FI = new double[size];
        this.jsonObject = iOobject.getJsonObject().getJSONObject("Field").getJSONObject(name);

        boundaryConditionType = new int[size];
        boundaryConditionParameter = new RobinBC[size];
        this.setInitialValue();
        this.setOutPut();
    }

    public final void setOutPut() {
        boolean ioField
                = iOobject.getJsonObject().getJSONObject("IO").
                        getJSONObject("Field").getBoolean(name);
        if (ioField) {
            iOobject.getField().add(this);
        }
    }

    public final void setInitialValue() {
        double interValue = jsonObject.getDoubleValue("internal value");
        for (int i = 0; i < FI.length; i++) {
            FI[i] = interValue;

        }

    }

    public void setBoundaryCondition() {
        JSONObject jsonObjectBoundary = jsonObject.getJSONObject("boundary condition");

        String type = jsonObjectBoundary.getString("type");

        System.out.println(this.name + " boundary condition type = " + type);
        RobinBC bc = null;
        if ("noSlip".equals(type)) {
            bc = new NoSlip();

        } else if ("ZeroGradient".equals(type)) {
            bc = new ZeroGradient();

        }
        if (bc != null) {
            try {
                this.loopAllBoundary(1, bc);

            } catch (CloneNotSupportedException ex) {
                Logger.getLogger(VolScalarField.class
                        .getName()).log(Level.SEVERE, null, ex);
            }
        }

    }

    /**
     *
     * @param type boundary type:noSlip =1
     * @param robinBc basic parameter
     * @throws CloneNotSupportedException
     */
    public void loopAllBoundary(int type, RobinBC robinBc) throws CloneNotSupportedException {
        for (int I = 0; I <= mesh.getNX() + 1; I++) {
            int IJS = mesh.getCellIndex(I, 0);
            int IJN = mesh.getCellIndex(I, mesh.getNY() + 1);
            boundaryConditionType[IJS] = type;
            boundaryConditionType[IJN] = type;
            boundaryConditionParameter[IJS] = robinBc.clone();
            boundaryConditionParameter[IJN] = robinBc.clone();
        }
        for (int J = 0; J <= mesh.getNY() + 1; J++) {
            int IJW = mesh.getCellIndex(0, J);
            int IJE = mesh.getCellIndex(mesh.getNX() + 1, J);
            boundaryConditionType[IJW] = type;
            boundaryConditionType[IJE] = type;
            boundaryConditionParameter[IJW] = robinBc.clone();
            boundaryConditionParameter[IJE] = robinBc.clone();
        }
    }

    public double[] getFI() {
        return FI;
    }

    public int[] getBoundaryConditionType() {
        return boundaryConditionType;
    }

    public RobinBC[] getBoundaryConditionParameter() {
        return boundaryConditionParameter;
    }

    public Structed2D getMesh() {
        return mesh;
    }

    public String getName() {
        return name;
    }

}
