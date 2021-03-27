package com.winswe.mesh.element;

import com.winswe.math.basic.var.Vector;

public class Surface {

    private final double area;
    private final Vector centroid;
    private Vector unitNormal;
    private final Vector sf;
    private Vector unitTangent1;
    private Vector unitTangent2;

    public Surface(double area, Vector centroid, Vector unitNormal) {
        this.area = area;
        this.centroid = centroid;
        setUnitNormal(unitNormal);
        sf = unitNormal.mult(area);
    }

    public final void setUnitNormal(Vector unitNormal) {
        this.unitNormal = unitNormal;
        setUnitTangents();
    }

    /**
     *
     * @return 面法向量
     */
    public Vector getSf() {
        return sf;
    }

    public Vector unitNormal() {
        return unitNormal;
    }

    public Vector unitTangent1() {
        return unitTangent1;
    }

    public Vector unitTangent2() {
        return unitTangent2;
    }

    private void setUnitTangents() {
        Vector aVector;
        if (Math.abs(unitNormal.comtValue(0)) < 0.6) {
            aVector = new Vector(1, 0, 0);
        } else if (Math.abs(unitNormal.comtValue(1)) < 0.6) {
            aVector = new Vector(0, 1, 0);
        } else {
            aVector = new Vector(0, 0, 1);
        }

        Vector projection = unitNormal.mult(unitNormal.dotDouble(aVector));
        this.unitTangent1 = aVector.sub(projection).unit();

        this.unitTangent2 = unitNormal.cross(unitTangent1).unit();
    }

    public double getArea() {
        return area;
    }

    public Vector getCentroid() {
        return centroid;
    }

    @Override
    public String toString() {
        return "Surface{"
                + "area=" + area
                + ", centroid=" + centroid
                + ", unitNormal=" + unitNormal
                + '}';
    }
}
