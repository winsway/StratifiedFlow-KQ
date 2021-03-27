package com.winswe.mesh.element;

import com.winswe.math.basic.var.Vector;
import java.util.Arrays;
import java.util.List;

public class Face {

    private int index;
    public final Node[] nodes;
    public final Surface surface;
    /**
     * 位置矢量的权重
     */
    private double gFaceWeight = 1.0;
    /**
     * 数值格式的权重
     */
    private double gf;
    /**
     * 左永远代表内部点
     */
    public Cell left;
    /**
     * 右表示外部，特别是边界处，基于右手定则
     */
    public Cell right;
    /**
     *
     */
    public double maxAbsEigenvalue;

    public Face(
            Node[] nodes,
            Surface surface
    ) {
        this.index = -1;
        this.nodes = nodes;
        this.surface = surface;
    }

    public void setIndex(
            int index
    ) {
        if (this.index == -1) {
            this.index = index;

        } else {
            throw new IllegalStateException("The index can be set only once.");

        }
    }

    public int index() {
        return this.index;
    }

    /**
     * 这个很重要直接得到面的关系
     *
     * @param o 对象
     */
    @Override
    public boolean equals(Object o) {
        //比较的是对象的地址
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        Face otherFace = (Face) o;
        Face thisFace = this;

        boolean same = hasSameNodes(thisFace.nodes, otherFace.nodes);
        if (same) {
            // Set the right nodes from the other face
            thisFace.right = otherFace.left;
            otherFace.right = thisFace.left;

            // Average of the normal from the other face 
            //(subtract since it is pointing in opposite direction)
            Vector avgNormal = thisFace.surface.unitNormal()
                    .sub(otherFace.surface.unitNormal())
                    .unit();
            thisFace.surface.setUnitNormal(avgNormal);
            otherFace.surface.setUnitNormal(avgNormal.mult(-1));
        }
        return same;
    }

    @Override
    public int hashCode() {
        return Arrays.stream(this.nodes)
                .mapToInt(Object::hashCode)
                .sum();
    }

    /**
     *
     *
     * @param a1 点集a1
     * @param a2 点集a2
     * @return 如果此集合包含指定集合中的所有元素，则为true;
     */
    private boolean hasSameNodes(Node[] a1, Node[] a2) {
        List<Node> l1 = Arrays.asList(a1);
        List<Node> l2 = Arrays.asList(a2);

        if (l1.size() != l2.size()) {
            return false;
        }
        return l1.containsAll(l2);
    }

    @Override
    public String toString() {
        return "Face{"
                + ", surface=" + surface
                + '}';
    }

    public double getFaceWeight() {
        return gFaceWeight;
    }

    public void setgFaceWeight(double gFaceWeight) {
        this.gFaceWeight = gFaceWeight;
    }

    public double getGf() {
        return gf;
    }

    public void setGf(double gf) {
        this.gf = gf;
    }

}
