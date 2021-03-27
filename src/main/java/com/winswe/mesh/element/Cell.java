package com.winswe.mesh.element;

import com.winswe.finiteVolume.geom.GeometryHelper;
import com.winswe.finiteVolume.geom.GeometryHelper.TriGeom;
import com.winswe.finiteVolume.geom.Point;
import com.winswe.math.basic.var.Vector;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class Cell {

    /**
     * This value must be equal to the index in Mesh.cells() List. This value
     * must be equal to -1 for ghost cells.
     */
    private int index;
    private int IZONE;
    private Node[] nodes;
    private ArrayList<Face> faces;
    private Shape shape;

    /**
     * Local time step or pseudo-time step
     */
    public double dt;

    public Cell() {
        this.index = -1;
        this.IZONE = -1;
        this.faces = new ArrayList<>();
    }

    public Cell(Node[] nodes, Shape shape) {
        this.nodes = nodes;
        this.shape = shape;
//        
        this.index = -1;
        this.IZONE = -1;
        this.faces = new ArrayList<>();
    }

    public int index() {
        return index;
    }

    public void setIndex(int index) {
        if (this.index != -1) {
            throw new IllegalStateException("The index of a cell can be set only once.");
        }
        this.index = index;
    }

    public int Izone() {
        return IZONE;
    }

    public void setIzone(int Izone) {
        if (this.IZONE != -1) {
            throw new IllegalStateException("The index of a cell can be set only once.");
        }
        this.IZONE = Izone;
    }

    /**
     * 设定形状
     */
    public void setShape() {
//      分割面和三角形
        List<TriGeom> surfaceTriangles = new ArrayList<>();
//
        for (Face face : faces) {
            //得到点集  
//            <editor-fold>得到一个面点集，且均为逆时针方向  
            Point[] po = new Point[face.nodes.length];
            //需要判断方向：顺时针/逆时针
            if (face.left.index != index) {
                for (int I = 0; I < face.nodes.length; I++) {
                    po[I]
                            = new Point(
                                    face.nodes[face.nodes.length - 1 - I]
                                            .location()
                            );
                }
            } else {
                for (int I = 0; I < face.nodes.length; I++) {
                    po[I] = new Point(face.nodes[I].location());
                }
            }
//               </editor-fold>        
            //分割三角形
            Point p0 = po[0];
            for (int i = 2; i < po.length; i++) {
                Point p1 = po[i - 1];
                Point p2 = po[i];
                TriGeom e = new TriGeom(p0, p1, p2);
                surfaceTriangles.add(e);
            }
        }

        TriGeom[] geo = surfaceTriangles.toArray(new TriGeom[surfaceTriangles.size()]);
        double volume = GeometryHelper.volume(geo);
        Vector centroid = GeometryHelper.centroid(geo).toVector();

        this.shape = new Shape(volume, centroid);
    }

    /**
     * 设定点，排列顺序没有要求。
     */
    public void setNodes() {
        List<Node> no = new ArrayList<>();
        for (Face face : faces) {
            for (Node e : face.nodes) {
                no.add(e);
            }
        }
        no = no.stream().distinct().collect(Collectors.toList());
        this.nodes = no.toArray(new Node[no.size()]);
    }

    @Override
    public String toString() {
        return "Cell{"
                + "\nindex=" + index
                + "\nnodes=" + (nodes == null ? null : Arrays.asList(nodes))
                + "\nfaces=" + faces
                + "\nshape=" + shape
                + "\n}";
    }

    @Override
    public boolean equals(Object o) {
        //比较的是对象的地址
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        Cell othercCell = (Cell) o;
        Cell thisCell = this;

        boolean same = hasSameNodes(thisCell.nodes, othercCell.nodes);
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

    public ArrayList<Face> getFaces() {
        return faces;
    }

    public Node[] getNodes() {
        return nodes;
    }

    public Shape getShape() {
        return shape;
    }

}
