package com.winswe.solver;

import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.solve.SimpleControl;
import com.winswe.mesh.Structed2D;
import com.winswe.mesh.factory.BipolarCoordianteMesh;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午2:58:23
 */
public class TurbulenceSolver {

    /**
     * Finite volume mesh
     */
    private Structed2D mesh;

    /**
     * velocity
     */
    private VolScalarField U;

    /**
     * effective viscosity
     */
    private VolScalarField mueff;

    /**
     * density
     */
    private VolScalarField density;

    /**
     * configure file
     */
    private final IOobject iOobject;

    /**
     * solve control
     */
    SimpleControl simple;

    /**
     *
     * @param path case path
     * @param caseName
     */
    public TurbulenceSolver(String path, String caseName) {
        iOobject = new IOobject(path, caseName);
        simple = new SimpleControl(iOobject);
    }

    public void readConfigure() {
        System.out.println("read configure file.");
        iOobject.readConfigure();
    }

    public void createMesh() {
        System.out.println("create mesh.");
        mesh = new BipolarCoordianteMesh(iOobject);
        mesh.blockMesh();
    }

    public void outPutMesh() {
        System.out.println("out put mesh file.");
        iOobject.outPutMesh(mesh);
    }

    public void createField() {
        U = new VolScalarField("Velocity", mesh, iOobject);
        U.setBoundaryCondition();
        mueff = new VolScalarField("Viscosity", mesh, iOobject);
        mueff.setBoundaryCondition();
        density = new VolScalarField("Density", mesh, iOobject);
        density.setBoundaryCondition();
    }

    public void outPutFields() {
        System.out.println("out put the fields.");
        iOobject.outPutField();
    }

    /**
     * solve equation
     */
    public void solve() {
        UEquation UEqu = new UEquation(2.2301, mesh, U, mueff, density, iOobject);
        simple.getResult().add(UEqu.getSolverPerformance());
        do {
            UEqu.discrete();
            UEqu.solve();
        } while (simple.loop());

        double temp = 0;
        for (int i = 0; i < U.getFI().length; i++) {
            temp = temp > U.getFI()[i] ? temp : U.getFI()[i];
        }
        System.out.println("max     velocity = " + temp + " m/s");
        System.out.println("average velocity = " + temp / 2.0 + " m/s");
    }

}
