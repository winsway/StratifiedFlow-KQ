package com.winswe.solver;

import com.winswe.turbulence.Turbulence;
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
     * Finite Volume Mesh
     */
    private Structed2D mesh;

    /**
     * Velocity
     */
    private VolScalarField U;

    /**
     * Effective Viscosity
     */
    private VolScalarField mueff;

    /**
     * Dynamic Viscosity
     */
    private VolScalarField mum;

    /**
     * Density
     */
    private VolScalarField density;

    /**
     * Configure File
     */
    private final IOobject iOobject;

    /**
     * Turbulence Model
     */
    Turbulence turbulence;

    /**
     * Solve Control
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
        mueff = new VolScalarField("Effective Viscosity", mesh, iOobject);
        mueff.setBoundaryCondition();
        mum = new VolScalarField("Dynamic Viscosity", mesh, iOobject);
        mum.setBoundaryCondition();
        density = new VolScalarField("Density", mesh, iOobject);
        density.setBoundaryCondition();

        turbulence
                = Turbulence.factory(
                        U, mum, density,
                        mesh,
                        iOobject
                );

    }

    public void outPutFields() {
        System.out.println("out put the fields.");
        iOobject.outPutField();
    }

    /**
     * solve equation
     */
    public void solve() {
        UEquation UEqu = new UEquation(88, mesh, U, mueff, density, iOobject);
        simple.getResult().add(UEqu.getSolverPerformance());
        turbulence.setSolverPerformance(simple.getResult());
//        do {
//            UEqu.discrete();
//            UEqu.solve();
//        } while (simple.loop());
//        turbulence.setTurPar();
////        iOobject.outPutField();

        do {
            UEqu.discrete();
            UEqu.solve();
            if (simple.getCount() == 0) {
                turbulence.setTurPar();
            }
//            iOobject.outPutField();
            turbulence.solve();
            this.modifiedViscosity();

//            iOobject.outPutField();
        } while (simple.loop());

        double temp = 0;
        for (int i = 0; i < U.getFI().length; i++) {
            temp = temp > U.getFI()[i] ? temp : U.getFI()[i];
        }
        System.out.println("max     velocity = " + temp + " m/s");
        System.out.println("average velocity = " + temp / 2.0 + " m/s");
    }

    public void modifiedViscosity() {
        turbulence.updateMut();
        int IJ;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);
                mueff.getFI()[IJ]
                        = (mum.getFI()[IJ]
                        + 0.20 * turbulence.getEddyViscosity().getFI()[IJ]);
            }
        }
    }

}
