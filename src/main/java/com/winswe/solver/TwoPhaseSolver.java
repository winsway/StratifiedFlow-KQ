package com.winswe.solver;

import com.alibaba.fastjson.JSONObject;
import com.winswe.turbulence.Turbulence;
import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.solve.SimpleControl;
import com.winswe.mesh.Structed2D;
import com.winswe.mesh.factory.BipolarCoordianteMesh;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.abs;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午2:58:23
 */
public class TwoPhaseSolver {

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
     * Density
     */
    private VolScalarField phi;

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
    public TwoPhaseSolver(String path, String caseName) {
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
        phi = new VolScalarField("Phi", mesh, iOobject);
        phi.setBoundaryCondition();
        mueff = new VolScalarField("Effective Viscosity", mesh, iOobject);
        mueff.setBoundaryCondition();
        mum = new VolScalarField("Dynamic Viscosity", mesh, iOobject);
        mum.setBoundaryCondition();
        density = new VolScalarField("Density", mesh, iOobject);
        density.setBoundaryCondition();

        JSONObject jSONObject = iOobject.getJsonObject().getJSONObject("transportProperties");
        Qoil = jSONObject.getJSONObject("Oil").getDoubleValue("FlowRate");
        Qwater = jSONObject.getJSONObject("Water").getDoubleValue("FlowRate");

        int IJ;
        for (int Y = 0; Y <= mesh.getNY() + 1; ++Y) {
            for (int X = 0; X <= mesh.getNX() + 1; ++X) {
                IJ = mesh.getCellIndex(X, Y);

                if (Y <= mesh.getNY() / 2) {
                    phi.getFI()[IJ] = -1;
                    density.getFI()[IJ] = jSONObject.getJSONObject("Oil").getDoubleValue("Density");
                    mum.getFI()[IJ] = jSONObject.getJSONObject("Oil").getDoubleValue("Viscosity");
                    mueff.getFI()[IJ] = jSONObject.getJSONObject("Oil").getDoubleValue("Viscosity");
                } else {
                    phi.getFI()[IJ] = 1;
                    density.getFI()[IJ] = jSONObject.getJSONObject("Water").getDoubleValue("Density");
                    mum.getFI()[IJ] = jSONObject.getJSONObject("Water").getDoubleValue("Viscosity");
                    mueff.getFI()[IJ] = jSONObject.getJSONObject("Water").getDoubleValue("Viscosity");
                }
            }
        }

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
     *
     * @param pressureDrop
     * @param hl liquid level
     */
    public void solve(double pressureDrop, double hl) {
        simple = new SimpleControl(iOobject);
        mesh.move(hl);

        UEquation UEqu
                = new UEquation(
                        pressureDrop,
                        mesh,
                        U,
                        mueff,
                        density,
                        iOobject
                );

        simple.getResult().add(UEqu.getSolverPerformance());

        if (turbulence != null) {
            turbulence.setSolverPerformance(simple.getResult());
        }

        do {
            UEqu.discrete();
            UEqu.solve();

            if (turbulence != null) {
                if (simple.getCount() == 0) {
                    turbulence.setTurPar();
                }

                turbulence.solve();
                this.modifiedViscosity();
            }

        } while (simple.loop());

        double temp = 0;
        for (int i = 0; i < U.getFI().length; i++) {
            temp = temp > U.getFI()[i] ? temp : U.getFI()[i];
        }
        System.out.println("max     velocity = " + temp + " m/s");
        System.out.println("average velocity = " + temp / 2.0 + " m/s");
    }

    public void modifiedViscosity() {
        double gamma
                = iOobject.getJsonObject().
                        getJSONObject("turbulence").
                        getDoubleValue("gamma");

        turbulence.updateMut();
        int IJ;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);
                mueff.getFI()[IJ]
                        = (mum.getFI()[IJ]
                        + gamma * turbulence.getEddyViscosity().getFI()[IJ]);
            }
        }
    }

    private double oilFlowRate, waterFlowRate;

    public void getFlowRate() {
        oilFlowRate = 0;
        waterFlowRate = 0;
        int IJ;
        double volume;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);
                volume = mesh.getVolume(X, Y);

                if (phi.getFI()[IJ] == 1) {
                    waterFlowRate += U.getFI()[IJ] * volume;
                }

                if (phi.getFI()[IJ] == -1) {
                    oilFlowRate += U.getFI()[IJ] * volume;
                }

            }
        }
    }

    private double Qoil, Qwater;

    public void startIteration() {

        double errorWater, errorOil;
        double xold = 1, yold = 0.5;
        double xnew, ynew;
        double water, oil, waterY, waterX, oilY, oilX;
        PrintWriter pw;

        String dirName = iOobject.getPath() + "/" + iOobject.getCaseName();
        String title = "pressure and liquid high";
        String fileName = dirName + "/" + title + "." + "txt";

        try {
            pw = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
            pw.printf("Qwater   = %.4e\t  Qoil = %.4e\t\n", Qwater, Qoil);
            pw.println();
            pw.flush();
            do {

                this.solve(xold, yold);
                this.getFlowRate();
                oil = this.oilError(Qoil);
                water = this.waterError(Qwater);

                pw.printf("F dpdz   = %.4e\t  hl   = %.4e\t\n", xold, yold);
                pw.printf("Qwater   = %.4e\t  Qoil = %.4e\t\n", waterFlowRate, oilFlowRate);
                pw.flush();

                this.solve(xold * 1.01, yold);
                this.getFlowRate();
                oilX = (this.oilError(Qoil) - oil) / (0.01 * xold);
                waterX = (this.waterError(Qwater) - water) / (0.01 * xold);

                pw.printf("Fx dpdz  = %.4e\t  hl    = %.4e\t\n", xold * 1.01, yold);
                pw.printf("Qwater   = %.4e\t  Qoil  = %.4e\t\n", waterFlowRate, oilFlowRate);
                pw.flush();

                this.solve(xold, yold * 1.01);
                this.getFlowRate();
                oilY = (this.oilError(Qoil) - oil) / (0.01 * yold);
                waterY = (this.waterError(Qwater) - water) / (0.01 * yold);

                pw.printf("Fy dpdz  = %.4e\t  hl    = %.4e\t\n", xold, yold * 1.01);
                pw.printf("Qwater   = %.4e\t  Qoil  = %.4e\t\n", waterFlowRate, oilFlowRate);
                pw.flush();
                pw.println();

                xnew
                        = xold + (oil * waterY - water * oilY)
                        / (waterX * oilY - oilX * waterY);
                ynew
                        = yold + (water * oilX - oil * waterX)
                        / (waterX * oilY - oilX * waterY);
//          
                xnew = Modified(xnew, xold, 1e-3, 10000000);
                ynew = Modified(ynew, yold, 0.1, 0.9);

                xold = xnew;
                yold = ynew;
//
                errorWater = abs((water) / Qwater);
                errorOil = abs((oil) / Qoil);

            } while (errorWater >= 1e-5 || errorOil >= 1e-5);

            pw.close();

        } catch (FileNotFoundException ex) {
            Logger.getLogger(TwoPhaseSolver.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public double Modified(
            double value,
            double old,
            double minClip,
            double maxClip
    ) {
        double temp;
        temp
                = value < minClip ? old * 0.8
                        : value > maxClip ? old * 1.05
                                : value;
        return temp;
    }

    double waterError(double value) {
        return value - waterFlowRate;
    }

    double oilError(double value) {
        return value - oilFlowRate;
    }

}
