/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.kepsilon;

import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.solve.SolverPerformance;
import com.winswe.mesh.Structed2D;
import com.winswe.turbulence.Turbulence;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.List;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年3月25日 下午9:01:09
 */
public class KEpsilonModel extends Turbulence {

    private final KEquation kEquation;
    private final EpsilonEquation epsilonEquation;

    private final VolScalarField k;
    private final VolScalarField epsilon;
    private final VolScalarField S;

    protected static double sigmak;
    protected static double sigmae;
    protected static double Cu;
    protected static double Cu25;
    protected static double Cu75;
    protected static double CTRANS;
    protected static double Ce1;
    protected static double Ce2;

    public KEpsilonModel(
            VolScalarField velocity,
            VolScalarField dynamicViscosity,
            VolScalarField density,
            Structed2D mesh,
            IOobject iOobject
    ) {
        super(velocity, dynamicViscosity, density, mesh, iOobject);
        k = new VolScalarField("K", mesh, iOobject);
        epsilon = new VolScalarField("Epsilon", mesh, iOobject);
        S = new VolScalarField("S", mesh, iOobject);
        sigmak = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmak");
        sigmae = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmae");
        Cu = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("Cu");
        Cu25 = sqrt(sqrt(Cu));
        Cu75 = pow(Cu25, 3.0);
        CTRANS = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("CTRANS");
        Ce1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("Ce1");
        Ce2 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("Ce2");

        kEquation = new KEquation(
                mesh,
                k,
                epsilon,
                this.velocity,
                this.dynamicViscosity,
                mut,
                density,
                iOobject,
                S);

        epsilonEquation = new EpsilonEquation(
                mesh,
                k,
                epsilon,
                this.velocity,
                dynamicViscosity,
                mut,
                density,
                iOobject,
                S);
    }

    @Override
    public void solve() {

        kEquation.discrete();
        kEquation.solve();
        epsilonEquation.discrete();
        epsilonEquation.solve();
    }

    /**
     *
     * @param K Turbulent Kinetic Energy
     * @return Cu^0.25 * sqrt(max(0, K));
     */
    static public double getCK(double K) {
        return Cu25 * sqrt(max(0, K));
    }

    @Override
    public void setTurPar() {
        int IJ;
        double diameter
                = iOobject.getJsonObject().
                        getJSONObject("geometric")
                        .getDoubleValue("diameter");
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);
                k.getFI()[IJ] = setK(velocity.getFI()[IJ]);
                epsilon.getFI()[IJ] = setEpsilon(k.getFI()[IJ], diameter);
                mut.getFI()[IJ]
                        = cacMut(
                                density.getFI()[IJ],
                                k.getFI()[IJ],
                                epsilon.getFI()[IJ]
                        );
            }
        }
        this.calculateS(mesh, velocity, S);
    }

    @Override
    public void updateMut() {
        int IJ;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);
                mut.getFI()[IJ]
                        = cacMut(
                                density.getFI()[IJ],
                                k.getFI()[IJ],
                                epsilon.getFI()[IJ]
                        );
            }
        }
    }

    public double cacMut(double rho, double K, double epsilon) {
        return Cu * rho * K * K / (epsilon + 1e-30);
    }

    private double setK(double velocity) {
        return 1e-4 * pow(velocity, 2.0);
    }

    /**
     *
     * @param k
     * @param diameter
     * @return
     */
    public double setEpsilon(double k, double diameter) {
        return pow(0.09, 0.75) * pow(k, 3.0 / 2.0)
                / turbulenceLengthScale(diameter);
    }

    @Override
    public void setSolverPerformance(
            List<SolverPerformance> solverPerformances
    ) {
        solverPerformances.add(kEquation.getSolverPerformance());
        solverPerformances.add(epsilonEquation.getSolverPerformance());
    }

}
