/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.turbulence.sstkomega;

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
public class SSTKOmegaModel extends Turbulence {

    private final KEquation kEquation;
    private final OmegaEquation omegaEquation;

    private final VolScalarField k;
    private final VolScalarField omega;
    private final VolScalarField S;
    private final VolScalarField DOMEGA;
    private final VolScalarField F1SST;
    private final VolScalarField distance;
    private final VolScalarField GEN;

    protected static double a1;
    protected static double alpha1;
    protected static double alpha2;
    protected static double beta1;
    protected static double beta2;
    protected static double sigmak1;
    protected static double sigmak2;
    protected static double sigmaOmega1;
    protected static double sigmaOmega2;
    protected static double Cu;
    protected static double Cu25;
    protected static double Cu75;
    protected static double CTRANS;

    protected static double betaStar;

    public SSTKOmegaModel(
            VolScalarField velocity,
            VolScalarField dynamicViscosity,
            VolScalarField density,
            Structed2D mesh,
            IOobject iOobject
    ) {
        super(velocity, dynamicViscosity, density, mesh, iOobject);
        k = new VolScalarField("K", mesh, iOobject);
        omega = new VolScalarField("Omega", mesh, iOobject);
        S = new VolScalarField("S", mesh, iOobject);
        DOMEGA = new VolScalarField("DOMEGA", mesh, iOobject);
        F1SST = new VolScalarField("F1SST", mesh, iOobject);
        distance = new VolScalarField("distance", mesh, iOobject);
        GEN = new VolScalarField("GEN", mesh, iOobject);

        a1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("a1");

        alpha1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("alpha1");

        alpha2 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("alpha2");

        beta1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("beta1");

        beta2 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("beta2");

        sigmak1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmak1");

        sigmak2 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmak2");

        sigmaOmega1 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmaOmega1");

        sigmaOmega2 = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("sigmaOmega2");

        betaStar = iOobject.
                getJsonObject().
                getJSONObject("turbulence").
                getDoubleValue("betaStar");

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

        kEquation = new KEquation(
                mesh,
                k,
                omega,
                this.velocity,
                this.dynamicViscosity,
                mut,
                density,
                iOobject,
                S,
                F1SST,
                GEN);

        omegaEquation = new OmegaEquation(
                mesh,
                k,
                omega,
                this.velocity,
                this.dynamicViscosity,
                mut,
                density,
                iOobject,
                S,
                F1SST,
                DOMEGA,
                GEN);

    }

    @Override
    public void solve() {

        this.calculateS(mesh, velocity, S);

        this.sstBlendingFunction();

        kEquation.discrete();
        kEquation.solve();

        omegaEquation.discrete();
        omegaEquation.solve();

        this.updateBound();

    }

    public void updateBound() {

        for (int J = 0; J <= mesh.getNY() + 1; J++) {
            int IJW = mesh.getCellIndex(1, J);
            omega.getFI()[mesh.getCellIndex(0, J)]
                    = omega.getFI()[IJW];

            int IJE = mesh.getCellIndex(mesh.getNX(), J);
            omega.getFI()[mesh.getCellIndex(mesh.getNX() + 1, J)]
                    = omega.getFI()[IJE];
        }

        //        
        for (int I = 0; I <= mesh.getNX() + 1; I++) {
            int IJS = mesh.getCellIndex(I, 1);
            omega.getFI()[mesh.getCellIndex(I, 0)]
                    = omega.getFI()[IJS];

            int IJN = mesh.getCellIndex(I, mesh.getNY());
            omega.getFI()[mesh.getCellIndex(mesh.getNX() + 1, 0)]
                    = omega.getFI()[IJN];
        }

    }

    private void sstBlendingFunction() {

        cacDKW(k, omega, mesh, DOMEGA);

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {

                int IJ = mesh.getCellIndex(X, Y);

                double CDKW = getCDKW(
                        density.getFI()[IJ],
                        omega.getFI()[IJ],
                        DOMEGA.getFI()[IJ]
                );

                double KSI = getKSI(
                        k.getFI()[IJ],
                        omega.getFI()[IJ],
                        dynamicViscosity.getFI()[IJ],
                        density.getFI()[IJ],
                        distance.getFI()[IJ],
                        CDKW);

                F1SST.getFI()[IJ] = getF1(KSI);

            }

        }

    }

    /**
     *
     * @param KSI
     * @return
     */
    private double getF1(double KSI) {
        return Math.tanh(pow(KSI, 4.0));
    }

    private double getKSI(
            double K,
            double Omega,
            double mu,
            double rho,
            double DN,
            double CDKW
    ) {
        double A, B, C;
        A = Math.sqrt(K) / (betaStar * Omega * DN);
        B = (500.0 * mu) / (rho * Omega * DN * DN);
        C = 4.0 * rho * sigmaOmega2 * K / (CDKW * DN * DN);
        return Math.min(Math.max(A, B), C);
    }

    private double getCDKW(double rho, double Omega, double DKW) {
        return Math.max(
                2.0 * rho * sigmaOmega2 / Omega * DKW,
                1.E-10);
    }

    /**
     * phi1 * F1 + phi2 * (1.0 - F1)
     *
     * @param phi1
     * @param phi2
     * @param F1
     * @return
     */
    public static double getPhi(double phi1, double phi2, double F1) {
        return phi1 * F1 + phi2 * (1.0 - F1);
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

        this.cacDistance(mesh, distance);

        int IJ;
        double diameter
                = iOobject.getJsonObject().
                        getJSONObject("geometric")
                        .getDoubleValue("diameter");

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);

                k.getFI()[IJ] = setK(velocity.getFI()[IJ]);

                omega.getFI()[IJ] = setOmega(k.getFI()[IJ], diameter);

                mut.getFI()[IJ]
                        = cacMut(
                                density.getFI()[IJ],
                                k.getFI()[IJ],
                                omega.getFI()[IJ],
                                0,
                                0
                        );
                F1SST.getFI()[IJ] = 1.0;
            }
        }

    }

    @Override
    public void updateMut() {
        int IJ;
        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                IJ = mesh.getCellIndex(X, Y);

                double F2
                        = getF2(
                                k.getFI()[IJ],
                                omega.getFI()[IJ],
                                dynamicViscosity.getFI()[IJ],
                                density.getFI()[IJ],
                                distance.getFI()[IJ]
                        );

                mut.getFI()[IJ]
                        = cacMut(
                                density.getFI()[IJ],
                                k.getFI()[IJ],
                                omega.getFI()[IJ],
                                sqrt(S.getFI()[IJ]),
                                F2
                        );
            }
        }
    }

    public double getF2(
            double K,
            double Omega,
            double mu,
            double rho,
            double DN
    ) {
        double A, B;
        A = 2.0 * Math.sqrt(K) / (betaStar * DN * Omega);
        B = (500.0 * mu / rho) / (DN * DN * Omega);
        return Math.tanh(pow(Math.max(A, B), 2.0));
    }

    /**
     *
     * @param rho density
     * @param K turbulence kinetic energy
     * @param omega specific dissipation rate
     * @param smod the magnitude of the strain rate
     * @param F2 blend factor
     * @return eddy viscosity
     */
    public double cacMut(
            double rho,
            double K,
            double omega,
            double smod,
            double F2
    ) {
        return a1 * rho * K
                / Math.max(a1 * omega, smod * F2);
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
    public double setOmega(double k, double diameter) {
        return pow(0.09, 0.75) * pow(k, 3.0 / 2.0)
                / turbulenceLengthScale(diameter)
                / k;
    }

    @Override
    public void setSolverPerformance(
            List<SolverPerformance> solverPerformances
    ) {
        solverPerformances.add(kEquation.getSolverPerformance());
        solverPerformances.add(omegaEquation.getSolverPerformance());
    }

    /**
     * calculate cross term
     *
     * @param k turbulent kinetic energy
     * @param omega specific dissipation rate
     * @param mesh mesh
     * @param DOMEGA cross diffusion term
     */
    private void cacDKW(
            VolScalarField k,
            VolScalarField omega,
            Structed2D mesh,
            VolScalarField DOMEGA
    ) {
        double alpha;
        double a, b, c, d;

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                int IJ = mesh.getCellIndex(X, Y);
                int IJW = mesh.getCellIndex(X - 1, Y);
                int IJE = mesh.getCellIndex(X + 1, Y);
                int IJS = mesh.getCellIndex(X, Y - 1);
                int IJN = mesh.getCellIndex(X, Y + 1);

                alpha = mesh.alpha(mesh.getPointX1()[X], mesh.getPointX2()[Y]);

                a = (k.getFI()[IJE] - k.getFI()[IJW])
                        / (mesh.getDXP()[X - 1] + mesh.getDXP()[X]);

                b = (omega.getFI()[IJE] - omega.getFI()[IJW])
                        / (mesh.getDXP()[X - 1] + mesh.getDXP()[X]);

                c = (k.getFI()[IJN] - k.getFI()[IJS])
                        / (mesh.getDYP()[Y - 1] + mesh.getDYP()[Y]);

                d = (omega.getFI()[IJN] - omega.getFI()[IJS])
                        / (mesh.getDYP()[Y - 1] + mesh.getDYP()[Y]);

                DOMEGA.getFI()[IJ] = Math.max(0, (a * b + c * d) / alpha);

            }
        }
    }

    /**
     *
     * @param mesh
     * @param distance
     */
    private void cacDistance(
            Structed2D mesh,
            VolScalarField distance
    ) {
        double temp0, temp1;
        double aw, ae, as, an;

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {

                double aw1 = 1E40, ae1 = 1E40, as1 = 1E40, an1 = 1E40;

                int IJ = mesh.getCellIndex(X, Y);
                //
                for (int J = 0; J <= mesh.getNY() + 1; J++) {
                    aw = distanceBetweenTwoPoints(
                            mesh.X(mesh.getPointX1()[0], mesh.getPointX2()[J]),
                            mesh.Y(mesh.getPointX1()[0], mesh.getPointX2()[J]),
                            mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                            mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y])
                    );
                    ae = distanceBetweenTwoPoints(
                            mesh.X(mesh.getPointX1()[mesh.getNX() + 1], mesh.getPointX2()[J]),
                            mesh.Y(mesh.getPointX1()[mesh.getNX() + 1], mesh.getPointX2()[J]),
                            mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                            mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y])
                    );
                    aw1 = Math.min(aw, aw1);
                    ae1 = Math.min(ae, ae1);
                }

                temp0 = Math.min(aw1, ae1);

                for (int I = 0; I <= mesh.getNX() + 1; I++) {
                    as = distanceBetweenTwoPoints(
                            mesh.X(mesh.getPointX1()[I], mesh.getPointX2()[0]),
                            mesh.Y(mesh.getPointX1()[I], mesh.getPointX2()[0]),
                            mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                            mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y])
                    );
                    an = distanceBetweenTwoPoints(
                            mesh.X(mesh.getPointX1()[I], mesh.getPointX2()[mesh.getNY() + 1]),
                            mesh.Y(mesh.getPointX1()[I], mesh.getPointX2()[mesh.getNY() + 1]),
                            mesh.X(mesh.getPointX1()[X], mesh.getPointX2()[Y]),
                            mesh.Y(mesh.getPointX1()[X], mesh.getPointX2()[Y])
                    );
                    as1 = Math.min(as, as1);
                    an1 = Math.min(an, an1);
                }

                temp1 = Math.min(as1, an1);
                distance.getFI()[IJ] = Math.min(temp0, temp1);
            }
        }

    }

}
