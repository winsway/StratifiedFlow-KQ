package com.winswe.turbulence;

import static com.cup.turbulence.Turbulence.CAPPA;
import static com.cup.turbulence.Turbulence.ELOG;
import com.winswe.field.VolScalarField;
import com.winswe.io.IOobject;
import com.winswe.matrix.solve.SolverPerformance;
import com.winswe.mesh.Structed2D;
import com.winswe.turbulence.kepsilon.KEpsilonModel;
import com.winswe.turbulence.komega.KOmegaModel;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.List;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午2:59:42
 */
public abstract class Turbulence {

    /**
     * Von Karman Constant
     */
    final public static double CAPPA = 0.41;

    /**
     * A Constant
     */
    final public static double ELOG = 8.342;

    /**
     * Solve Turbulence Model
     */
    public abstract void solve();

    /**
     * velocity
     */
    protected final VolScalarField velocity;

    /**
     * dynamic viscosity
     */
    protected final VolScalarField dynamicViscosity;

    /**
     * density
     */
    protected final VolScalarField density;

    /**
     * mesh
     */
    protected final Structed2D mesh;

    /**
     * io object
     */
    protected final IOobject iOobject;

    /**
     * eddy viscosity
     */
    protected final VolScalarField mut;

    public Turbulence(
            VolScalarField velocity,
            VolScalarField dynamicViscosity,
            VolScalarField density,
            Structed2D mesh,
            IOobject iOobject
    ) {
        this.velocity = velocity;
        this.dynamicViscosity = dynamicViscosity;
        this.density = density;
        this.mesh = mesh;
        this.iOobject = iOobject;

        this.mut = new VolScalarField("Eddy Viscosity", mesh, iOobject);
    }

    /**
     *
     * @param velocity
     * @param dynamicViscosity
     * @param density
     * @param mesh
     * @param iOobject
     * @return turbulence model
     */
    public static Turbulence factory(
            VolScalarField velocity,
            VolScalarField dynamicViscosity,
            VolScalarField density,
            Structed2D mesh,
            IOobject iOobject
    ) {

        String name
                = iOobject.
                        getJsonObject().
                        getJSONObject("turbulence").
                        getString("Model");

        if ("KEpsilon".equals(name)) {
            return new KEpsilonModel(velocity, dynamicViscosity, density, mesh, iOobject);
        } else if ("KOmega".equals(name)) {
            return new KOmegaModel(velocity, dynamicViscosity, density, mesh, iOobject);
        } else if ("SSTKOmega".equals(name)) {
            return null;
        } else {
            throw new ArithmeticException("without the solver name");
        }
    }

    /**
     * sqrt(pow(dx, 2.0) + pow(dy, 2.0));
     *
     * @param x0
     * @param y0
     * @param x1
     * @param y1
     * @return sqrt(pow(dx, 2.0) + pow(dy, 2.0))
     */
    public static double distanceBetweenTwoPoints(
            double x0,
            double y0,
            double x1,
            double y1
    ) {
        double dx, dy;
        dx = x1 - x0;
        dy = y1 - y0;
        return sqrt(pow(dx, 2.0) + pow(dy, 2.0));
    }

    /**
     *
     * @param density density
     * @param mu dynamic viscosity
     * @param dw distance to wall
     * @param tauw strain stress at wall
     * @return dw * sqrt(density * tauw) / mu
     */
    static public double yPlusCore(
            double density,
            double mu,
            double dw,
            double tauw
    ) {
        return dw * sqrt(density * tauw) / mu;
    }

    /**
     *
     * @param density density
     * @param mu dynamic viscosity
     * @param dw distance to wall
     * @param CK A Constant
     * @return density * CK * dw / mu;
     */
    static public double yPlusSubLay(
            double density,
            double mu,
            double dw,
            double CK
    ) {
        return density * CK * dw / mu;
    }

    /**
     *
     * @param yPlus
     * @param mu dynamic viscosity
     * @return yPlus * mu * CAPPA / Math.log(ELOG * yPlus);
     */
    static public double getViscosityWall(double yPlus, double mu) {
        return yPlus * mu * CAPPA / Math.log(ELOG * yPlus);
    }

    /**
     *
     * @param vis viscosity
     * @param uP velocity at near wall
     * @param uWall velocity at wall
     * @param distance distance to wall
     * @return vis * (uP - uWall) / distance;
     */
    static public double getTauWall(
            double vis,
            double uP,
            double uWall,
            double distance
    ) {
        return vis * (uP - uWall) / distance;
    }

    /**
     *
     * @param tau strain stress
     * @param K turbulent kinetic energy
     * @param distance distance to wall
     * @return abs(tau) * Cu^0.25 * sqrt(max(0, K)) / (distance * CAPPA);
     */
    static public double getProduceTerm(
            double tau,
            double K,
            double distance
    ) {
        double CMU25 = sqrt(sqrt(0.09));
        return abs(tau) * CMU25 * sqrt(max(0, K)) / (distance * CAPPA);
    }

    /**
     *
     * @param mesh mesh
     * @param velocity velocity m/s
     * @param S strain rate tensor
     */
    public void calculateS(
            Structed2D mesh,
            VolScalarField velocity,
            VolScalarField S
    ) {
        double alpha;
        double aw, ae, as, an;
        double a, b;

        for (int Y = 1; Y <= mesh.getNY(); ++Y) {
            for (int X = 1; X <= mesh.getNX(); ++X) {
                int IJ = mesh.getCellIndex(X, Y);
                int IJW = mesh.getCellIndex(X - 1, Y);
                int IJE = mesh.getCellIndex(X + 1, Y);
                int IJS = mesh.getCellIndex(X, Y - 1);
                int IJN = mesh.getCellIndex(X, Y + 1);
                alpha = mesh.alpha(mesh.getPointX1()[X], mesh.getPointX2()[Y]);
                aw = (velocity.getFI()[IJ] - velocity.getFI()[IJW])
                        / (mesh.getDXP()[X - 1]);
                ae = (velocity.getFI()[IJE] - velocity.getFI()[IJ])
                        / (mesh.getDXP()[X]);
                as = (velocity.getFI()[IJ] - velocity.getFI()[IJS])
                        / (mesh.getDYP()[Y - 1]);
                an = (velocity.getFI()[IJN] - velocity.getFI()[IJ])
                        / (mesh.getDYP()[Y]);

                a = 0.5 * (aw + ae);

                b = 0.5 * (as + an);

                S.getFI()[IJ] = (a * a + b * b) / alpha;

            }
        }

    }

    public double turbulenceLengthScale(double dh) {
        return 0.07 * dh;
    }

    /**
     * set initial turbulent parameter
     */
    public abstract void setTurPar();

    /**
     * update eddy viscosity
     */
    public abstract void updateMut();

    public VolScalarField getEddyViscosity() {
        return mut;
    }

    public abstract void setSolverPerformance(
            List<SolverPerformance> solverPerformances
    );

}
