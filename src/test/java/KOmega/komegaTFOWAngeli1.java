/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package KOmega;

import com.cup.system.SystemControl;
import com.cup.field.Scalarfield;
import com.cup.field.Field;
import com.cup.boundary.factory.Label;
import com.cup.util.MatrixSolver;
import static com.cup.boundary.Boundary.*;
import com.cup.system.ControlDict;
import com.cup.system.FvSolution;
import com.cup.field.Coefficient;
import com.cup.io.Output;
import com.winswe.io.Writer;
import static com.cup.log.Residual.getRes;
import com.cup.mesh.Mesh;
import com.cup.mesh.factory.BiSameV1;
import com.cup.mesh.factory.BiSingleV1;
import com.cup.turbulence.Turbulence;
import static com.cup.util.Flux.faceValue;
import com.cup.util.Tool;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.exp;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import java.util.Objects;

/**
 *
 * @see de Sampaio, P.A.B., J.L.H. Faccini and J. Su, Modelling of stratified
 * gas–liquid two-phase flow in horizontal circular pipes. International Journal
 * of Heat and Mass Transfer, 2008. 51(11-12): p. 2752-2761. 液-液
 * <n>1.引入界面摩擦因素</n>
 * <n>2.解决压力对不上的问题</n>
 * <n>3.解决速度场对不上的问题</n>
 * @author 齐雪宇
 */
public final class komegaTFOWAngeli1 {

    double[][][] S, rho, mu, mut, phi;
    double[][][] Soil, rhoOil, muOil, mutOil, lambdaOil, cpOil;
    double[][][] Swater, rhoWater, muWater, mutWater, lambdaWater, cpWater;
    // 其他参量
    FvSolution fvsolution;
    ControlDict controldict;
    SystemControl sys, sysOil, sysWater;
    Mesh mesh;
    BiSameV1 bm;
    BiSingleV1 OilMesh;
    BiSingleV1 WaterMesh;
    //场变量
    Field W, K, Omega, T;
    Field WOil, KOil, OmegaOil, TOil;
    Field WWater, KWater, OmegaWater, TWater;
    //系数
    Coefficient coeW, coeK, coeO, coeT;
    Coefficient coeWoil, coeKoil, coeOoil, coeToil;
    Coefficient coeWwater, coeKwater, coeOwater, coeTwater;
    double dpdz, hl;
    double Qso, Qsw;
    double R;
    double g = 0;

    public boolean createFields() {
        int x, y, z;
        x = mesh.numPx();
        y = mesh.numPy();
        z = mesh.numPz();
        //建立场
        W = new Scalarfield(x, y, z);
        K = new Scalarfield(x, y, z);
        Omega = new Scalarfield(x, y, z);
        T = new Scalarfield(x, y, z);
        S = new double[x][y][z];
        rho = new double[x][y][z];
        mu = new double[x][y][z];
        mut = new double[x][y][z];
        phi = new double[x][y][z];
        //方程系数
        coeW = new Coefficient(x, y);
        coeK = new Coefficient(x, y);
        coeO = new Coefficient(x, y);
        coeT = new Coefficient(x, y);
        //   油相建立场
        x = OilMesh.numPx();
        y = OilMesh.numPy();
        z = OilMesh.numPz();
        WOil = new Scalarfield(x, y, z);
        TOil = new Scalarfield(x, y, z);
        KOil = new Scalarfield(x, y, z);
        OmegaOil = new Scalarfield(x, y, z);
        Soil = new double[x][y][z];
        rhoOil = new double[x][y][z];
        muOil = new double[x][y][z];
        mutOil = new double[x][y][z];
        lambdaOil = new double[x][y][z];
        cpOil = new double[x][y][z];
        coeWoil = new Coefficient(x, y);
        coeKoil = new Coefficient(x, y);
        coeOoil = new Coefficient(x, y);
        coeToil = new Coefficient(x, y);
        //    水相建立场
        x = WaterMesh.numPx();
        y = WaterMesh.numPy();
        z = WaterMesh.numPz();
        WWater = new Scalarfield(x, y, z);
        TWater = new Scalarfield(x, y, z);
        KWater = new Scalarfield(x, y, z);
        OmegaWater = new Scalarfield(x, y, z);
        Swater = new double[x][y][z];
        rhoWater = new double[x][y][z];
        muWater = new double[x][y][z];
        mutWater = new double[x][y][z];
        lambdaWater = new double[x][y][z];
        cpWater = new double[x][y][z];

        coeWwater = new Coefficient(x, y);
        coeKwater = new Coefficient(x, y);
        coeOwater = new Coefficient(x, y);
        coeTwater = new Coefficient(x, y);
        return true;
    }

    public boolean initField() {
//<editor-fold>
        int k = 1;
        for (int i = 0; i < mesh.numPx(); i++) {
            for (int j = 0; j < mesh.numPy(); j++) {
                if (i == 0 || i == (mesh.numPx() - 1)
                        || j == 0 || j == (mesh.numPy() - 1)) {
                    W.getBoundary().setType(i, j, k, 1);
                    K.getBoundary().setType(i, j, k, 1);
                    Omega.getBoundary().setType(i, j, k, 1);
                }
                if (j <= bm.interhL) {
                    phi[i][j][k] = -1;
                    rho[i][j][k] = 1000;
                    mu[i][j][k] = 1.0e-3;
                } else {
                    phi[i][j][k] = 1;
                    rho[i][j][k] = 1000;
                    mu[i][j][k] = 1.0e-3;
                }
            }
        }
        //</editor-fold>
//        油相设定初场
//<editor-fold>
        k = 1;
        for (int i = 0; i < OilMesh.numPx(); i++) {
            for (int j = 0; j < OilMesh.numPy(); j++) {

                if (i == 0 || i == (OilMesh.numPx() - 1)
                        || j == 0 || j == (OilMesh.numPy() - 1)) {
                    WOil.getBoundary().setType(i, j, k, 1);
                    KOil.getBoundary().setType(i, j, k, 1);
                    OmegaOil.getBoundary().setType(i, j, k, 1);
                    TOil.getBoundary().setType(i, j, k, 1);
                }
                rhoOil[i][j][k] = 801;
                muOil[i][j][k] = 1.6e-3;
                lambdaOil[i][j][k] = 0.1260;
                cpOil[i][j][k] = 2000;
            }
        }
        //</editor-fold>
//        水相设定初场
//<editor-fold>
        k = 1;
        for (int i = 0; i < WaterMesh.numPx(); i++) {
            for (int j = 0; j < WaterMesh.numPy(); j++) {

                if (i == 0 || i == (WaterMesh.numPx() - 1)
                        || j == 0 || j == (WaterMesh.numPy() - 1)) {
                    WWater.getBoundary().setType(i, j, k, 1);
                    KWater.getBoundary().setType(i, j, k, 1);
                    OmegaWater.getBoundary().setType(i, j, k, 1);
                    TWater.getBoundary().setType(i, j, k, 1);
                }
                rhoWater[i][j][k] = 996;
                muWater[i][j][k] = 8.6e-4;
                lambdaWater[i][j][k] = 0.6;
                cpWater[i][j][k] = 4186;
            }
        }
        //</editor-fold>
        return true;
    }

    public boolean blockMesh() {
        bm = new BiSameV1(12, Math.PI, 100, 100);
        bm.radia = R;
        double D = bm.radia * 2.0;
        hl = D * hl;
        double gamma = acos(1.0 - 2 * hl / D);
        bm.theta1 = gamma;
        bm.theta2 = PI;
        mesh = bm;
        mesh.blockMesh();
//       设定两相网格
//       油相
        System.out.println("begin OilMesh");
        OilMesh = new BiSingleV1(12, bm.theta2 - bm.theta1,
                bm.getNX(), bm.interhL);
        OilMesh.theta1 = bm.theta1;
        OilMesh.theta2 = bm.theta2;
        OilMesh.a = bm.a;
        OilMesh.radia = bm.radia;
        OilMesh.blockMesh();
//       水相
        System.out.println("begin LiquidMesh");
        WaterMesh = new BiSingleV1(12, bm.theta3 - bm.theta2,
                bm.getNX(), bm.getNY() - bm.interhL);
        WaterMesh.theta1 = bm.theta2;
        WaterMesh.theta2 = bm.theta3;
        WaterMesh.a = bm.a;
        WaterMesh.radia = bm.radia;
        WaterMesh.blockMesh();
        return true;
    }

    public boolean Solver() {
        scheme5();
        return true;
    }

    /**
     * 单圆管层流条件下的问题
     */
    void scheme1(double dpdz) {
        double res = 0;
        int count = 0;
        Boolean judge = Boolean.FALSE;
//        进入大循环
        do {
            count++;
            System.out.println("计算速度场W");
            cacW(dpdz, mesh, W, mu, rho, mut, sys);
            res = getRes(W.getNewField(), W.getOldField());
            W.copyNew2Old();
            System.out.println("res=" + res);
            outPut(count, judge, "W", W.getNewField(), mesh);
            System.out.println("计算总流量=" + Tool.getTotalVol(W, mesh));
        } while (res > 1e-5);
        judge = Boolean.TRUE;
        outPut(count, judge, "W", W.getNewField(), mesh);
    }

    /**
     * 单圆管K-Omega 引入湍流粘度
     */
    void scheme2(double dpdz) {
        double res = 0, res1;
        int count = 0;
        Boolean judge = Boolean.TRUE;
//假定速度场
        cacW(dpdz, mesh, W, mu, rho, mut, sys);
        W.copyNew2Old();
        outPut(count, judge, "W", W.getNewField(), mesh);
        System.out.println("初始化湍流动能！");
        setK(K, W);
        K.copyNew2Old();
        outPut(count, judge, "K", K.getNewField(), mesh);
        System.out.println("初始化湍流Omega！");
        setOmega();
        Omega.copyNew2Old();
        outPut(count, judge, "Omega", Omega.getNewField(), mesh);
//进入循环
        do {
            judge = Boolean.FALSE;
            count++;
            Turbulence.cacScenter(W, mesh, S);
            outPut(count, judge, "S", S, mesh);
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMut(rho, mut, K, Omega);
            outPut(count, judge, "Mut", mut, mesh);
            //</editor-fold>
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacK();
            K.copyNew2Old();
            outPut(count, judge, "K", K.getNewField(), mesh);
            //</editor-fold>
            //<editor-fold>能量耗散率
            System.out.println("计算Omega！");
            cacOmega();
            outPut(count, judge, "Omega", Omega.getNewField(), mesh);
            res1 = getRes(Omega.getNewField(), Omega.getOldField());
            System.out.println("湍流残差res1=" + res1);
            Omega.copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, mesh, W, mu, rho, mut, sys);
            res = getRes(W.getNewField(), W.getOldField());
            W.copyNew2Old();
            System.out.println("resW=" + res);
            System.out.println("计算总流量=" + Tool.getTotalVol(W, mesh));
            outPut(count, judge, "W", W.getNewField(), mesh);
            //</editor-fold>
        } while (res > 1e-4);
        judge = Boolean.TRUE;
        outPut(count, judge, "W", W.getNewField(), mesh);
    }

    /**
     * 单相双流体方程求解
     */
    void scheme3() {
        double resw = 0, reso = 0;
        int count = 0;
        double Qw, Qo;
        Boolean judge = Boolean.FALSE;
//        进入大循环
        do {
            count++;
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
            cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
            updateInterfaceW();
            resw = getRes(WWater.getNewField(), WWater.getOldField());
            reso = getRes(WOil.getNewField(), WOil.getOldField());
            WWater.copyNew2Old();
            WOil.copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
            outPut(count, judge, "WWater", WWater.getNewField(), WaterMesh);
            Qw = Tool.getTotalVol(WWater, WaterMesh);
            Qo = Tool.getTotalVol(WOil, OilMesh);
            System.out.println("计算总流量=Qw" + Qw);
            System.out.println("计算总流量=Qo" + Qo);
            System.out.println("总流量=" + (Qw + Qo));
        } while (max(resw, reso) > 1e-3 || count < 500);
        judge = Boolean.TRUE;
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
    }

    /**
     * 双流体方程+湍流模型求解
     */
    void scheme4() {
        double resw = 0, reso = 0;
        double reswt = 0, resot = 0;
        double Qw, Qo;
        int count = 0;
        Boolean judge = Boolean.TRUE;
//      假定速度场
        cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
        cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
        WWater.copyNew2Old();
        WOil.copyNew2Old();
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流动能！");
        setK(KOil, WOil);
        KOil.copyNew2Old();
        setK(KWater, WWater);
        KWater.copyNew2Old();
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流能量耗散率！");
        setOmega(OmegaOil, KOil, OilMesh);
        OmegaOil.copyNew2Old();
        setOmega(OmegaWater, KWater, WaterMesh);
        OmegaWater.copyNew2Old();
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
//        进入大循环
        do {
            judge = Boolean.FALSE;
            count++;
            Turbulence.cacScenter(WOil, OilMesh, Soil);
            Turbulence.cacScenter(WWater, WaterMesh, Swater);
            outPut(count, judge, "SOil", Soil, OilMesh);
            outPut(count, judge, "Swater", Swater, WaterMesh);
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMut(rhoOil, mutOil, KOil, OmegaOil);
            cacMut(rhoWater, mutWater, KWater, OmegaWater);
            outPut(count, judge, "MutOil", mutOil, OilMesh);
            outPut(count, judge, "Mutwater", mutWater, WaterMesh);
            //</editor-fold>
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacKOil();
            cacKWater();
            KOil.copyNew2Old();
            KWater.copyNew2Old();
            outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
            outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
            //</editor-fold>
            //<editor-fold>湍流Omega
            System.out.println("计算Omega！");
            cacOmegaOil();
            cacOmegaWater();
            outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
            outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
            resot = getRes(OmegaOil.getNewField(), OmegaOil.getOldField());
            reswt = getRes(OmegaWater.getNewField(), OmegaWater.getOldField());
            System.out.println("湍流残差resot = " + resot + "湍流残差reswt = " + reswt);
            OmegaOil.copyNew2Old();
            OmegaWater.copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
            cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
            updateInterfaceW();
            resw = getRes(WWater.getNewField(), WWater.getOldField());
            reso = getRes(WOil.getNewField(), WOil.getOldField());
            WWater.copyNew2Old();
            WOil.copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
            outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
//            Qw = Mesh.getTotalVol(WWater, WaterMesh);
//            Qo = Mesh.getTotalVol(WOil, OilMesh);
//            System.out.println("计算总流量=Qw" + Qw);
//            System.out.println("计算总流量=Qo" + Qo);
//            System.out.println("总流量=" + (Qw + Qo));
            //</editor-fold>
        } while (max(resw, reso) > 1e-3 || count < 500);
        judge = Boolean.TRUE;
        outPut(count, judge, "SOil", Soil, OilMesh);
        outPut(count, judge, "Swater", Swater, WaterMesh);
        outPut(count, judge, "MutOil", mutOil, OilMesh);
        outPut(count, judge, "Mutwater", mutWater, WaterMesh);
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
    }

    /**
     * 双流体方程+湍流模型求解+界面摩擦
     */
    void scheme5() {
        double resw = 0, reso = 0;
        double reswt = 0, resot = 0;
        int count = 0;
        Boolean judge = Boolean.TRUE;
//      假定速度场
        cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
        cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
        WWater.copyNew2Old();
        WOil.copyNew2Old();
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流动能！");
        setK(KOil, WOil);
        KOil.copyNew2Old();
        setK(KWater, WWater);
        KWater.copyNew2Old();
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流能量耗散率！");
        setOmega(OmegaOil, KOil, OilMesh);
        OmegaOil.copyNew2Old();
        setOmega(OmegaWater, KWater, WaterMesh);
        OmegaWater.copyNew2Old();
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
//        进入大循环
        do {
            judge = Boolean.FALSE;
            count++;
            Turbulence.cacScenter1(WOil, OilMesh, Soil);
            Turbulence.cacScenter1(WWater, WaterMesh, Swater);
            outPut(count, judge, "SOil", Soil, OilMesh);
            outPut(count, judge, "Swater", Swater, WaterMesh);
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMut(rhoOil, mutOil, KOil, OmegaOil);
            cacMut(rhoWater, mutWater, KWater, OmegaWater);
            outPut(count, judge, "MutOil", mutOil, OilMesh);
            outPut(count, judge, "Mutwater", mutWater, WaterMesh);
            //</editor-fold>
            cacVelocity();
            cacTauInter();
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacKOil();
            cacKWater();
            KOil.copyNew2Old();
            KWater.copyNew2Old();
            outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
            outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
            //</editor-fold>
            //<editor-fold>湍流Omega
            System.out.println("计算Omega！");
            cacOmegaOil();
            cacOmegaWater();
            outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
            outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
            resot = getRes(OmegaOil.getNewField(), OmegaOil.getOldField());
            reswt = getRes(OmegaWater.getNewField(), OmegaWater.getOldField());
            System.out.println("湍流残差resot = " + resot + " 湍流残差reswt = " + reswt);
            OmegaOil.copyNew2Old();
            OmegaWater.copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
            cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
            updateInterfaceW();
            resw = getRes(WWater.getNewField(), WWater.getOldField());
            reso = getRes(WOil.getNewField(), WOil.getOldField());
            WWater.copyNew2Old();
            WOil.copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
            outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
            //</editor-fold>
        } while (max(resw, reso) > 1e-3 && count < 10000);
        judge = Boolean.TRUE;
        outPut(count, judge, "SOil", Soil, OilMesh);
        outPut(count, judge, "Swater", Swater, WaterMesh);
        outPut(count, judge, "MutOil", mutOil, OilMesh);
        outPut(count, judge, "Mutwater", mutWater, WaterMesh);
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
    }

    /**
     * 双流体方程+湍流模型求解+界面摩擦+非牛顿
     */
    void scheme6() {
        double resw = 0, reso = 0;
        double reswt = 0, resot = 0;
        int count = 0;
        Boolean judge = Boolean.TRUE;
//      假定速度场
        cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
        cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
        WWater.copyNew2Old();
        WOil.copyNew2Old();
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流动能！");
        setK(KOil, WOil);
        KOil.copyNew2Old();
        setK(KWater, WWater);
        KWater.copyNew2Old();
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        System.out.println("初始化湍流能量耗散率！");
        setOmega(OmegaOil, KOil, OilMesh);
        OmegaOil.copyNew2Old();
        setOmega(OmegaWater, KWater, WaterMesh);
        OmegaWater.copyNew2Old();
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
//        进入大循环
        do {
            judge = Boolean.FALSE;
            count++;
            Turbulence.cacScenter(WOil, OilMesh, Soil);
            Turbulence.cacScenter(WWater, WaterMesh, Swater);
            outPut(count, judge, "SOil", Soil, OilMesh);
            outPut(count, judge, "Swater", Swater, WaterMesh);
            //<editor-fold>湍流粘度
            System.out.println("计算湍流粘度！");
            cacMu();
            cacMut(rhoOil, mutOil, KOil, OmegaOil);
            cacMut(rhoWater, mutWater, KWater, OmegaWater);
            outPut(count, judge, "MutOil", mutOil, OilMesh);
            outPut(count, judge, "Mutwater", mutWater, WaterMesh);
            //</editor-fold>
            cacTauInter();
            //<editor-fold>湍流动能
            System.out.println("计算湍流动能！");
            cacKOil();
            cacKWater();
            KOil.copyNew2Old();
            KWater.copyNew2Old();
            outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
            outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
            //</editor-fold>
            //<editor-fold>湍流Omega
            System.out.println("计算Omega！");
            cacOmegaOil();
            cacOmegaWater();
            outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
            outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
            resot = getRes(OmegaOil.getNewField(), OmegaOil.getOldField());
            reswt = getRes(OmegaWater.getNewField(), OmegaWater.getOldField());
            System.out.println("湍流残差resot = " + resot + "湍流残差reswt = " + reswt);
            OmegaOil.copyNew2Old();
            OmegaWater.copyNew2Old();
            //</editor-fold>
            //<editor-fold>计算速度
            System.out.println("计算速度场W");
            cacW(dpdz, WaterMesh, WWater, muWater, rhoWater, mutWater, sysWater);
            cacW(dpdz, OilMesh, WOil, muOil, rhoOil, mutOil, sysOil);
            updateInterfaceW();
            resw = getRes(WWater.getNewField(), WWater.getOldField());
            reso = getRes(WOil.getNewField(), WOil.getOldField());
            WWater.copyNew2Old();
            WOil.copyNew2Old();
            System.out.println("resw=" + resw + " reso=" + reso);
            outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
            outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
            //</editor-fold>
        } while (max(resw, reso) > 1e-3 || count < 500);
        judge = Boolean.TRUE;
        outPut(count, judge, "SOil", Soil, OilMesh);
        outPut(count, judge, "Swater", Swater, WaterMesh);
        outPut(count, judge, "MutOil", mutOil, OilMesh);
        outPut(count, judge, "Mutwater", mutWater, WaterMesh);
        outPut(count, judge, "KOil", KOil.getNewField(), OilMesh);
        outPut(count, judge, "KWater", KWater.getNewField(), WaterMesh);
        outPut(count, judge, "OmegaOil", OmegaOil.getNewField(), OilMesh);
        outPut(count, judge, "OmegaWater", OmegaWater.getNewField(), WaterMesh);
        outPut(count, judge, "WOil", WOil.getNewField(), OilMesh);
        outPut(count, judge, "Wwater", WWater.getNewField(), WaterMesh);
    }

    /**
     * @param Qw Qw
     * @return
     */
    double F(double Qw) {
        double Qcwater = 0;
        Qcwater = Tool.getTotalVol(WWater, WaterMesh);
        System.out.println("Qcwate = " + Qcwater);
        return Qcwater - Qw;
    }

    /**
     *
     * @param Qoil Qoil
     * @return
     */
    double G(double Qoil) {
        double Qcoil = 0;
        Qcoil = Tool.getTotalVol(WOil, OilMesh);
        System.out.println("Qcoil = " + Qcoil);
        return Qcoil - Qoil;
    }

    /**
     * 计算湍流动能场
     */
    void cacK() {
        double URF = 0.6;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simga2 = 0.5;
        double betastar = 0.09;
        Label flagK = new Label("K");
        for (int z = 1; z < K.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < K.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < K.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagK.setFlag(K.getBoundary(), x, y);
                    double vol
                            = mesh.getDXU()[x - 1] * mesh.getDYV()[y - 1]
                            * mesh.J(mesh.gettPx()[x], mesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (mu[x][y][z] + mut[x][y][z] * simga2);
                    muw = (mu[x - 1][y][z] + mut[x - 1][y][z] * simga2);
                    mue = (mu[x + 1][y][z] + mut[x + 1][y][z] * simga2);
                    mus = (mu[x][y - 1][z] + mut[x][y - 1][z] * simga2);
                    mun = (mu[x][y + 1][z] + mut[x][y + 1][z] * simga2);
//得到gamma
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[x - 1], mesh.getTx()[x - 1], mesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[x], mesh.getTx()[x], mesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[y - 1], mesh.getTy()[y - 1], mesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[y], mesh.getTy()[y], mesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * mesh.getDYV()[y - 1] / mesh.getDXP()[x - 1];
                    De = game * mesh.getDYV()[y - 1] / mesh.getDXP()[x];
                    Ds = gams * mesh.getDXU()[x - 1] / mesh.getDYP()[y - 1];
                    Dn = gamn * mesh.getDXU()[x - 1] / mesh.getDYP()[y];
//系数
                    coeK.aw[x][y][z] = (1 - flagK.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeK.ae[x][y][z] = (1 - flagK.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeK.as[x][y][z] = (1 - flagK.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeK.an[x][y][z] = (1 - flagK.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -betastar * rho[x][y][z] * Omega.getNewField()[x][y][z];
                    Spad
                            = flagK.w
                            * Spad1(K.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagK.e
                            * Spad2(K.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagK.s
                            * Spad1(K.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagK.n
                            * Spad2(K.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, 1);
//                    得到ap系数
                    coeK.ap[x][y][z]
                            = coeK.aw[x][y][z] + coeK.ae[x][y][z]
                            + coeK.as[x][y][z] + coeK.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;

//                    得到源项
                    Sc = (mut[x][y][z] * S[x][y][z]);
                    Scad
                            = flagK.w
                            * Scad1(K.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, K.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagK.e
                            * Scad2(K.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, K.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagK.s
                            * Scad1(K.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, K.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagK.n
                            * Scad2(K.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, K.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeK.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
//Under-relaxizaion
                    coeK.b[x][y][z]
                            = coeK.b[x][y][z] + (URFFI - 1.0) * coeK.ap[x][y][z]
                            * K.getOldField()[x][y][z];
                    coeK.ap[x][y][z] = (coeK.ap[x][y][z]) * URFFI;
                }
            }
        }

        K.getLog().rl = MatrixSolver.gaussSeidel(K.getNewField(), coeK, sys);
        K.getLog().log("K");
    }

    /**
     * 计算油相湍流动能场
     */
    void cacKOil() {
        double URF = 0.8;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simga2 = 0.5;
        double betastar = 0.09;
        double GEN = 0;
        double CTRANS = 11.63;
        double CAPPA = 0.41;
        Label flagK = new Label("KOil");
        for (int z = 1; z < KOil.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < KOil.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < KOil.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagK.setFlag(KOil.getBoundary(), x, y);
                    double vol
                            = OilMesh.getDXU()[x - 1] * OilMesh.getDYV()[y - 1]
                            * OilMesh.J(OilMesh.gettPx()[x], OilMesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (muOil[x][y][z] + mutOil[x][y][z] * simga2);
                    muw = (muOil[x - 1][y][z] + mutOil[x - 1][y][z] * simga2);
                    mue = (muOil[x + 1][y][z] + mutOil[x + 1][y][z] * simga2);
                    mus = (muOil[x][y - 1][z] + mutOil[x][y - 1][z] * simga2);
                    mun = (muOil[x][y + 1][z] + mutOil[x][y + 1][z] * simga2);
//得到gamma
                    gamw = faceValue(muw, mup - muw, OilMesh.gettPx()[x - 1], OilMesh.getTx()[x - 1], OilMesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, OilMesh.gettPx()[x], OilMesh.getTx()[x], OilMesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, OilMesh.gettPy()[y - 1], OilMesh.getTy()[y - 1], OilMesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, OilMesh.gettPy()[y], OilMesh.getTy()[y], OilMesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * OilMesh.getDYV()[y - 1] / OilMesh.getDXP()[x - 1];
                    De = game * OilMesh.getDYV()[y - 1] / OilMesh.getDXP()[x];
                    Ds = gams * OilMesh.getDXU()[x - 1] / OilMesh.getDYP()[y - 1];
                    Dn = gamn * OilMesh.getDXU()[x - 1] / OilMesh.getDYP()[y];
//系数
                    coeKoil.aw[x][y][z] = (1 - flagK.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeKoil.ae[x][y][z] = (1 - flagK.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeKoil.as[x][y][z] = (1 - flagK.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeKoil.an[x][y][z] = (1 - flagK.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -betastar * rhoOil[x][y][z] * OmegaOil.getNewField()[x][y][z];
                    Spad
                            = flagK.w
                            * Spad1(KOil.getBoundary().getType(x - 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagK.e
                            * Spad2(KOil.getBoundary().getType(x + 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagK.s
                            * Spad1(KOil.getBoundary().getType(x, y - 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagK.n
                            * Spad2(KOil.getBoundary().getType(x, y + 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y],
                                    gamn, Fn, 1);
//                    得到ap系数
                    coeKoil.ap[x][y][z]
                            = coeKoil.aw[x][y][z] + coeKoil.ae[x][y][z]
                            + coeKoil.as[x][y][z] + coeKoil.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;

//                    得到源项
                    //<editor-fold>
                    if ((flagK.w == 1) || (flagK.e == 1)
                            || (flagK.s == 1) || (flagK.n == 1)) {
                        double VISS = muOil[x][y][z];
                        double YPL, DN = 0;
                        double TAU = 0;
                        double CMU = 0.09;
                        double CMU25 = sqrt(sqrt(CMU));
                        double ELOG = 8.342;
                        double CK = CMU25 * sqrt(max(0, KOil.getOldField()[x][y][z]));
                        double VISCW, VISW;
                        double YPL1, epl, es, Uint;
                        if (flagK.w == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x - 1], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x - 1], OilMesh.gettPy()[y]));
                            YPL = rhoOil[x][y][z] * CK * DN / muOil[x][y][z];
                            VISCW = YPL * muOil[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muOil[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WOil.getNewField()[x][y][z] - WOil.getNewField()[x - 1][y][z]) / DN;
                        }
                        if (flagK.e == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x + 1], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x + 1], OilMesh.gettPy()[y]));
                            YPL = rhoOil[x][y][z] * CK * DN / muOil[x][y][z];
                            VISCW = YPL * muOil[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muOil[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WOil.getNewField()[x][y][z] - WOil.getNewField()[x + 1][y][z]) / DN;
                        }
                        if (flagK.s == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y - 1]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y - 1]));
                            YPL = rhoOil[x][y][z] * CK * DN / muOil[x][y][z];
                            VISCW = YPL * muOil[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muOil[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WOil.getNewField()[x][y][z] - WOil.getNewField()[x][y - 1][z]) / DN;
                        }
                        if (flagK.n == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y + 1]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y + 1]));
                            YPL = rhoOil[x][y][z] * CK * DN / muOil[x][y][z];
//
                            Uint = sqrt(abs(Tauo) / rhoOil[x][y][z]);
                            es = cacFintES(Tauo, rhoOil[x][y][z], x, y, z);
                            epl = rhoOil[x][y][z] * es * Uint / muOil[x][y][z];
                            if (epl <= 4.535) {
                                YPL1 = 0;
                            } else {
                                YPL1 = max(0, 0.9 * sqrt(epl) - epl * exp(-epl / 6.0));
                            }
                            YPL = YPL + YPL1;
//
                            VISCW = YPL * muOil[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muOil[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WOil.getNewField()[x][y][z] - WOil.getNewField()[x][y + 1][z]) / DN;
                        }
                        GEN
                                = abs(TAU) * CMU25 * sqrt(max(0, KOil.getOldField()[x][y][z]))
                                / (DN * CAPPA);
                    } else {
                        GEN = (mutOil[x][y][z] * Soil[x][y][z]);
                    }
                    //</editor-fold>

                    Sc = GEN;
                    Scad
                            = flagK.w
                            * Scad1(KOil.getBoundary().getType(x - 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x - 1],
                                    gamw, Fw, KOil.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagK.e
                            * Scad2(KOil.getBoundary().getType(x + 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x],
                                    game, Fe, KOil.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagK.s
                            * Scad1(KOil.getBoundary().getType(x, y - 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y - 1],
                                    gams, Fs, KOil.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagK.n
                            * Scad2(KOil.getBoundary().getType(x, y + 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y],
                                    gamn, Fn, KOil.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeKoil.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
//Under-relaxizaion
                    coeKoil.b[x][y][z]
                            = coeKoil.b[x][y][z] + (URFFI - 1.0) * coeKoil.ap[x][y][z]
                            * KOil.getOldField()[x][y][z];
                    coeKoil.ap[x][y][z] = (coeKoil.ap[x][y][z]) * URFFI;
                }
            }
        }
        KOil.getLog().rl = MatrixSolver.gaussSeidel(KOil.getNewField(), coeKoil, sysOil);
        KOil.getLog().log("Koil");
    }

    /**
     * 计算水相湍流动能场
     */
    void cacKWater() {
        double URF = 0.8;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double simga2 = 0.5;
        double betastar = 0.09;
        double GEN = 0;
        double CTRANS = 11.63;
        double CAPPA = 0.41;
        Label flagK = new Label("KOil");
        for (int z = 1; z < KWater.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < KWater.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < KWater.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagK.setFlag(KWater.getBoundary(), x, y);
                    double vol
                            = WaterMesh.getDXU()[x - 1] * WaterMesh.getDYV()[y - 1]
                            * WaterMesh.J(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (muWater[x][y][z] + mutWater[x][y][z] * simga2);
                    muw = (muWater[x - 1][y][z] + mutWater[x - 1][y][z] * simga2);
                    mue = (muWater[x + 1][y][z] + mutWater[x + 1][y][z] * simga2);
                    mus = (muWater[x][y - 1][z] + mutWater[x][y - 1][z] * simga2);
                    mun = (muWater[x][y + 1][z] + mutWater[x][y + 1][z] * simga2);
//得到gamma
                    gamw = faceValue(muw, mup - muw, WaterMesh.gettPx()[x - 1], WaterMesh.getTx()[x - 1], WaterMesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, WaterMesh.gettPx()[x], WaterMesh.getTx()[x], WaterMesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, WaterMesh.gettPy()[y - 1], WaterMesh.getTy()[y - 1], WaterMesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, WaterMesh.gettPy()[y], WaterMesh.getTy()[y], WaterMesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * WaterMesh.getDYV()[y - 1] / WaterMesh.getDXP()[x - 1];
                    De = game * WaterMesh.getDYV()[y - 1] / WaterMesh.getDXP()[x];
                    Ds = gams * WaterMesh.getDXU()[x - 1] / WaterMesh.getDYP()[y - 1];
                    Dn = gamn * WaterMesh.getDXU()[x - 1] / WaterMesh.getDYP()[y];
//系数
                    coeKwater.aw[x][y][z] = (1 - flagK.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeKwater.ae[x][y][z] = (1 - flagK.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeKwater.as[x][y][z] = (1 - flagK.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeKwater.an[x][y][z] = (1 - flagK.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -betastar * rhoWater[x][y][z] * OmegaWater.getNewField()[x][y][z];
                    Spad
                            = flagK.w
                            * Spad1(KWater.getBoundary().getType(x - 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagK.e
                            * Spad2(KWater.getBoundary().getType(x + 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagK.s
                            * Spad1(KWater.getBoundary().getType(x, y - 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagK.n
                            * Spad2(KWater.getBoundary().getType(x, y + 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y],
                                    gamn, Fn, 1);
//                    得到ap系数
                    coeKwater.ap[x][y][z]
                            = coeKwater.aw[x][y][z] + coeKwater.ae[x][y][z]
                            + coeKwater.as[x][y][z] + coeKwater.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//                    得到源项
//<editor-fold>
                    if ((flagK.w == 1) || (flagK.e == 1)
                            || (flagK.s == 1) || (flagK.n == 1)) {
                        double VISS = muWater[x][y][z];
                        double YPL, DN = 0;
                        double TAU = 0;
                        double CMU = 0.09;
                        double CMU25 = sqrt(sqrt(CMU));
                        double ELOG = 8.342;
                        double CK = CMU25 * sqrt(max(0, KWater.getOldField()[x][y][z]));
                        double VISCW, VISW;
                        double YPL1, epl, es, Uint;
                        if (flagK.w == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x - 1], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x - 1], WaterMesh.gettPy()[y]));
                            YPL = rhoWater[x][y][z] * CK * DN / muWater[x][y][z];
                            VISCW = YPL * muWater[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muWater[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WWater.getNewField()[x][y][z] - WWater.getNewField()[x - 1][y][z]) / DN;
                        }
                        if (flagK.e == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x + 1], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x + 1], WaterMesh.gettPy()[y]));
                            YPL = rhoWater[x][y][z] * CK * DN / muWater[x][y][z];
                            VISCW = YPL * muWater[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muWater[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WWater.getNewField()[x][y][z] - WWater.getNewField()[x + 1][y][z]) / DN;
                        }
                        if (flagK.s == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y - 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y - 1]));
                            YPL = rhoWater[x][y][z] * CK * DN / muWater[x][y][z];
//                            
                            Uint = sqrt(abs(Tauw) / rhoWater[x][y][z]);
                            es = cacFintES(Tauw, rhoWater[x][y][z], x, y, z);
                            epl = rhoWater[x][y][z] * es * Uint / muWater[x][y][z];
                            if (epl <= 4.535) {
                                YPL1 = 0;
                            } else {
                                YPL1 = max(0, 0.9 * sqrt(epl) - epl * exp(-epl / 6.0));
                            }
                            YPL = YPL + YPL1;
//
                            VISCW = YPL * muWater[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muWater[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WWater.getNewField()[x][y][z] - WWater.getNewField()[x][y - 1][z]) / DN;

                        }
                        if (flagK.n == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y + 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y + 1]));
                            YPL = rhoWater[x][y][z] * CK * DN / muWater[x][y][z];
                            VISCW = YPL * muWater[x][y][z] * CAPPA / Math.log(ELOG * YPL);
                            VISW = max(muWater[x][y][z], VISCW);
                            if (YPL >= CTRANS) {
                                VISS = VISW;
                            }
                            TAU = VISS * (WWater.getNewField()[x][y][z] - WWater.getNewField()[x][y + 1][z]) / DN;
                        }
                        GEN
                                = abs(TAU) * CMU25 * sqrt(max(0, KWater.getOldField()[x][y][z]))
                                / (DN * CAPPA);
                    } else {
                        GEN = (mutWater[x][y][z] * Swater[x][y][z]);
                    }
//</editor-fold>
                    Sc = GEN;
                    Scad
                            = flagK.w
                            * Scad1(KWater.getBoundary().getType(x - 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x - 1],
                                    gamw, Fw, KWater.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagK.e
                            * Scad2(KWater.getBoundary().getType(x + 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x],
                                    game, Fe, KWater.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagK.s
                            * Scad1(KWater.getBoundary().getType(x, y - 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y - 1],
                                    gams, Fs, KWater.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagK.n
                            * Scad2(KWater.getBoundary().getType(x, y + 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y],
                                    gamn, Fn, KWater.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeKwater.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
//Under-relaxizaion
                    coeKwater.b[x][y][z]
                            = coeKwater.b[x][y][z] + (URFFI - 1.0) * coeKwater.ap[x][y][z]
                            * KWater.getOldField()[x][y][z];
                    coeKwater.ap[x][y][z] = (coeKwater.ap[x][y][z]) * URFFI;
                }
            }
        }
        KWater.getLog().rl = MatrixSolver.gaussSeidel(KWater.getNewField(), coeKwater, sysWater);
        KWater.getLog().log("Kwater");
    }

    double Tauo = 0, Tauw = 0;
    double eso = 0, esw = 0;
    double wso, wsw;

    /**
     * 计算表观速度
     */
    void cacVelocity() {
        wso = Tool.getTotalVol(WOil, OilMesh) / (PI * this.R * this.R);
        wsw = Tool.getTotalVol(WWater, WaterMesh) / (PI * this.R * this.R);
    }

    /**
     * 求解界面的摩擦
     */
    void cacTauInter() {
        double yp, leta;
        Tauo = 0;
        Tauw = 0;
//      K油相
        for (int x = 1; x < KOil.getNewField().length - 1; ++x) {
            leta = sqrt(OilMesh.alpha(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]));
            yp = Tool.distance(
                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]),
                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]),
                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]));
            Tauo
                    = Tauo + 1.0 / (2.0 * OilMesh.a) * OilMesh.getDXU()[x - 1] * leta
                    * (muOil[x][OilMesh.getNY()][1] + mutOil[x][OilMesh.getNY()][1])
                    * (WOil.getNewField()[x][OilMesh.getNY()][1] - WOil.getNewField()[x][OilMesh.getNY() + 1][1])
                    / yp;
        }
//      K水相
        for (int x = 1; x < KWater.getNewField().length - 1; ++x) {
            leta = sqrt(WaterMesh.alpha(WaterMesh.gettPx()[x], WaterMesh.gettPy()[WaterMesh.getNY()]));
            yp = Tool.distance(
                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]),
                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]),
                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1])
            );
            Tauw
                    = Tauw + 1.0 / (2.0 * WaterMesh.a) * WaterMesh.getDXU()[x - 1] * leta
                    * (muWater[x][1][1] + mutWater[x][1][1])
                    * (WWater.getNewField()[x][1][1] - WWater.getNewField()[x][0][1])
                    / yp;
        }
    }

    /**
     * 计算粗糙度
     *
     * @param Tau
     * @param rho
     * @param x
     * @param y
     * @param z
     * @return
     */
    double cacFintES(double Tau, double rho, int x, int y, int z) {
        double hw = this.hl;
        double H = hw / this.R - 1.0;
        double Si = 2 * this.R * sqrt(1.0 - H * H);
        double So = 2 * this.R * Math.acos(H);
        double Sw = 2 * this.R - So;
        double Ao = this.R * 2.0 / 4.0 * (So - Si * H);
        double Aw = PI * this.R * this.R - Ao;
        double Dint, Reint;
        if (wso > wsw) {
            Dint = 4.0 * Ao / (So + Si);
            Reint = rhoOil[x][OilMesh.getNY()][z] * wso * Dint / muOil[x][OilMesh.getNY()][z];
        } else if (wsw > wso) {
            Dint = 4.0 * Aw / (Sw + Si);
            Reint = rhoWater[x][1][z] * wsw * Dint / muWater[x][1][z];
        } else {
            Dint = 4.0 * Aw / Sw;
            Reint = rhoWater[x][1][z] * wsw * Dint / muWater[x][1][z];
        }
        double fint = cacFint(Tau, rho, wso, wsw);
        double es = cacEs(abs(fint), Reint, Dint);
        return es;
    }

    /**
     * 计算能量耗散率场
     */
    void cacOmega() {
        double URF = 0.8;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double sigma1 = 0.5;
        double alpha = 0.555;
        double beta = 0.075;
        Label flagW = new Label("Omega");
        for (int z = 1; z < Omega.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < Omega.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < Omega.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagW.setFlag(Omega.getBoundary(), x, y);
                    double vol
                            = mesh.getDXU()[x - 1] * mesh.getDYV()[y - 1]
                            * mesh.J(mesh.gettPx()[x], mesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (mu[x][y][z] + mut[x][y][z] * sigma1);
                    muw = (mu[x - 1][y][z] + mut[x - 1][y][z] * sigma1);
                    mue = (mu[x + 1][y][z] + mut[x + 1][y][z] * sigma1);
                    mus = (mu[x][y - 1][z] + mut[x][y - 1][z] * sigma1);
                    mun = (mu[x][y + 1][z] + mut[x][y + 1][z] * sigma1);
//得到gamma
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[x - 1], mesh.getTx()[x - 1], mesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[x], mesh.getTx()[x], mesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[y - 1], mesh.getTy()[y - 1], mesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[y], mesh.getTy()[y], mesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * mesh.getDYV()[y - 1] / mesh.getDXP()[x - 1];
                    De = game * mesh.getDYV()[y - 1] / mesh.getDXP()[x];
                    Ds = gams * mesh.getDXU()[x - 1] / mesh.getDYP()[y - 1];
                    Dn = gamn * mesh.getDXU()[x - 1] / mesh.getDYP()[y];
//系数
                    coeO.aw[x][y][z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeO.ae[x][y][z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeO.as[x][y][z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeO.an[x][y][z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -Omega.getOldField()[x][y][z] * rho[x][y][z] * beta;
                    Spad
                            = flagW.w
                            * Spad1(Omega.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagW.e
                            * Spad2(Omega.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagW.s
                            * Spad1(Omega.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagW.n
                            * Spad2(Omega.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, 1);
//得到ap系数
                    coeO.ap[x][y][z]
                            = coeO.aw[x][y][z] + coeO.ae[x][y][z]
                            + coeO.as[x][y][z] + coeO.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs)
                            - (Sp + Spad) * vol;
//得到源项
                    Sc
                            = alpha * (mut[x][y][z] * S[x][y][z])
                            * Omega.getOldField()[x][y][z] / K.getOldField()[x][y][z];
                    Scad
                            = flagW.w
                            * Scad1(Omega.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, Omega.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagW.e
                            * Scad2(Omega.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, Omega.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagW.s
                            * Scad1(Omega.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, Omega.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagW.n
                            * Scad2(Omega.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, Omega.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeO.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>

                }
            }
        }
        double yp = 0;
        for (int z = 1; z < Omega.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < Omega.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < Omega.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    if ((x == 1) || (x == Omega.getNewField().length - 2)
                            || (y == 1) || (y == Omega.getNewField()[0].length - 2)) {
                        if (x == 1) {
                            yp = Tool.distance(
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realX(mesh.gettPx()[x - 1], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x - 1], mesh.gettPy()[y]));
                        }
                        if (x == Omega.getNewField().length - 2) {
                            yp = Tool.distance(
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realX(mesh.gettPx()[x + 1], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x + 1], mesh.gettPy()[y]));
                        }
                        if (y == 1) {
                            yp = Tool.distance(
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y - 1]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y - 1]));

                        }
                        if (y == Omega.getNewField()[0].length - 2) {
                            yp = Tool.distance(
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y]),
                                    mesh.realX(mesh.gettPx()[x], mesh.gettPy()[y + 1]),
                                    mesh.realY(mesh.gettPx()[x], mesh.gettPy()[y + 1]));
                        }
                        coeO.ap[x][y][z] = 1.0;
                        coeO.aw[x][y][z] = 0;//w
                        coeO.ae[x][y][z] = 0;//e
                        coeO.as[x][y][z] = 0;//s
                        coeO.an[x][y][z] = 0;//n
                        coeO.b[x][y][z]
                                = 6.0 * mu[x][y][z] / (rho[x][y][z] * beta * yp * yp);
                    }
                    //</editor-fold>
                    //under-relaxization
                    coeO.b[x][y][z]
                            = coeO.b[x][y][z] + (URFFI - 1.0) * coeO.ap[x][y][z]
                            * Omega.getOldField()[x][y][z];
                    coeO.ap[x][y][z] = (coeO.ap[x][y][z]) * URFFI;
                }
            }
        }

        Omega.getLog().rl = MatrixSolver.gaussSeidel(Omega.getNewField(), coeO, sys);
        Omega.getLog().log("Omega");
    }

    /**
     * 计算油相Omega
     */
    void cacOmegaOil() {
        double URF = 0.8;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double sigma1 = 0.5;
        double alpha = 0.555;
        double beta = 0.075;
        Label flagW = new Label("OmegaOil");
        for (int z = 1; z < OmegaOil.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < OmegaOil.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < OmegaOil.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagW.setFlag(OmegaOil.getBoundary(), x, y);
                    double vol
                            = OilMesh.getDXU()[x - 1] * OilMesh.getDYV()[y - 1]
                            * OilMesh.J(OilMesh.gettPx()[x], OilMesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (muOil[x][y][z] + mutOil[x][y][z] * sigma1);
                    muw = (muOil[x - 1][y][z] + mutOil[x - 1][y][z] * sigma1);
                    mue = (muOil[x + 1][y][z] + mutOil[x + 1][y][z] * sigma1);
                    mus = (muOil[x][y - 1][z] + mutOil[x][y - 1][z] * sigma1);
                    mun = (muOil[x][y + 1][z] + mutOil[x][y + 1][z] * sigma1);
//得到gamma
                    gamw = faceValue(muw, mup - muw, OilMesh.gettPx()[x - 1], OilMesh.getTx()[x - 1], OilMesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, OilMesh.gettPx()[x], OilMesh.getTx()[x], OilMesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, OilMesh.gettPy()[y - 1], OilMesh.getTy()[y - 1], OilMesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, OilMesh.gettPy()[y], OilMesh.getTy()[y], OilMesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * OilMesh.getDYV()[y - 1] / OilMesh.getDXP()[x - 1];
                    De = game * OilMesh.getDYV()[y - 1] / OilMesh.getDXP()[x];
                    Ds = gams * OilMesh.getDXU()[x - 1] / OilMesh.getDYP()[y - 1];
                    Dn = gamn * OilMesh.getDXU()[x - 1] / OilMesh.getDYP()[y];
//系数
                    coeOoil.aw[x][y][z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeOoil.ae[x][y][z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeOoil.as[x][y][z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeOoil.an[x][y][z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -OmegaOil.getOldField()[x][y][z] * rhoOil[x][y][z] * beta;
                    Spad
                            = flagW.w
                            * Spad1(OmegaOil.getBoundary().getType(x - 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagW.e
                            * Spad2(OmegaOil.getBoundary().getType(x + 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagW.s
                            * Spad1(OmegaOil.getBoundary().getType(x, y - 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagW.n
                            * Spad2(OmegaOil.getBoundary().getType(x, y + 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y],
                                    gamn, Fn, 1);
//得到ap系数
                    coeOoil.ap[x][y][z]
                            = coeOoil.aw[x][y][z] + coeOoil.ae[x][y][z]
                            + coeOoil.as[x][y][z] + coeOoil.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc
                            = alpha * (mutOil[x][y][z] * Soil[x][y][z])
                            * OmegaOil.getOldField()[x][y][z] / KOil.getOldField()[x][y][z];
                    Scad
                            = flagW.w
                            * Scad1(OmegaOil.getBoundary().getType(x - 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x - 1],
                                    gamw, Fw, OmegaOil.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagW.e
                            * Scad2(OmegaOil.getBoundary().getType(x + 1, y, z),
                                    OilMesh.getDYV()[y - 1], vol, OilMesh.getDXP()[x],
                                    game, Fe, OmegaOil.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagW.s
                            * Scad1(OmegaOil.getBoundary().getType(x, y - 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y - 1],
                                    gams, Fs, OmegaOil.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagW.n
                            * Scad2(OmegaOil.getBoundary().getType(x, y + 1, z),
                                    OilMesh.getDXU()[x - 1], vol, OilMesh.getDYP()[y],
                                    gamn, Fn, OmegaOil.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeOoil.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
                }
            }
        }
        double DN = 1E-30;
        double YPL1, YP, epl, es, Uint;
        for (int z = 1; z < OmegaOil.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < OmegaOil.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < OmegaOil.getNewField().length - 1; ++x) {

                    //<editor-fold>
                    flagW.setFlag(OmegaOil.getBoundary(), x, y);
                    if ((flagW.w == 1) || (flagW.e == 1)
                            || (flagW.s == 1) || (flagW.n == 1)) {
                        if (flagW.w == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x - 1], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x - 1], OilMesh.gettPy()[y]));
                        }
                        if (flagW.e == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x + 1], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x + 1], OilMesh.gettPy()[y]));
                        }
                        if (flagW.s == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y - 1]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y - 1]));

                        }
                        if (flagW.n == 1) {
                            DN = Tool.distance(
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y]),
                                    OilMesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[y + 1]),
                                    OilMesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[y + 1]));

                            Uint = sqrt(abs(Tauo) / rhoOil[x][y][z]);
                            es = cacFintES(Tauo, rhoOil[x][y][z], x, y, z);
                            epl = rhoOil[x][y][z] * es * Uint / muOil[x][y][z];
                            if (epl <= 4.535) {
                                YPL1 = 0;
                            } else {
                                YPL1 = max(0, 0.9 * sqrt(epl) - epl * exp(-epl / 6.0));
                            }
                            YP = YPL1 * muOil[x][y][z] / rhoOil[x][y][z] / Uint;
                            DN = DN + YP;
                        }
                        coeOoil.ap[x][y][z] = 1.0;
                        coeOoil.aw[x][y][z] = 0;//w
                        coeOoil.ae[x][y][z] = 0;//e
                        coeOoil.as[x][y][z] = 0;//s
                        coeOoil.an[x][y][z] = 0;//n
                        coeOoil.b[x][y][z]
                                = 6.0 * muOil[x][y][z] / (rhoOil[x][y][z] * beta * DN * DN);
                    }
                    //</editor-fold>
                    //under-relaxization
                    coeOoil.b[x][y][z]
                            = coeOoil.b[x][y][z] + (URFFI - 1.0) * coeOoil.ap[x][y][z]
                            * OmegaOil.getOldField()[x][y][z];
                    coeOoil.ap[x][y][z] = (coeOoil.ap[x][y][z]) * URFFI;
                }
            }
        }
        OmegaOil.getLog().rl = MatrixSolver.gaussSeidel(OmegaOil.getNewField(), coeOoil, sysOil);
        OmegaOil.getLog().log("OmegaOil");
    }

    /**
     * 计算水相Omega
     */
    void cacOmegaWater() {
        double URF = 0.8;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double sigma1 = 0.5;
        double alpha = 0.555;
        double beta = 0.075;
        Label flagW = new Label("OmegaWater");
        for (int z = 1; z < OmegaWater.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < OmegaWater.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < OmegaWater.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagW.setFlag(OmegaWater.getBoundary(), x, y);
                    double vol
                            = WaterMesh.getDXU()[x - 1] * WaterMesh.getDYV()[y - 1]
                            * WaterMesh.J(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (muWater[x][y][z] + mutWater[x][y][z] * sigma1);
                    muw = (muWater[x - 1][y][z] + mutWater[x - 1][y][z] * sigma1);
                    mue = (muWater[x + 1][y][z] + mutWater[x + 1][y][z] * sigma1);
                    mus = (muWater[x][y - 1][z] + mutWater[x][y - 1][z] * sigma1);
                    mun = (muWater[x][y + 1][z] + mutWater[x][y + 1][z] * sigma1);
//得到gamma
                    gamw = faceValue(muw, mup - muw, WaterMesh.gettPx()[x - 1], WaterMesh.getTx()[x - 1], WaterMesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, WaterMesh.gettPx()[x], WaterMesh.getTx()[x], WaterMesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, WaterMesh.gettPy()[y - 1], WaterMesh.getTy()[y - 1], WaterMesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, WaterMesh.gettPy()[y], WaterMesh.getTy()[y], WaterMesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * WaterMesh.getDYV()[y - 1] / WaterMesh.getDXP()[x - 1];
                    De = game * WaterMesh.getDYV()[y - 1] / WaterMesh.getDXP()[x];
                    Ds = gams * WaterMesh.getDXU()[x - 1] / WaterMesh.getDYP()[y - 1];
                    Dn = gamn * WaterMesh.getDXU()[x - 1] / WaterMesh.getDYP()[y];
//系数
                    coeOwater.aw[x][y][z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeOwater.ae[x][y][z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeOwater.as[x][y][z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeOwater.an[x][y][z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = -OmegaWater.getOldField()[x][y][z] * rhoWater[x][y][z] * beta;
                    Spad
                            = flagW.w
                            * Spad1(OmegaWater.getBoundary().getType(x - 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagW.e
                            * Spad2(OmegaWater.getBoundary().getType(x + 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagW.s
                            * Spad1(OmegaWater.getBoundary().getType(x, y - 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagW.n
                            * Spad2(OmegaWater.getBoundary().getType(x, y + 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y],
                                    gamn, Fn, 1);
//得到ap系数
                    coeOwater.ap[x][y][z]
                            = coeOwater.aw[x][y][z] + coeOwater.ae[x][y][z]
                            + coeOwater.as[x][y][z] + coeOwater.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc
                            = alpha * (mutWater[x][y][z] * Swater[x][y][z])
                            * OmegaWater.getOldField()[x][y][z] / KWater.getOldField()[x][y][z];
                    Scad
                            = flagW.w
                            * Scad1(OmegaWater.getBoundary().getType(x - 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x - 1],
                                    gamw, Fw, OmegaWater.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagW.e
                            * Scad2(OmegaWater.getBoundary().getType(x + 1, y, z),
                                    WaterMesh.getDYV()[y - 1], vol, WaterMesh.getDXP()[x],
                                    game, Fe, OmegaWater.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagW.s
                            * Scad1(OmegaWater.getBoundary().getType(x, y - 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y - 1],
                                    gams, Fs, OmegaWater.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagW.n
                            * Scad2(OmegaWater.getBoundary().getType(x, y + 1, z),
                                    WaterMesh.getDXU()[x - 1], vol, WaterMesh.getDYP()[y],
                                    gamn, Fn, OmegaWater.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeOwater.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
                }
            }
        }
        double DN = 1E-30;
        double Uint, es, epl, YPL1, YP;
        for (int z = 1; z < OmegaWater.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < OmegaWater.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < OmegaWater.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    if ((x == 1) || (x == OmegaWater.getNewField().length - 2)
                            || (y == 1) || (y == OmegaWater.getNewField()[0].length - 2)) {
                        if (x == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x - 1], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x - 1], WaterMesh.gettPy()[y]));
                        }
                        if (x == OmegaWater.getNewField().length - 2) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x + 1], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x + 1], WaterMesh.gettPy()[y]));
                        }
                        if (y == 1) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y - 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y - 1]));

                            Uint = sqrt(abs(Tauw) / rhoWater[x][y][z]);
                            es = cacFintES(Tauw, rhoWater[x][y][z], x, y, z);
                            epl = rhoWater[x][y][z] * es * Uint / muWater[x][y][z];
                            if (epl <= 4.535) {
                                YPL1 = 0;
                            } else {
                                YPL1 = max(0, 0.9 * sqrt(epl) - epl * exp(-epl / 6.0));
                            }
                            YP = YPL1 * muWater[x][y][z] / rhoWater[x][y][z] / Uint;
                            DN = DN + YP;
                        }
                        if (y == OmegaWater.getNewField()[0].length - 2) {
                            DN = Tool.distance(
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y]),
                                    WaterMesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y + 1]),
                                    WaterMesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[y + 1]));
                        }
                        coeOwater.ap[x][y][z] = 1.0;
                        coeOwater.aw[x][y][z] = 0;//w
                        coeOwater.ae[x][y][z] = 0;//e
                        coeOwater.as[x][y][z] = 0;//s
                        coeOwater.an[x][y][z] = 0;//n
                        coeOwater.b[x][y][z]
                                = 6.0 * muWater[x][y][z] / (rhoWater[x][y][z] * beta * DN * DN);
                    }
                    //</editor-fold>
                    //under-relaxization
                    coeOwater.b[x][y][z]
                            = coeOwater.b[x][y][z] + (URFFI - 1.0) * coeOwater.ap[x][y][z]
                            * OmegaWater.getOldField()[x][y][z];
                    coeOwater.ap[x][y][z] = (coeOwater.ap[x][y][z]) * URFFI;
                }
            }
        }
        OmegaWater.getLog().rl = MatrixSolver.gaussSeidel(OmegaWater.getNewField(), coeOwater, sysWater);
        OmegaWater.getLog().log("OmegaWater");
    }

    /**
     * 计算速度场
     *
     * @param dpdz
     */
    void cacW(double dpdz, Mesh mesh, Field W,
            double[][][] mu, double[][][] rho,
            double[][][] Mut, SystemControl sys
    ) {
        double URF = 0.9;
        double URFFI = 1. / URF;
        double Fe, Fw, Fn, Fs;
        double De, Dw, Dn, Ds;
        double Spad, Scad, Sp, Sc;
        double muw, mue, mus, mun, mup;
        double gamw, game, gams, gamn;
        double alpha = 1.0;
        Label flagW = new Label("W");
        int i = mesh.numPx();
        int j = mesh.numPy();
        Coefficient coeW = new Coefficient(i, j);
        for (int z = 1; z < W.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < W.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < W.getNewField().length - 1; ++x) {
                    //<editor-fold>
                    flagW.setFlag(W.getBoundary(), x, y);
                    double vol
                            = mesh.getDXU()[x - 1] * mesh.getDYV()[y - 1]
                            * mesh.J(mesh.gettPx()[x], mesh.gettPy()[y]);
                    Fw = 0;
                    Fe = 0;
                    Fs = 0;
                    Fn = 0;
//
                    mup = (mu[x][y][z] + Mut[x][y][z] / alpha);
                    muw = (mu[x - 1][y][z] + Mut[x - 1][y][z] / alpha);
                    mue = (mu[x + 1][y][z] + Mut[x + 1][y][z] / alpha);
                    mus = (mu[x][y - 1][z] + Mut[x][y - 1][z] / alpha);
                    mun = (mu[x][y + 1][z] + Mut[x][y + 1][z] / alpha);
//得到gamma
                    gamw = faceValue(muw, mup - muw, mesh.gettPx()[x - 1], mesh.getTx()[x - 1], mesh.gettPx()[x]);
                    game = faceValue(mup, mue - mup, mesh.gettPx()[x], mesh.getTx()[x], mesh.gettPx()[x + 1]);
                    gams = faceValue(mus, mup - mus, mesh.gettPy()[y - 1], mesh.getTy()[y - 1], mesh.gettPy()[y]);
                    gamn = faceValue(mup, mun - mup, mesh.gettPy()[y], mesh.getTy()[y], mesh.gettPy()[y + 1]);
//GET d
                    Dw = gamw * mesh.getDYV()[y - 1] / mesh.getDXP()[x - 1];
                    De = game * mesh.getDYV()[y - 1] / mesh.getDXP()[x];
                    Ds = gams * mesh.getDXU()[x - 1] / mesh.getDYP()[y - 1];
                    Dn = gamn * mesh.getDXU()[x - 1] / mesh.getDYP()[y];
//系数
                    coeW.aw[x][y][z] = (1 - flagW.w) * (Math.max(Fw, 0.0) + Dw);//w
                    coeW.ae[x][y][z] = (1 - flagW.e) * (Math.max(-Fe, 0.0) + De);//e
                    coeW.as[x][y][z] = (1 - flagW.s) * (Math.max(Fs, 0.0) + Ds);//s
                    coeW.an[x][y][z] = (1 - flagW.n) * (Math.max(-Fn, 0.0) + Dn);//n
//设定附加源项
                    Sp = 0;
                    Spad
                            = flagW.w
                            * Spad1(W.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, 1)
                            + flagW.e
                            * Spad2(W.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, 1)
                            + flagW.s
                            * Spad1(W.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, 1)
                            + flagW.n
                            * Spad2(W.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, 1);
//得到ap系数
                    coeW.ap[x][y][z]
                            = coeW.aw[x][y][z] + coeW.ae[x][y][z]
                            + coeW.as[x][y][z] + coeW.an[x][y][z]
                            + (Fe - Fw) + (Fn - Fs) - (Sp + Spad) * vol;
//得到源项
                    Sc = dpdz - rho[x][y][z] * g * Math.sin(0);
                    Scad
                            = flagW.w
                            * Scad1(W.getBoundary().getType(x - 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x - 1],
                                    gamw, Fw, W.getNewField()[x - 1][y][z], 0, 1, 1)
                            + flagW.e
                            * Scad2(W.getBoundary().getType(x + 1, y, z),
                                    mesh.getDYV()[y - 1], vol, mesh.getDXP()[x],
                                    game, Fe, W.getNewField()[x + 1][y][z], 0, 1, 1)
                            + flagW.s
                            * Scad1(W.getBoundary().getType(x, y - 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y - 1],
                                    gams, Fs, W.getNewField()[x][y - 1][z], 0, 1, 1)
                            + flagW.n
                            * Scad2(W.getBoundary().getType(x, y + 1, z),
                                    mesh.getDXU()[x - 1], vol, mesh.getDYP()[y],
                                    gamn, Fn, W.getNewField()[x][y + 1][z], 0, 1, 1);
                    coeW.b[x][y][z] = (Sc + Scad) * vol;
                    //</editor-fold>
//under-relaxization
                    coeW.b[x][y][z]
                            = coeW.b[x][y][z] + (URFFI - 1.0) * coeW.ap[x][y][z]
                            * W.getOldField()[x][y][z];
                    coeW.ap[x][y][z] = (coeW.ap[x][y][z]) * URFFI;
                }
            }
        }
        W.getLog().rl = MatrixSolver.gaussSeidel(W.getNewField(), coeW, sys);
        W.getLog().log("W");
    }

    /**
     * 计算湍流粘度
     *
     * @param rho
     * @param Mut
     * @param K
     * @param dpdz 压降
     * @param omega
     */
    public void cacMut(double[][][] rho, double[][][] Mut, Field K, Field omega) {
        double alpha = 4.0;
        for (int z = 1; z < K.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < K.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < K.getNewField().length - 1; ++x) {
                    Mut[x][y][z]
                            = alpha * rho[x][y][z] * K.getNewField()[x][y][z]
                            / (omega.getNewField()[x][y][z] + 1e-30);
                }
            }
        }
    }

    /**
     * 引入温度因素改变粘度
     */
    void cacMu() {
        double K = 5.0e-3, n = 0.5;
        for (int z = 1; z < WOil.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < WOil.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < WOil.getNewField().length - 1; ++x) {
                    if (TOil.getNewField()[x][y][z] > 10) {
                        muOil[x][y][z]
                                = 0.058174
                                * exp(-0.0371
                                        * (TOil.getNewField()[x][y][z] - 273.15));
                    } else {
                        muOil[x][y][z] = K * pow(Soil[x][y][z], n - 1.0);
                    }
                }
            }
        }
    }

    public komegaTFOWAngeli1(double Qso, double Qsw, double Radia) {
        this.Qso = Qso;
        this.Qsw = Qsw;
        this.R = Radia;
    }

    public static void main(String[] args) throws FileNotFoundException {
        double dpdz = 252.0, yint = 0.50;
        double Radia = 0.0243 / 2.0;
        double Qwater = 0.22 * PI * Radia * Radia;
        double Qoil = 0.32 * PI * Radia * Radia;
//        
        double F, G;
        double Fy, Fx, Gy, Gx;
        double xold = dpdz, yold = yint;
        double xnew, ynew;
        double errF, errG;

        String dirName = "/D:/winsway/" + "komegaTFOWAngeli-othermut/";
        String title = "pressure and liquid high";
        Writer.createDir(dirName);
        String fileName = dirName + title + "." + "txt";
        Writer.createFile(fileName);
        PrintWriter pw
                = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        pw.printf("Qwater = %e\t Qoil = %e\t\n", Qwater, Qoil);
        pw.flush();
        do {
            System.out.println("Qwater = " + Qwater + " Qoil = " + Qoil);
            System.out.println("F1   " + "xold =" + xold + " yold =" + yold);
            pw.printf("F1 xold =%e\t yold =%e\t\n", xold, yold);
            pw.flush();
            komegaTFOWAngeli1 F1 = new komegaTFOWAngeli1(Qoil, Qwater, Radia);
            F1.application(xold, yold);
            F = F1.F(Qwater);
            G = F1.G(Qoil);
            pw.printf("Qcwater = %e\t Qcoil = %e\t\n", F + Qwater, G + Qoil);
            pw.flush();
            System.out.println("Fx1   " + "xold * 1.01 =" + xold * 1.01 + " yold =" + yold);
            pw.printf("Fx1 xold * 1.01 =%e\t yold =%e\t\n", xold * 1.01, yold);
            pw.flush();
            komegaTFOWAngeli1 Fx1 = new komegaTFOWAngeli1(Qoil, Qwater, Radia);
            Fx1.application(xold * 1.01, yold);
            System.out.println("Fy1   " + " xold =" + xold + " yold * 1.01 =" + yold * 1.01);
            pw.printf("Fy1 xold =%e\t yold * 1.01 =%e\t\n", xold, yold * 1.01);
            pw.flush();
            komegaTFOWAngeli1 Fy1 = new komegaTFOWAngeli1(Qoil, Qwater, Radia);
            Fy1.application(xold, yold * 1.01);
//            

            Fy = (Fy1.F(Qwater) - F) / (0.01 * yold);
            Gy = (Fy1.G(Qoil) - G) / (0.01 * yold);
            Fx = (Fx1.F(Qwater) - F) / (0.01 * xold);
            Gx = (Fx1.G(Qoil) - G) / (0.01 * xold);
            xnew = xold + (G * Fy - F * Gy) / (Fx * Gy - Gx * Fy);
            ynew = yold + (F * Gx - G * Fx) / (Fx * Gy - Gx * Fy);
            xold = xnew;
            yold = ynew;
            System.out.println("xnew = " + xnew + " ynew = " + ynew);
            errF = abs((F) / Qwater);
            errG = abs((G) / Qoil);
            System.out.println("errF = " + errF + " errG = " + errG);
        } while (errF >= 1e-3 || errG >= 1e-3);
        pw.close();
    }

    public boolean controlDict() {
        this.controldict = new ControlDict();
        this.controldict.deltaT = 1e-3;
        this.controldict.endTime = 40;
        this.controldict.writeInterval = 1000;
        this.controldict.position = "/D:/winsway/" + this.toString();
        return true;
    }

    public boolean fvSolution() {
        this.fvsolution = new FvSolution();
        return true;
    }

    @Override
    public String toString() {
        return "komegaTFOWAngeli-othermut";
    }

    public void application(double dpdz, double hl) {
        this.dpdz = dpdz;
        this.hl = hl;
        blockMesh();
        createFields();
        controlDict();
        fvSolution();
        sys = new SystemControl(mesh, controldict, fvsolution);
        sysOil = new SystemControl(OilMesh, controldict, fvsolution);
        sysWater = new SystemControl(WaterMesh, controldict, fvsolution);
        initField();
        Solver();
        System.out.println("求解完成！");
    }

    /**
     * 设定湍流动能
     *
     * @param K
     * @param W
     */
    public void setK(Field K, Field W) {
        for (int z = 1; z < K.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < K.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < K.getNewField().length - 1; ++x) {
                    double value = 1e-4 * pow(W.getNewField()[x][y][z], 2.0);
                    K.setNewField(x, y, z, value);
                }
            }
        }
    }

    /**
     * 设定能量耗散率
     */
    public void setOmega() {
        for (int z = 1; z < Omega.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < Omega.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < Omega.getNewField().length - 1; ++x) {
                    double value
                            = pow(0.09, 0.75)
                            * pow(K.getNewField()[x][y][z], 3.0 / 2.0)
                            / Turbulence.turbulenceLengthscale(bm.radia * 2);
                    value = value / K.getNewField()[x][y][z];
                    Omega.setNewField(x, y, z, value);
                }
            }
        }
    }

    /**
     * 设定能量耗散率
     *
     * @param omega
     * @param K
     * @param bm
     */
    public void setOmega(Field omega, Field K, BiSingleV1 bm) {
        for (int z = 1; z < omega.getNewField()[0][0].length - 1; ++z) {
            for (int y = 1; y < omega.getNewField()[0].length - 1; ++y) {
                for (int x = 1; x < omega.getNewField().length - 1; ++x) {
                    double value
                            = pow(0.09, 0.75)
                            * pow(K.getNewField()[x][y][z], 3.0 / 2.0)
                            / Turbulence.turbulenceLengthscale(bm.radia * 2);
                    value = value / K.getNewField()[x][y][z];
                    omega.setNewField(x, y, z, value);
                }
            }
        }
    }

    /**
     * 更新界面速度
     */
    void updateInterfaceW() {
        for (int x = 1; x < W.getNewField().length - 1; ++x) {
            double yOil = Tool.distance(
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY() + 1]),
                    mesh.realX(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]),
                    mesh.realY(OilMesh.gettPx()[x], OilMesh.gettPy()[OilMesh.getNY()]));
            double yWater = Tool.distance(
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[1]),
                    mesh.realX(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]),
                    mesh.realY(WaterMesh.gettPx()[x], WaterMesh.gettPy()[0]));
            double w = (muWater[x][1][1] + mutWater[x][1][1]) / yWater;
            double o = (muOil[x][OilMesh.getNY()][1] + mutOil[x][OilMesh.getNY()][1]) / yOil;
            double u = (WOil.getNewField()[x][OilMesh.getNY()][1] * o + WWater.getNewField()[x][1][1] * w)
                    / (o + w);
            WOil.setNewField(x, OilMesh.getNY() + 1, 1, u);
            WWater.setNewField(x, 0, 1, u);
        }
    }

    void outPut(int count, Boolean judge, String var, double[][][] F, Mesh mesh) {
        if (count % 1000 == 0 || Objects.equals(judge, Boolean.TRUE)) {
            String dir = sys.conDict.position + "/" + var + "/";
            String title = var + "_" + count;
            Output write = new Output(F, mesh, dir, title);
        }
    }

    double cacEs(double fint, double Reint, double Dint) {
        double A = sqrt(1.0 / (fint + 1e-30));
        double B = pow(10, -A / 3.6);
        double C = 3.7 * Dint;
        B = max(0, B - 6.9 / Reint);
        double temp = C * pow(B, 1.0 / 1.11);
        return temp;
    }

    double cacFint(double tau, double rho, double wso, double wsw) {
        return 2.0 * tau / (rho * (wso - wsw) * abs((wso - wsw)));
    }

}
