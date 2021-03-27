/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.io;

import com.alibaba.fastjson.JSONObject;
import com.winswe.mesh.Structed2D;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午9:28:49
 */
public class IOMesh {

    /**
     * 文件路径
     */
    private final String path;

    /**
     * 网格类
     */
    private final Structed2D mesh;

    /**
     * 网格配置文件
     */
    private final JSONObject ioMesh;

    /**
     * 输出对象
     */
    private PrintWriter out;

    public IOMesh(String path, Structed2D mesh, JSONObject ioMesh) {
        this.path = path;
        this.mesh = mesh;
        this.ioMesh = ioMesh;
    }

    /**
     *
     * @return title名字
     */
    private String Title() {
        return "Title=" + "\"" + "Mesh" + "\"";
    }

    /**
     *
     * @return 变量
     */
    private String Variables() {
        return "Variables=\"X\",\"Y\"";
    }

    /**
     *
     * @return 域
     */
    private String Zone() {
        return "Zone"
                + " I=" + (mesh.getNX() + 2)
                + " J=" + (mesh.getNY() + 2)
                + " F=POINT";
    }

    /**
     * 输出物理几何
     *
     * @param out
     * @param mesh
     */
    private void OutPutPhysicalGeo(PrintWriter out, Structed2D mesh) {
        out.println(this.Title());
        out.println(this.Variables());
        out.println(this.Zone());
        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                out.printf("%16.6E\t", mesh.X(mesh.getPointX1()[I], mesh.getPointX2()[J]));
                out.printf("%16.6E\t", mesh.Y(mesh.getPointX1()[I], mesh.getPointX2()[J]));
                out.println();
            }
        }
        out.close();
    }

    /**
     * 输出计算几何
     *
     * @param out
     * @param mesh
     */
    private void OutPutCalculateGeo(PrintWriter out, Structed2D mesh) {
        out.println(this.Title());
        out.println(this.Variables());
        out.println(this.Zone());
        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                out.printf("%16.6E\t", mesh.getPointX1()[I]);
                out.printf("%16.6E\t", mesh.getPointX2()[J]);
                out.println();
            }
        }
        out.close();
    }

    /**
     * 输出网格
     */
    public void outPutMesh() {
        String fileName;
        if (ioMesh.getBoolean("PhysicalGeo")) {
            fileName = path + "/" + "PhysicalGeo" + "." + "dat";
            Writer.createFile(fileName);
            try {
                out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
            } catch (FileNotFoundException ex) {
                Logger.getLogger(IOMesh.class.getName()).log(Level.SEVERE, null, ex);
            }
            this.OutPutPhysicalGeo(out, mesh);
        }

        //        
        if (ioMesh.getBoolean("PhysicalGeo")) {
            fileName = path + "/" + "CalculateGeo" + "." + "dat";
            Writer.createFile(fileName);
            try {
                out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
            } catch (FileNotFoundException ex) {
                Logger.getLogger(IOMesh.class.getName()).log(Level.SEVERE, null, ex);
            }
            this.OutPutCalculateGeo(out, mesh);
        }
    }

}
