/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.io;

import com.alibaba.fastjson.JSONObject;
import com.winswe.field.VolScalarField;
import com.winswe.mesh.Structed2D;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午9:28:59
 */
public class IOField {

    /**
     * 文件路径
     */
    private final String path;

    /**
     * 网格类
     */
    private final Structed2D mesh;

    /**
     * 输出对象
     */
    private PrintWriter out;

    private List<VolScalarField> fields;

    public IOField(String path, List<VolScalarField> fields) {
        this.path = path;
        this.fields = fields;
        this.mesh = fields.get(0).getMesh();
    }

    public void outPutField() {
        String fileName;

        fileName = path + "/" + "field" + "." + "dat";
        Writer.createFile(fileName);
        try {
            out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(fileName)));
        } catch (FileNotFoundException ex) {
            Logger.getLogger(IOMesh.class.getName()).log(Level.SEVERE, null, ex);
        }
        out.println(this.Title());
        out.println(this.Variables());
        out.println(this.Zone());

        for (int J = 0; J <= mesh.getNY() + 1; ++J) {
            for (int I = 0; I <= mesh.getNX() + 1; ++I) {
                out.printf("%16.6E\t", mesh.X(mesh.getPointX1()[I], mesh.getPointX2()[J]));
                out.printf("%16.6E\t", mesh.Y(mesh.getPointX1()[I], mesh.getPointX2()[J]));
                int IJ = mesh.getCellIndex(I, J);
                for (VolScalarField i : fields) {
                    out.printf("%16.6E\t", i.getFI()[IJ]);
                }
                out.println();
            }
        }
        out.close();
    }

    private String Title() {
        return "Title=" + "\"" + "Fields" + "\"";
    }

    private String Variables() {
        String name = "Variables=\"X\"\"Y\"";
        for (VolScalarField i : fields) {
            name = name + "\"" + i.getName() + "\"";
        }
        return name;
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

}
