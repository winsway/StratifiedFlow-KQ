package com.winswe.io;

import com.alibaba.fastjson.JSON;
import com.alibaba.fastjson.JSONObject;
import com.winswe.field.VolScalarField;
import com.winswe.mesh.Structed2D;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.rmi.StubNotFoundException;
import java.util.ArrayList;
import java.util.List;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月12日 下午3:22:22
 */
public class IOobject {

    /**
     * 读取JSON配置文件
     */
    private JSONObject jsonObject;

    /**
     * 网格输出
     */
    private IOMesh ioMesh;

    /**
     * 场输出
     */
    private IOField ioField;

    /**
     * 路径
     */
    private final String path;

    /**
     * 案例名字
     */
    private final String caseName;

    /**
     *
     */
    private final List<VolScalarField> fields = new ArrayList<>();

    public IOobject(String path, String caseName) {
        this.path = path;
        this.caseName = caseName;
    }

    /**
     * 读取配置文件
     */
    public void readConfigure() {
        String systemPath = path + "/" + caseName + "/" + "system" + "/" + "controlDict.json";
        String temp = readJsonFile(systemPath);
        //
        jsonObject = JSON.parseObject(temp);
        //
        System.out.println(path + "/" + caseName);
    }

    /**
     * 读取json文件
     *
     * @param fileName 包含路径的文件名
     * @return 读取的字符串
     */
    public String readJsonFile(String fileName) {
        String jsonStr = "";
        try {
            //创建一个文件
            File jsonFile = new File(fileName);
            //根据文件创建文件读取器
            FileReader fileReader = new FileReader(jsonFile);
            //读取文件
            Reader reader = new InputStreamReader(new FileInputStream(jsonFile), "utf-8");
            int ch = 0;
            //建立字符串缓存
            StringBuffer sb = new StringBuffer();
            while ((ch = reader.read()) != -1) {
                sb.append((char) ch);
            }
            //关闭文件读取器
            fileReader.close();
            //
            reader.close();
            //转化字符串
            jsonStr = sb.toString();
            //返回字符串
            return jsonStr;
        } catch (IOException e) {
            //输出错误信息
            e.printStackTrace();
            return null;
        }
    }

    /**
     *
     * @param mesh 网格
     */
    public void outPutMesh(Structed2D mesh) {

        JSONObject meshJsonObject = jsonObject.getJSONObject("IO").getJSONObject("Mesh");
        ioMesh = new IOMesh(path + "/" + caseName + "/" + "constant", mesh, meshJsonObject);
        ioMesh.outPutMesh();
    }

    public void outPutField() {
        ioField = new IOField(path + "/" + caseName, fields);
        ioField.outPutField();
    }

    public JSONObject getJsonObject() {
        return jsonObject;
    }

    public IOMesh getIoMesh() {
        return ioMesh;
    }

    public IOField getIoField() {
        return ioField;
    }

    public String getPath() {
        return path;
    }

    public String getCaseName() {
        return caseName;
    }

    public List<VolScalarField> getField() {
        return fields;
    }

}
