/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.io;

import com.alibaba.fastjson.JSONObject;
import org.junit.Test;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 */
public class IOobjectTest {

    public IOobjectTest() {
    }

    /**
     * Test of readConfigure method, of class IOobject.
     */
    @Test
    public void testReadConfigure() {
        System.out.println("readConfigure");
        String path = "./tutorials/io";
        String caseName = "test";
        IOobject instance = new IOobject(path, caseName);
        instance.readConfigure();
        //得到整体的Json对象
        JSONObject globalJsonObject = instance.getJsonObject();
        System.out.println("name json =" + globalJsonObject.toJSONString());
    }

}
