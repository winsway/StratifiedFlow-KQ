/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.cup.io;

/**
 * 用于Case输出
 *
 * @author winsway
 */
public interface CaseWriter {

    default void describe() {
        System.out.println("用于文件描述！");
    }

}
