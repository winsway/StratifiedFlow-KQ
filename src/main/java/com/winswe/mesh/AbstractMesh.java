/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.winswe.mesh;

import com.winswe.mesh.element.Cell;
import com.winswe.mesh.element.Face;
import com.winswe.mesh.element.Node;
import com.winswe.mesh.factory.FvPatch;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author winswe <halo.winswe@gmail.com>
 * @date 2021年2月13日 上午9:57:01
 */
public abstract class AbstractMesh {

    /**
     * 网格节点
     */
    private List<Node> nodes = new ArrayList<>();

    /**
     * 网格面
     */
    private List<Face> faces = new ArrayList<>();

    /**
     * 网格
     */
    private List<Cell> cells = new ArrayList<>();

    /**
     * 边界面
     */
    private List<FvPatch> fvPatchs = new ArrayList<>();

}
