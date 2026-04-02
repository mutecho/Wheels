# 概述

该项目用于从蒙特卡洛事件生成器产生的信息中提取出想要的粒子对（目前主要是pi+pi+对）的3D-femto的相空间（根据中心度，粒子对mt和粒子对的phi-psi切片）；随后再拟合得到六个$R^2$随着phi-psi角的变化图。



## 语言

C++

## 输入（暂定假设）

一个ROOT文件，里面包含了许多事件；每个事件下产生的粒子又有对应的能动量和时空信息，以及对应的粒子种类（PDGCode）



## 输出

一个ROOT文件，内部以中心度和mt和phi为目录，每个目录下有对应切片的3D-femto的相空间图，六个方向上相空间的一维投影图和对应的拟合图，以及拟合参数，不确定度等数据。文件内另外的一个单独$R ^2$文件夹以中心度和mt为区分，画出每个区间下的六个$R^2$随着phi-psi角的变化图。



## 注意事项

- 一些函数的写法可以参考本路径下其他的md文件

- fit部分逻辑按照我Zotero库中Femtoscopic signatures of unique nuclear structures in relativistic collisions这篇文章
- fit的代码可以部分参考/Users/allenzhou/ALICE/scripts/macros/femtoep/3d_cf_from_exp.cpp，不过该程序的fit方法和对象都和此处不同；可以参考的是输出的参数，形式等