语言：C++

函数：

- 目标
  - 计算一对粒子对于3Dfemto源的贡献
- 输入
  - 两个粒子的px，py，pz，mass，x，y，z，t （double类型）。
- 输出
  - 该粒子对在LCMS坐标系下的3D-femto相空间三方向坐标：$(\rho_{\rm out},\rho_{\rm side},\rho_{\rm long})$，同样以double类型，可使用vector或者自定义类

- 过程
  - 通过两粒子的四动量获得该粒子对四动量（以及横动量等后续计算需要的量）
  - 根据粒子对四动量定义out，side，long方向，以及LCMS的beta和gamma
  - 将该粒子对的四位移投影到out，side，long方向，并boost到LCMS下
  - 得到$(\rho_{\rm out},\rho_{\rm side},\rho_{\rm long})$，并返回结果
- 注意事项
  - 在需要时（比如需要矢量操作时）调用ROOT库
  - 检查粒子对距离是否过近，如果过近则拒绝该对粒子（close pair rejection）
  - 不要使用高于C++20的语法
  - 如果对过程有不确定之处，或者发现其中存在错误（逻辑上，物理上等等）马上提出



以下代码仅供参考，功能与逻辑仅部分相同：

```C++
/// Compute the 3d components of the pair momentum in LCMS and PRF
  /// Copy from femto universe
  /// \tparam T type of tracks
  /// \param part1 Particle 1
  /// \param mass1 Mass of particle 1
  /// \param part2 Particle 2
  /// \param mass2 Mass of particle 2
  /// \param isiden Identical or non-identical particle pair
  template <typename T>
  static std::vector<double> newpairfunc(const T& part1, const float mass1, const T& part2, const float mass2, bool isiden){
    //计算能量
    const double e1 = std::sqrt(std::pow(part1.px(), 2) + std::pow(part1.py(), 2) + std::pow(part1.pz(), 2) + std::pow(mass1, 2));
    const double e2 = std::sqrt(std::pow(part2.px(), 2) + std::pow(part2.py(), 2) + std::pow(part2.pz(), 2) + std::pow(mass2, 2));
		//构建两粒子和粒子对的矢量
    const ROOT::Math::PxPyPzEVector vecpart1(part1.px(), part1.py(), part1.pz(), e1);
    const ROOT::Math::PxPyPzEVector vecpart2(part2.px(), part2.py(), part2.pz(), e2);
    const ROOT::Math::PxPyPzEVector trackSum = vecpart1 + vecpart2;
    
    //pair符号约定
    const double tPx = trackSum.px();
    const double tPy = trackSum.py();
    const double tPz = trackSum.pz();
    const double tE = trackSum.E();

    const double tPtSq = (tPx * tPx + tPy * tPy);
    const double tMtSq = (tE * tE - tPz * tPz);
    const double tM = std::sqrt(tMtSq - tPtSq);
    const double tMt = std::sqrt(tMtSq);
    const double tPt = std::sqrt(tPtSq);
    
    
    // Boost to LCMS

    const double beta_LCMS = tPz / tE;
    const double gamma_LCMS = tE / tMt;

    const double fDKOut = (part1.px() * tPx + part1.py() * tPy) / tPt;
    const double fDKSide = (-part1.px() * tPy + part1.py() * tPx) / tPt;
    const double fDKLong = gamma_LCMS * (part1.pz() - beta_LCMS * e1);
    const double fDE = gamma_LCMS * (e1 - beta_LCMS * part1.pz());

    const double px1LCMS = fDKOut;
    const double py1LCMS = fDKSide;
    const double pz1LCMS = fDKLong;
    const double pE1LCMS = fDE;

    const double px2LCMS = (part2.px() * tPx + part2.py() * tPy) / tPt;
    const double py2LCMS = (part2.py() * tPx - part2.px() * tPy) / tPt;
    const double pz2LCMS = gamma_LCMS * (part2.pz() - beta_LCMS * e2);
    const double pE2LCMS = gamma_LCMS * (e2 - beta_LCMS * part2.pz());

    const double fDKOutLCMS = px1LCMS - px2LCMS;
    const double fDKSideLCMS = py1LCMS - py2LCMS;
    const double fDKLongLCMS = pz1LCMS - pz2LCMS;
    
    //根据是否全同粒子处理
    if (isiden) {
      vect.push_back(qinv);
      vect.push_back(fDKOutLCMS);
      vect.push_back(fDKSideLCMS);
      vect.push_back(fDKLongLCMS);
      vect.push_back(qlcms);
    } else {
      vect.push_back(kstar);
      vect.push_back(fDKOut);
      vect.push_back(fDKSide);
      vect.push_back(fDKLong);
    }
    return vect;
  }
```

