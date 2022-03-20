# 矩阵形式

## 符号约定

$$
\begin{aligned}
D_{l} &= e^{-\Delta\tau \mathbf{c}^{\dagger} K \mathbf{c}} e^{-\Delta \tau \mathbf{c}^{\dagger} V(l) \mathbf{c}}\\
B_{l} &=  e^{-\Delta\tau K} e^{-\Delta \tau V(l)}
\end{aligned}
$$

$$
\begin{aligned}
c(l) &= (D_{l-1}\cdots D_{1})^{-1} c D_{l-1}\cdots D_{1}\\
c^{\dagger}(l) &= (D_{l-1}\cdots D_{1})^{-1} c^{\dagger} D_{l-1}\cdots D_{1}
\end{aligned}
$$

## 产生湮灭算符和$D$算符的对易关系

考虑算符 $c_{x}$和 $e^{\mathbf{c}^{\dagger} A \mathbf{c}}$ ,其中$A$为任意方阵。先考虑$c_{x}$与$\mathbf{c}^{\dagger} A \mathbf{c}$
$$
\begin{aligned}
c_{x} \mathbf{c}^{\dagger} A \mathbf{c} &= c_{x} c^{\dagger}_{y} A_{yz} c_{z} \\
&= \{c_{x}, c_{y}^{\dagger}\} A_{yz}c_{z} - c^{\dagger}_{y}c_{x}A_{yz}c_{z}\\
&= \delta_{xy} A_{yz}c_{z} + c_{y}^{\dagger} A_{yz} c_{z} c_{x} \\
&= \mathbf{c}^{\dagger} A \mathbf{c} c_{x} + A_{xy}c_{y}\\
&= \left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A \right)_{xy}c_{y} \\
%&= \left[\left(\mathbf{c}^{\dagger} A \mathbf{c} I + A \right) \mathbf{c} \right]_{x}
\end{aligned}
$$
其中$\mathbf{1}$是与$A$维度相同的单位矩阵。故
$$
c_{x} \left( \mathbf{c}^{\dagger} A \mathbf{c} \right)^{n} = \left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A \right)_{xy_{1}} \cdots\left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A \right)_{y_{n-1}y_{n}} c_{y_{n}} = \left[\left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A \right)^{n}\right]_{xy} c_{y}
$$
故
$$
\begin{aligned}
c_{x} e^{\mathbf{c}^{\dagger} A \mathbf{c}} &= c_{x} \sum_{n = 0}^{\infty} \frac{1}{n!} \left(\mathbf{c}^{\dagger} A \mathbf{c}\right)^{n} \\
&= \sum_{n = 0}^{\infty} \frac{1}{n!} \left[\left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A \right)^{n}\right]_{xy} c_{y} \\
&= \left[ e^{\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} + A} \right]_{xy} c_{y} \\
&= \left[e^{\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1}}e^{A} \right]_{xy} c_{y} \\
&= e^{\mathbf{c}^{\dagger} A \mathbf{c}} \left[e^{A}\right]_{xy} c_{y}
\end{aligned}
$$
即对于任意方阵，有
$$
c_{x} e^{\mathbf{c}^{\dagger} A \mathbf{c}} = e^{\mathbf{c}^{\dagger} A \mathbf{c}} \left[e^{A}\right]_{xy} c_{y}
$$
所以我们有
$$
c_{x} D_{l} = \left[B_{l}\right]_{xy} D_{l} c_{y}
$$
同理：
$$
\begin{aligned}
c^{\dagger}_{x} \mathbf{c}^{\dagger} A \mathbf{c} &= c^{\dagger}_{x} c^{\dagger}_{y} A_{yz} c_{z} \\
&= -c^{\dagger}_{y} c^{\dagger}_{x} A_{yz} c_{z} \\
&= -c^{\dagger}_{y} A_{yz} \{c^{\dagger}_{x}, c_{z}\} + c^{\dagger}_{y} A_{yz}c_{z} c^{\dagger}_{x}\\
&= -c^{\dagger}_{y} A_{yx} + \mathbf{c}^{\dagger} A \mathbf{c} c^{\dagger}_{x}\\
&= \left(\mathbf{c}^{\dagger} A \mathbf{c}\mathbf{1} - A\right)_{yx} c^{\dagger}_{y}
\end{aligned}
$$
且
$$
c^{\dagger}_{x} \left(\mathbf{c}^{\dagger} A \mathbf{c} \right)^{n} = \left(\mathbf{c}^{\dagger} A \mathbf{c}\mathbf{1} - A\right)_{y_{1}x} \cdots \left(\mathbf{c}^{\dagger} A \mathbf{c}\mathbf{1} - A\right)_{y_{n}y_{n-1}} c^{\dagger}_{y_{n}} = \left[\left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} - A \right)^{n}\right]_{yx} c^{\dagger}_{y}
$$
故
$$
\begin{aligned}
c^{\dagger}_{x} e^{\mathbf{c}^{\dagger} A \mathbf{c}} &= c^{\dagger}_{x} \sum_{n = 0}^{\infty} \frac{1}{n!} \left(\mathbf{c}^{\dagger} A \mathbf{c}\right)^{n} \\
&= \sum_{n = 0}^{\infty} \frac{1}{n!} \left[\left(\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} - A \right)^{n}\right]_{yx} c^{\dagger}_{y} \\
&= \left[ e^{\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1} - A} \right]_{yx} c^{\dagger}_{y} \\
&= \left[e^{\mathbf{c}^{\dagger} A \mathbf{c} \mathbf{1}}e^{-A} \right]_{yx} c_{y} \\
&= e^{\mathbf{c}^{\dagger} A \mathbf{c}} \left[e^{-A}\right]_{yx} c_{y}
\end{aligned}
$$
综上可得
$$
\begin{aligned}
c_{x} D_{l} &= \left[B_{l}\right]_{xy} D_{l} c_{y}\\
D_{l} c_{x} &= \left[B_{l}^{-1}\right]_{xy} c_{y} D_{l}\\
c^{\dagger}_{x} D_{l} &= D_{l} c^{\dagger}_{y} \left[B^{-1}_{l}\right]_{yx} \\
D_{l} c^{\dagger}_{x} &= c^{\dagger}_{y} D_{l} \left[B_{l}\right]_{yx}
\end{aligned}
$$

## 格林函数

定义
$$
\langle c_{x}(l) c^{\dagger}_{y}(l)\rangle = \frac{\mathrm{Tr} \left[D_{M}\cdots D_{l}c_{x} c^{\dagger}_{y} D_{l - 1} \cdots D_{1}\right]}{\mathrm{Tr} \left[D_{M} \cdots D_{1}\right]}
$$
其中
$$
\begin{aligned}
&\mathrm{Tr} \left[D_{M}\cdots D_{l}c_{x} c^{\dagger}_{y} D_{l - 1} \cdots D_{1}\right] \\
=& \mathrm{Tr} \left[D_{M}\cdots D_{l}\left(\delta_{xy} - c^{\dagger}_{y}c_{x}\right) D_{l - 1} \cdots D_{1}\right] \\
=& \mathrm{Tr} \left[D_{M} \cdots D_{1}\right] \delta_{xy} - \left[B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l}\right]_{xz} \mathrm{Tr} \left[D_{M}\cdots D_{l}c_{z} c^{\dagger}_{y} D_{l - 1} \cdots D_{1}\right]
\end{aligned}
$$
故
$$
\begin{aligned}
\langle c_{x}(l) c^{\dagger}_{y}(l)\rangle = \delta_{xy} - \left[B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l}\right]_{xz} \langle c_{z}(l) c^{\dagger}_{y}(l)\rangle
\end{aligned}
$$
即
$$
\left[\mathbf{1} + B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l}\right]_{xz} \langle c_{z}(l) c^{\dagger}_{y}(l)\rangle = \delta_{xy}
$$
故
$$
\langle c_{x}(l) c^{\dagger}_{y}(l)\rangle =\left[ \left(\mathbf{1} + B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l}\right)^{-1} \right]_{xy}
$$
定义格林函数
$$
G_{xy}(l_{1}, l_{2}) = \left\{ 
\begin{aligned}
&\langle c_{x}(l_{1})c_{y}^{\dagger}(l_{2}) \rangle & \mathrm{if} \quad l_{1} \ge l_{2}\\
-&\langle c_{y}^{\dagger}(l_{2}) c_{x}(l_{1})\rangle & \mathrm{if} \quad l_{1} \lt l_{2}
\end{aligned}
\right.
$$

$$
\tilde{G}_{xy}(l_{1}, l_{2}) = \left\{ 
\begin{aligned}
&\langle c_{x}(l_{1})c_{y}^{\dagger}(l_{2}) \rangle & \mathrm{if} \quad l_{1} \gt l_{2}\\
-&\langle c_{y}^{\dagger}(l_{2}) c_{x}(l_{1})\rangle & \mathrm{if} \quad l_{1} \le l_{2}
\end{aligned}
\right.
$$

接下来我们把$G$当成矩阵，下标$x$，$y$当成矩阵的指标来看待。$G$ 与 $\tilde{G}$的关系为
$$
\tilde{G}_{xy}(l_{1}, l_{2}) = \left\{ 
\begin{aligned}
&G     & \mathrm{if} \quad l_{1} \ne l_{2}\\
&G - \mathbf{1} & \mathrm{if} \quad l_{1} =  l_{2}
\end{aligned}
\right.
$$

$$
G\left(l, l\right) = \left[\mathbf{1} + B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l}\right]^{-1}
$$
且 
$$
G\left(l + 1, l + 1\right) = \left[\mathbf{1} + B_{l} \cdots B_{1} B_{M} \cdots B_{l + 1}\right]^{-1} = B_{l} G\left(l, l\right) B_{l}^{-1}
$$
若$l_{1} \gt l_{2}$
$$
G\left(l_{1},l_{2}\right) = B_{l_{1} - 1} \cdots B_{l_{2}} G\left(l_{2},l_{2}\right) = G\left(l_{1}, l_{1}\right) B_{l_{1} - 1} \cdots B_{l_{2}}
$$
若$l_{1} < l_{2}$
$$
G\left(l_{1},l_{2}\right) = B_{l_{1}}^{-1} \cdots B_{l_{2} - 1}^{-1} \tilde{G}\left(l_{2}, l_{2}\right) =  \tilde{G} \left( l_{1}, l_{1}\right) B^{-1}_{l_{1}} \cdots B^{-1}_{l_{2 - 1}}
$$
即
$$
G\left(l_{1},l_{2}\right) = 
\left\{
\begin{aligned}
& B_{l_{1} - 1} \cdots B_{l_{2}} G\left(l_{2},l_{2}\right) = G\left(l_{1}, l_{1}\right) B_{l_{1} - 1} \cdots B_{l_{2}} & \mathrm{if} \quad l_{1} \ge l_{2} \\
& B_{l_{1}}^{-1} \cdots B_{l_{2} - 1}^{-1} \tilde{G}\left(l_{2}, l_{2}\right) =  \tilde{G} \left( l_{1}, l_{1}\right) B^{-1}_{l_{1}} \cdots B^{-1}_{l_{2 - 1}} & \mathrm{if} \quad l_{1} \le l_{2}
\end{aligned}
\right.
$$

## 两点关联

$$
\langle c^{\dagger}_{x_{1}}(l) c_{y_{1}}(l) c^{\dagger}_{x_{2}}(l) c_{y_{2}}(l)\rangle = \frac{\mathrm{Tr} \left[D_{M}\cdots D_{l}c^{\dagger}_{x_{1}} c_{y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} D_{l - 1} \cdots D_{1}\right]}{\mathrm{Tr} \left[D_{M} \cdots D_{1}\right]}
$$

其中
$$
\begin{aligned}
& c^{\dagger}_{x_{1}} c_{y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}}\\
= &\delta_{x_{1}y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} - c_{y_{1}} c^{\dagger}_{x_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} \\
= & \delta_{x_{1}y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} + c_{y_{1}} c^{\dagger}_{x_{2}} c^{\dagger}_{x_{1}} c_{y_{2}} \\
= & \delta_{x_{1}y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} + c_{y_{1}} c^{\dagger}_{x_{2}} \delta_{x_{1}y_{2}} -  c_{y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} c^{\dagger}_{x_{1}}
\end{aligned}
$$
故
$$
\begin{aligned}
&\mathrm{Tr} \left[D_{M}\cdots D_{l} c^{\dagger}_{x_{1}} c_{y_{1}} c^{\dagger}_{x_{2}} c_{y_{2}} D_{l - 1} \cdots D_{1}\right]\\
= & \delta_{x_{1}y_{1}} \mathrm{Tr} \left[D_{M}\cdots D_{l} c^{\dagger}_{x_{2}} c_{y_{2}} D_{l - 1} \cdots D_{1}\right]
+ \delta_{x_{1}y_{2}} \mathrm{Tr} \left[D_{M}\cdots D_{l} c_{y_{1}} c^{\dagger}_{x_{2}}  D_{l - 1} \cdots D_{1}\right] \\
& - \mathrm{Tr} \left[D_{M}\cdots D_{l}c_{z} c^{\dagger}_{y_{1}} c_{x_{2}} c^{\dagger}_{y_{2}} D_{l - 1} \cdots D_{1}\right] \left[B^{-1}_{l} \cdots B^{-1}_{M} B^{-1}_{1}\cdots B^{-1}_{l-1}\right]_{zx_{1}}
\end{aligned}
$$


故
$$
\langle c^{\dagger}_{z}(l) c_{y_{1}}(l) c^{\dagger}_{x_{2}}(l) c_{y_{2}}(l)\rangle \left[\mathbf{1} + \left(B_{l - 1} \cdots B_{1} B_{M} \cdots B_{l} \right)^{-1}\right]_{zx_{1}}= \delta_{x_{1}y_{1}} \langle c^{\dagger}_{x_{2}}(l) c_{y_{2}} (l) \rangle + \delta_{x_{1}y_{2}} \langle c_{y_{1}}(l) c^{\dagger}_{x_{2}}(l) \rangle
$$
即
$$
\begin{aligned}
\langle c^{\dagger}_{x_{1}}(l) c_{y_{1}}(l) c^{\dagger}_{x_{2}}(l) c_{y_{2}}(l) \rangle 
&= \langle c^{\dagger}_{x_{1}}(l) c_{y_{1}}(l) \rangle \langle c^{\dagger}_{x_{2}}(l) c_{y_{2}}(l)\rangle
+ \langle c^{\dagger}_{x_{1}}(l) c_{y_{2}}(l)\rangle \langle c_{y_{1}}(l) c^{\dagger}_{x_{2}}(l) \rangle\\
&= \tilde{G}_{y_{1}x_{1}}\left(l,l\right) \tilde{G}_{y_{2}x_{2}}\left(l,l\right) - G_{y_{1}x_{2}} \left(l,l\right) \tilde{G}_{y_{2}x_{1}}\left(l,l\right)
\end{aligned}
$$
考察关联函数$\langle c^{\dagger}_{x}(l_{1})c_{x}(l_{1})c^{\dagger}_{y}(l_{2}) c_{y}(l_{2}) \rangle$,其中$l_{1} > l_{2}$
$$
\begin{aligned}
&\langle c^{\dagger}_{x}(l_{1})c_{x}(l_{1})c^{\dagger}_{y}(l_{2}) c_{y}(l_{2}) \rangle \\
=& \left[B_{l_{1} - 1} \cdots B_{l_{2}}\right]_{xz_{1}} \langle c^{\dagger}_{z_{2}}(l_{2})c_{z_{1}}(l_{2})c^{\dagger}_{y}(l_{2}) c_{y}(l_{2}) \rangle  \left[\left(B_{l_{1} - 1} \cdots B_{l_{2}}\right)^{-1}\right]_{z_{2}x} \\
=& \left[B_{l_{1} - 1} \cdots B_{l_{2}}\right]_{xz_{1}} \left( \langle c^{\dagger}_{z_{2}}\left(l_{2}\right) c_{z_{1}}(l_{2}) \rangle \langle c^{\dagger}_{y}\left(l_{2}\right) c_{y}\left(l_{2}\right)\rangle  + \langle c^{\dagger}_{z_{2}} \left(l_{2}\right) c_{y}\left(l_{2}\right) \rangle \langle c_{z_{1}}\left(l_{2}\right) c^{\dagger}_{y}\left(l_{2})\right)\rangle \right) \left[\left(B_{l_{1} - 1} \cdots B_{l_{2}}\right)^{-1}\right]_{z_{2}x} \\
=& \langle c^{\dagger}_{x}\left(l_{1}\right) c_{x}(l_{1}) \rangle \langle c^{\dagger}_{y}\left(l_{2}\right) c_{y}\left(l_{2}\right)\rangle + \langle c^{\dagger}_{x} \left(l_{1}\right) c_{y}\left(l_{2}\right) \rangle \langle c_{x}\left(l_{1}\right) c^{\dagger}_{y}\left(l_{2})\right)\rangle \\
=& \tilde{G}_{xx}\left(l_{1},l_{1}\right) \tilde{G}_{yy}\left(l_{2},l_{2}\right) - \tilde{G}_{yx}\left(l_{2},l_{1}\right) G_{xy} \left(l_{1},l_{2}\right)
\end{aligned}
$$
若 $l_{1} < l_{2}$
$$
\langle c^{\dagger}_{y}(l_{2}) c_{y}(l_{2}) c^{\dagger}_{x}(l_{1})c_{x}(l_{1}) \rangle =\tilde{G}_{yy}\left(l_{2},l_{2}\right)\tilde{G}_{xx}\left(l_{1},l_{1}\right) - \tilde{G}_{xy}\left(l_{1},l_{2}\right) G_{yx} \left(l_{2},l_{1}\right)
$$
注意到当$l_{1} \ne l_{2}$ 时， $\tilde{G} = G$ 故
$$
\langle T c^{\dagger}_{x}(l_{1})c_{x}(l_{1})c^{\dagger}_{y}(l_{2}) c_{y}(l_{2}) \rangle = \tilde{G}_{xx}\left(l_{1},l_{1}\right) \tilde{G}_{yy}\left(l_{2},l_{2}\right) - \tilde{G}_{yx}\left(l_{2},l_{1}\right) G_{xy} \left(l_{1},l_{2}\right)
$$
对于任意$l_{1}$和$l_{2}$都成立。