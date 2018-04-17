title: Model

# Model

![ROHSA](|media|/LogoMakr_0dTJ9B.png)
{: style="text-align: center" }

The two dimensional vector expressing the line of sight is denoted by \((l, b)\).

Lets start with a non linear least square formulation of the optimization problem. 
We have \(p\) empirical pairs \((x_u^{l,b}, y_u^{l,b})\) and we want to find the parameters \(\mathbf{\beta}\) of the model curve \(G(\mathbf{\beta})\)
so that the sum of the squares of deviations \(f(\mathbf{\beta})\) is minimized. 

\begin{align}
  \hat{\mathbf{\beta}} = \underset{\beta}{argmin} \, f(\mathbf{\beta})
\end{align}

Let's define the model G as a sum of Gaussian characterized by three parameters, i.e. amplitude \(A : \, q = 3k+1\), 
position \(\mu : \, q = 3k+2\) and dispersion \(\sigma : \, q = 3k+3\) with \(k \in \left\{ 0, ..., N \right\}\)

\begin{align}
  G(x_u^{l,b}, \mathbf{\beta_q^{l,b}}) = \sum_{k=1}^{N} \left ( \beta_{3k}^{l,b} \, exp \left \{  \frac{- \, (x_u^{l,b} \, - \,
    \beta_{3k+1}^{l,b})^2}{2 (\beta_{3k+2}^{l,b})^2}
    \right\} \right)
\end{align} 

We define the residual function \(F\) of each pairs \((x_u^{l,b}, y_u^{l,b})\) :

\begin{align}
  F(x_u^{l,b}, y_u^{l,b}, \mathbf{\beta^{l,b}}) = G(x_u^{l,b}, \mathbf{\beta^{l,b}}) \, - \, y_u^{l,b}
\end{align} 

For a non linear least square formulation, the objective function \(f\) is :

\begin{align}
  f(\mathbf{\beta}) = \frac{1}{2} \, \sum_{l=1}^L \, \sum_{b=1}^B \left [ \sum_{u=1}^{U} \, F^2(x_u^{l,b}, y_u^{l,b}, \mathbf{\beta^{l,b}}) \right ]
  \, + \, \underbrace{\sum_{q=1}^{3n} \, \lambda_q \, R_q(\mathbf{\beta_q})}_{regularization}
\end{align}

Note that we add a regularization term in order to minimize \(f(\mathbf{\beta_q})\) considering spatial constraints. 

If we regularize over 4 neighbors \(R_q(\mathbf{\beta})\) is equal to,

\begin{align}
  R_q(\mathbf{\beta_q}) = \frac{1}{2} \, \sum_{l=1}^L \, \sum_{b=1}^B \left[ Y_q^{l,b} \right]^2
\end{align}

with 
\begin{align}
	Y_q^{l,b} = \frac{1}{4} \, \left ( 4 \, \beta_q^{l,b} \, - \,
    \beta_q^{l+1,b} \, - \, \beta_q^{l,b+1} \, - \, \beta_q^{l,b-1} \, - \, \beta_q^{l-1,b} \right ) 
\end{align}

In order to consider the boundary of the cube, we use a mirror conditions.

Furthermore, it is convenient to express \(R_q(\mathbf{\beta_q})\) as a convolution product :

\begin{align}
  \mathbf{Y_q} = \mathbf{K} \, * \, \mathbf{\beta_q}
\end{align}

\noindent
with the kernel \(\textbf{K}\),

\begin{align}
  \mathbf{K} =
  \left[
  \begin{array}{ccc}
    0 & -1 & 0 \\
    -1 & 4 & -1 \\
    0 & -1 & 0 \\
  \end{array}
  \right] \, / \, 4
\end{align}

The gradient of the objective function that we need to use the optimization algorithm of our choice is then : 

\begin{align}
  \nabla f(\mathbf{\beta}) &= \frac{\partial f(\mathbf{\beta})}{\partial \beta_q^{l,b}} 
  &= \sum_{u=1}^U \, \frac{\partial F_u^{l,b}}{\partial \beta_q^{l,b}} \, \times \,
  F_u^{l,b} \, + \, \lambda_q \, \frac{\partial}{\partial \beta_q^{l,b}} \left[ \frac{1}{2} \, \sum_{l=1}^L \, \sum_{b=1}^B \left[ Y_q^{l,b} \right]^2
    \right]
\end{align}


Finally, the gradient of the objective function can be written as follows : 

\begin{align}
  \nabla f(\mathbf{\beta}) = \frac{\partial f(\mathbf{\beta})}{\partial \beta_q^{l,b}} = \sum_{u=1}^U \, \frac{\partial F_u^{l,b}}{\partial \beta_q^{l,b}} \, \times \,
  F_u^{l,b} \, + \, \lambda_q \, Z_q^{l,b}
\end{align}

with, 
\begin{align}
	\mathbf{Z_q} = \mathbf{K} \, * \, \mathbf{Y_q}
\end{align}