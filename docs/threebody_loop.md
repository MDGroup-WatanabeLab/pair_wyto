# Threebody loop

## Setting formulas and characters

$$ \Lambda_n(i,j,k)=\lambda_n(z_i) \Lambda_n^{\prime}(r_{ij},r_{ik}) $$

$$ j,k,l　\in　\text{neighbors of } i $$

## Loop
$\underset{n}{\Lambda}
  \underset{i}{\Lambda}
    \underset{j\neq i}{\Lambda}
      \underset{k \neq i}{\underset{k > j}{\Lambda}} \lbrace$

$$
\begin{eqnarray}
V \ &+=& \varepsilon \lambda_n(z_i) \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik}) \\
\boldsymbol{F_i} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{ij}}
  \frac{\boldsymbol{r_i}-\boldsymbol{r_j}}{r_{ij}} \\
\boldsymbol{F_j} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{ij}}
  \frac{\boldsymbol{r_j}-\boldsymbol{r_i}}{r_{ij}} \\
\boldsymbol{F_i} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{ik}}
  \frac{\boldsymbol{r_i}-\boldsymbol{r_k}}{r_{ik}} \\
\boldsymbol{F_k} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{ik}}
  \frac{\boldsymbol{r_k}-\boldsymbol{r_i}}{r_{ik}} \\
\boldsymbol{F_j} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{jk}}
  \frac{\boldsymbol{r_j}-\boldsymbol{r_k}}{r_{jk}} \\
\boldsymbol{F_k} \ &+=& -\varepsilon \lambda_n(z_i)
  \frac{\partial \Lambda_n^{\prime}(r_{ij},r_{ik}) \Theta_n(\theta_{jik})}{\partial r_{jk}}
  \frac{\boldsymbol{r_k}-\boldsymbol{r_j}}{r_{jk}} \\
\end{eqnarray}
$$

$\qquad\qquad\quad \underset{l \neq i}{\Lambda} \lbrace$

$$
\begin{eqnarray}
\boldsymbol{F_i} \ &+=& -\varepsilon
  \Lambda_n^{\prime}(r_{ij},r_{ik})\Theta_n(\theta_{jik})
  \frac{\partial\lambda_n(z_i)}{\partial z_i}
  \frac{\partial f_c(r_{il})}{\partial r_{il}}
  \frac{\boldsymbol{r_i}-\boldsymbol{r_l}}{r_{il}} \\
\boldsymbol{F_l} \ &+=& -\varepsilon
  \Lambda_n^{\prime}(r_{ij},r_{ik})\Theta_n(\theta_{jik})
  \frac{\partial\lambda_n(z_i)}{\partial z_i}
  \frac{\partial f_c(r_{il})}{\partial r_{il}}
  \frac{\boldsymbol{r_l}-\boldsymbol{r_i}}{r_{il}}
\end{eqnarray}
$$

$\qquad\qquad\quad \rbrace$


$\rbrace$
