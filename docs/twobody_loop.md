# Twobody loop

## Setting formulas and characters

$$f_2(r_{ij}) = A_{ij}(B_{ij}r_{ij}^{-p_{ij}}-r_{ij}^{-q_{ij}})\exp[(r_{ij}-a_{ij})^{-1}]$$

$$j,k \in \text{neighbors of } i $$

## Loop
### Regular force loop

$\underset{i}{\Lambda}\underset{j > i}{\Lambda}\lbrace$

$$
\begin{eqnarray}
V \ &+=& \varepsilon g_1(z_i)g_2(z_j)f_2(r_{ij}) \\
\boldsymbol{F_i} \ &+=&
    -\varepsilon g_1(z_i)g_2(z_j) 
    \frac{\partial f_2(r_{ij})}{\partial r_{ij}}
    \frac{\boldsymbol{r_i}-\boldsymbol{r_j}}{r_{ij}} \\
\boldsymbol{F_j} \ &+=&
    -\varepsilon g_1(z_i)g_2(z_j) 
    \frac{\partial f_2(r_{ij})}{\partial r_{ij}} 
    \frac{\boldsymbol{r_j}-\boldsymbol{r_i}}{r_{ij}}
\end{eqnarray}
$$

$\rbrace$

### External force loop

$\underset{i}{\Lambda}\underset{j\neq i}{\Lambda}\underset{k\neq i}{\Lambda}\lbrace$

$$
\begin{eqnarray}
\boldsymbol{F_i} \ &+=&
    -\varepsilon g_2(z_j)f_2(r_{ij})
    \frac{\partial g_1(z_i)}{\partial z_i}
    \frac{\partial f_c(r_{ik})}{\partial r_{ik}}
    \frac{\boldsymbol{r_i} - \boldsymbol{r_k}}{r_{ik}} \\
\boldsymbol{F_k} \ &+=&
    -\varepsilon g_2(z_j)f_2(r_{ij})
    \frac{\partial g_1(z_i)}{\partial z_i}
    \frac{\partial f_c(r_{ik})}{\partial r_{ik}}
    \frac{\boldsymbol{r_k} - \boldsymbol{r_i}}{r_{ik}}
\end{eqnarray}
$$

$\rbrace$
