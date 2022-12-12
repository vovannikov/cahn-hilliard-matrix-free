# Simple Cahn-Hilliard problem demo

## 

1. Checkout repo
```
git clone https://github.com/vovannikov/cahn-hilliard-matrix-free.git
```

2. Create directory for build and get in there
```
mkdir cahn-hilliard-matrix-free-build
cd cahn-hilliard-matrix-free-build
```

3. Configure the application providing `DEAL_II_DIR` and `CMAKE_BUILD_TYPE` (Release or Debug)
```
cmake ../cahn-hilliard-matrix-free -DDEAL_II_DIR=/path/to/dealii
```

4. Run the example
```
./cahn_hilliard_impl
```

## We solve the following equation
We solve the following problem:

```math
\dfrac{\partial c}{\partial t} = \nabla \cdot \left[ M \nabla \dfrac{\delta F}{\delta c} \right]
```
with the free energy defined as
```math
F \left(c\right) = \int_\Omega \left[ f \left(c\right) + \dfrac{1}{2} \kappa_c |c|^2 \right] \text{d} \Omega
```
where
```math
f \left(c\right) = c^2\left( c - 1 \right)^2
```
and mobility $M$ is a constant scalar.

By applying the functional derivative
```math
\dfrac{\delta J}{\delta g(\mathbf{x})} = \dfrac{\partial L}{\partial g} - \nabla \dfrac{\partial L}{\partial \nabla g}
```
and splitting the CH equation into two to lower the order of derivatives one obtains:
```math
\begin{empheq}[left = \empheqlbrace\,]{align}
	&\dfrac{\partial c}{\partial t} = \nabla \cdot \left[ M \nabla \mu \right], \\
	&\mu = \dfrac{\partial f}{\partial c} - \kappa_c \nabla^2 c.
\end{empheq}
```
whose weak form is written as
```math
\begin{empheq}[left = \empheqlbrace\,]{align}
	&\left( \dfrac{\partial c}{\partial t}, \varphi_c \right) + \left( M \nabla \mu, \nabla \varphi_c \right) - \langle M \nabla \mu \cdot \mathbf{n}, \varphi_c \rangle = 0,  \\
	&\left( \mu, \varphi_\mu \right) - \left( \dfrac{\partial f}{\partial c}, \varphi_\mu \right) - \left( \kappa_c \nabla c, \nabla \varphi_\mu \right) + \langle \kappa_c \nabla c \cdot \mathbf{n}, \varphi_\mu \rangle = 0.
\end{empheq}
```
where $\varphi$ denotes test functions for the variable identified buy the corresponding subscript and notation $\left( \cdot \right)$ and $\langle \cdot \rangle$ was used for domain and boundary integrals, respectively. For the no-flux boundary conditions the boundary integral vanish.

The matrix-free cell integral of the **right-hand side** is given by the following code:
```cpp
for (auto cell = range.first; cell < range.second; ++cell)
  {
    phi_old.reinit(cell);
    phi.reinit(cell);

    phi.gather_evaluate(src, EvaluationFlags::values | EvaluationFlags::gradients);

    // get values from old solution
    phi_old.read_dof_values_plain(old_solution);
    phi_old.evaluate(EvaluationFlags::values);

    for (unsigned int q = 0; q < phi.n_q_points; ++q)
      {
        const auto val  = phi.get_value(q);
        const auto grad = phi.get_gradient(q);
        // get old value
        const auto old_val = phi_old.get_value(q);
        
        // variables indices: 0 - c, 1 - mu

        Tensor<1, 2, VectorizedArrayType> value_result;
        value_result[0] = (val[0] - old_val[0]) / dt;
        value_result[1] = -val[1] + df_dc(val[0]);

        Tensor<1, 2, Tensor<1, dim, VectorizedArrayType>> gradient_result;
        gradient_result[0] = M * grad[1];
        gradient_result[1] = kappa * grad[0];

        phi.submit_value(value_result, q);
        phi.submit_gradient(gradient_result, q);
      }
    phi.integrate_scatter(EvaluationFlags::values | EvaluationFlags::gradients, dst);
  }
```
where an implicit backwards Euler time integration strategy was employed.
