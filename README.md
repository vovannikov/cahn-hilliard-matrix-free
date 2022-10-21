# Simple Cahn-Hilliard problem deom

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
