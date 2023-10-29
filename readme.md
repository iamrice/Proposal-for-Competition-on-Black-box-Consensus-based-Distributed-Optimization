
运行：
```
cd ./user_implement

mpic++ -std=c++11 user_implement.cpp ../framework/framework.cpp ../framework/main.cpp ../Benchmarks/Benchmarks.cpp -o test

mpirun -np 20 ./test ring-20n-1

```

参赛者只需修改./user_implement/user_implement.cpp中的agent_function函数。
上交代码时，只需打包./user_implement文件夹。