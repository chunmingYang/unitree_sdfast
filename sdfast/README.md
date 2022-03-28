# How to use this template?
1. win environment --- libraries generation

        input: sdfast.exe + sdfast.key + a1.sd
        tmd: sdfast -lc -ge -p a1 a1.sd a1
        output: a1_i + a1_d.c + a1_s.c + a1lib.c
2. linux environment

        input: a1_main.c + a1_d.c + a1_s.c + a1lib.c + run
        tmd: ./run
        output: a1 --- executable
