# OpenMP-Multithread-Framework
3 sequential algorithm (simplex, game-of-life and heatmap) were sped up using different number of threads.

Full specification of those algorithms and constraints about their speeding up can be seen [here](https://github.com/mdodovic/OpenMP-Multithread-Framework/blob/main/problems_description.pdf).

## [task1](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

Sequential code is modified and then both basic and modified codes are sped up using **manual scheduling** of all _unit-of-work_.

## [task2](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

Sequential code is modified and then both basic and modified codes are sped up using the **working sharing directive _for_**. This working sharing directive can divide the total amount of work using different techniques (_static_, _dynamic_, _guided_) with variable _unit-of-work_ sizes.

## [task3](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

Sequential code is sped up using the **working sharing directive _for_** with the _static_ technique of dividing the total amount of work into _unit-of-work_ with size 1.

## [task4](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

## [task5](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)
