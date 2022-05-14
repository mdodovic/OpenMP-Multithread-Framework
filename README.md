# OpenMP-Multithread-Framework
3 sequential algorithm (simplex, game-of-life and heatmap) were sped up using different number of threads.

Full specification of those algorithms and constraints about their speeding up can be seen [here](https://github.com/mdodovic/OpenMP-Multithread-Framework/blob/main/problems_description.pdf).

## [task1](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

Sequential code is modified and then both basic and modified codes are sped up using **manual scheduling** of all _unit-of-work_.

## [task2](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

Sequential code is modified and then both basic and modified codes are sped up using **working sharing directive _for_**. This working sharing directive can divide total amount using different techniques (_static_, _dynamic_, _guided_) with variable _unit-of-work_ size.

## [task3](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

## [task4](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)

## [task5](https://github.com/mdodovic/OpenMP-Multithread-Framework/tree/main/task1_simplex)
