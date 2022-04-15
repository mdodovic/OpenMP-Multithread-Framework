import matplotlib.pyplot as plt


def main():

    chunksize_for_scheduling = [1, 2, 5, 10, 30]

    simplex_modified_omp_static = [3.696957407, 3.530285261, 3.646845894, 2.894587966, 1.316416783]
    simplex_modified_omp_dynamic = [3.685526858, 3.586754597, 3.710329146, 2.896641349, 1.300100706]
    simplex_modified_omp_guided = [3.806597922, 3.641728315, 3.740723792, 2.90583124, 1.298719004, ]


    plt.plot(chunksize_for_scheduling, simplex_modified_omp_static)
    plt.plot(chunksize_for_scheduling, simplex_modified_omp_dynamic)
    plt.plot(chunksize_for_scheduling, simplex_modified_omp_guided)
    plt.legend(["static", "dynamic", "guided"], loc="best")
    plt.savefig("simplex_modified.png", dpi = 90)
    plt.show()


if __name__ == "__main__":
    main()